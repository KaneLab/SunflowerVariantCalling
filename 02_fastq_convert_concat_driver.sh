#!/bin/bash
# Driver for 02_fastq_convert_concat.pbs
# Converts all SRA files for one genotype to gzipped FASTQ, then concatenates them.
#
# Usage (called by 02_fastq_convert_concat.pbs — not normally invoked directly):
#   02_fastq_convert_concat_driver.sh GENOTYPE MANIFEST SRA_DIR FASTQ_DIR
#
# Respects these environment variables (set by the PBS script):
#   DELETE_SRA       — "true" to delete .sra after conversion (default: false)
#   MIN_FREE_GB      — minimum free GB required in FASTQ_DIR before each fasterq-dump (default: 100)
#   INFLATION_FACTOR — .sra → uncompressed FASTQ size multiplier for space estimation (default: 8)
#   NCPUS            — threads for fasterq-dump and pigz (default: 8)
#   SCRATCH          — PBS scratch directory for temporary uncompressed FASTQs
#
# Concatenation strategy:
#   Multiple .fastq.gz files per genotype are joined using `cat`, which produces a valid
#   multi-stream gzip file. This avoids any decompression/recompression overhead.
#   fastp, bwa-mem2, zcat, and bcftools all handle multi-stream gzip correctly.
#   Validation uses file-size checksums (no decompression) + gzip -t integrity check.

set -euo pipefail

GENOTYPE="$1"
MANIFEST="$2"
SRA_DIR="$3"
FASTQ_DIR="$4"

DELETE_SRA="${DELETE_SRA:-false}"
MIN_FREE_GB="${MIN_FREE_GB:-100}"
INFLATION_FACTOR="${INFLATION_FACTOR:-8}"
NCPUS="${NCPUS:-8}"
SCRATCH="${SCRATCH:-/tmp}"

DONE_DIR="${FASTQ_DIR}/.done"
DONE_MARKER="${DONE_DIR}/${GENOTYPE}.done"
LOG="${FASTQ_DIR}/${GENOTYPE}.convert_concat.log"

mkdir -p "$FASTQ_DIR" "$DONE_DIR"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

###############################################################################
# Skip if a previous run already completed this genotype
###############################################################################
if [ -f "$DONE_MARKER" ]; then
    log "$GENOTYPE: already complete (found $DONE_MARKER). Skipping."
    exit 0
fi

log "=== Starting $GENOTYPE ==="

###############################################################################
# Look up SRRs for this genotype in the MANIFEST
###############################################################################
mapfile -t SRRS < <(awk -v g="$GENOTYPE" '$1==g {print $2}' "$MANIFEST")
if [ "${#SRRS[@]}" -eq 0 ]; then
    log "ERROR: No SRRs found for genotype '$GENOTYPE' in $MANIFEST"
    exit 1
fi
log "Found ${#SRRS[@]} SRR(s): ${SRRS[*]}"

###############################################################################
# check_storage DIR NEEDED_GB
#   Exits non-zero if DIR has less than NEEDED_GB free.
###############################################################################
check_storage() {
    local dir="$1"
    local needed_gb="$2"
    local avail_kb avail_gb
    avail_kb=$(df -k "$dir" | awk 'NR==2 {print $4}')
    avail_gb=$(( avail_kb / 1048576 ))
    if [ "$avail_gb" -lt "$needed_gb" ]; then
        log "ERROR: Insufficient space in $dir."
        log "  Need ~${needed_gb}G, have ${avail_gb}G."
        log "  Free space or reduce parallelism (smaller -J range), then resubmit."
        return 1
    fi
    log "Storage OK in $dir: ${avail_gb}G available (need ~${needed_gb}G)"
    return 0
}

###############################################################################
# Convert each SRR to per-SRR .fastq.gz files in FASTQ_DIR
###############################################################################
shopt -s extglob nullglob
PER_SRR_R1=()
PER_SRR_R2=()

for SRR in "${SRRS[@]}"; do

    R1="${FASTQ_DIR}/${GENOTYPE}_${SRR}_1.fastq.gz"
    R2="${FASTQ_DIR}/${GENOTYPE}_${SRR}_2.fastq.gz"

    # Skip fasterq-dump if this SRR was already converted in a prior (interrupted) run
    if [ -f "$R1" ] && [ -f "$R2" ]; then
        log "$SRR: already converted, skipping fasterq-dump"
        PER_SRR_R1+=("$R1")
        PER_SRR_R2+=("$R2")
        continue
    fi

    # Locate the .sra or .sralite file
    sra_candidates=( "${SRA_DIR}/${SRR}/${SRR}".@(sra|sralite) )
    if [ "${#sra_candidates[@]}" -eq 0 ]; then
        log "ERROR: SRA file not found for $SRR in ${SRA_DIR}/${SRR}/"
        exit 1
    fi
    SRA_FILE="${sra_candidates[0]}"

    # Estimate space needed: uncompressed FASTQs go to $SCRATCH
    sra_bytes=$(stat -c%s "$SRA_FILE")
    scratch_needed_gb=$(( (sra_bytes * INFLATION_FACTOR / 1000000000) + 5 ))
    fastq_needed_gb=$(( (sra_bytes * 2 / 1000000000) + MIN_FREE_GB ))

    check_storage "$SCRATCH"    "$scratch_needed_gb"
    check_storage "$FASTQ_DIR"  "$fastq_needed_gb"

    log "$SRR: running fasterq-dump ($(( sra_bytes / 1000000 ))MB .sra → ~${scratch_needed_gb}G uncompressed in \$SCRATCH)"

    SCRATCH_TMP="${SCRATCH}/${GENOTYPE}_${SRR}_tmp"
    mkdir -p "$SCRATCH_TMP"

    # Dump to $SCRATCH to avoid inflating project quota with uncompressed FASTQs
    fasterq-dump \
        --temp    "$SCRATCH_TMP" \
        --split-3 \
        -e        "$NCPUS" \
        -O        "$SCRATCH_TMP" \
        "$SRA_FILE"

    # Compress each mate directly from $SCRATCH to $FASTQ_DIR; delete raw immediately
    for mate in 1 2; do
        RAW="${SCRATCH_TMP}/${SRR}_${mate}.fastq"
        DEST="${FASTQ_DIR}/${GENOTYPE}_${SRR}_${mate}.fastq.gz"
        if [ ! -f "$RAW" ]; then
            # Some runs are single-end; R1 must exist, R2 is optional
            if [ "$mate" -eq 1 ]; then
                log "ERROR: fasterq-dump did not produce ${RAW} for $SRR"
                rm -rf "$SCRATCH_TMP"
                exit 1
            else
                log "WARNING: No R2 produced for $SRR — this may be single-end data."
                continue
            fi
        fi
        log "  Compressing mate${mate}: $(( $(stat -c%s "$RAW") / 1000000 ))MB -> $DEST"
        pigz -p "$NCPUS" -c "$RAW" > "$DEST"
        rm -f "$RAW"
    done

    # Remove any remaining files (unpaired reads, temp index files, etc.)
    rm -rf "$SCRATCH_TMP"

    # Optionally reclaim .sra storage immediately after successful conversion
    if [ "$DELETE_SRA" = "true" ]; then
        log "  Deleting $SRA_FILE (DELETE_SRA=true)"
        rm -f "$SRA_FILE"
        rmdir "${SRA_DIR}/${SRR}" 2>/dev/null || true
    fi

    PER_SRR_R1+=("$R1")
    [ -f "$R2" ] && PER_SRR_R2+=("$R2")
    log "$SRR: conversion complete"
done

###############################################################################
# Concatenate per-SRR .fastq.gz files
#
# `cat` joins gzip streams end-to-end, producing a valid multi-stream gzip file.
# No decompression or recompression occurs — this is a pure byte copy.
# Validation: output byte size must equal the sum of input sizes, then gzip -t.
###############################################################################
SRR_TAG=$(IFS=_; echo "${SRRS[*]}")
CONCAT_R1="${FASTQ_DIR}/${GENOTYPE}_${SRR_TAG}_1.fastq.gz"
CONCAT_R2="${FASTQ_DIR}/${GENOTYPE}_${SRR_TAG}_2.fastq.gz"

concat_and_validate() {
    local mate="$1"; shift
    local out="$1"; shift
    local inputs=("$@")

    if [ "${#inputs[@]}" -eq 0 ]; then
        log "WARNING: No R${mate} inputs to concatenate for $GENOTYPE"
        return 0
    fi

    if [ "${#inputs[@]}" -eq 1 ]; then
        # Single file — rename in place
        mv "${inputs[0]}" "$out"
        log "mate${mate}: single file renamed -> $out"
    else
        log "mate${mate}: concatenating ${#inputs[@]} files -> $out"
        cat "${inputs[@]}" > "${out}.tmp"
        mv "${out}.tmp" "$out"

        # Validate: output size must equal sum of input sizes
        expected=0
        for f in "${inputs[@]}"; do
            expected=$(( expected + $(stat -c%s "$f") ))
        done
        actual=$(stat -c%s "$out")
        if [ "$actual" -ne "$expected" ]; then
            log "ERROR: Concatenated file size mismatch for mate${mate} (expected ${expected}B, got ${actual}B)"
            exit 1
        fi
        log "  mate${mate}: size check passed (${actual}B)"

        # Remove per-SRR files now that concatenation is verified
        rm -f "${inputs[@]}"
    fi

    # Gzip integrity check
    if ! gzip -t "$out"; then
        log "ERROR: gzip integrity check failed for $out"
        exit 1
    fi
    log "  mate${mate}: gzip -t passed"
}

concat_and_validate 1 "$CONCAT_R1" "${PER_SRR_R1[@]}"
if [ "${#PER_SRR_R2[@]}" -gt 0 ]; then
    concat_and_validate 2 "$CONCAT_R2" "${PER_SRR_R2[@]}"
fi

###############################################################################
# Mark this genotype as complete
###############################################################################
touch "$DONE_MARKER"
log "=== Done: $GENOTYPE ==="
log "  R1: $CONCAT_R1"
[ -f "$CONCAT_R2" ] && log "  R2: $CONCAT_R2"
