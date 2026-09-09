#!/usr/bin/env bash
#
# 00_prep_local_fastq.sh
#
# Alternative entry point when paired-end FASTQs already exist on the HPC.
# Creates symlinks with consistent naming in FASTQ_DIR and writes GENLIST.txt —
# the same outputs expected from steps 00-02 — so steps 03-06 can run unchanged.
#
# Symlinks are used instead of copies so no additional disk space is consumed.
# The source files can follow any naming convention; only the symlink names need
# to satisfy the glob used by 03_trim_align_markdup.pbs:
#     {sample}_R1.fastq.gz  (or .fq.gz)
#     {sample}_R2.fastq.gz  (or .fq.gz)
#
# REQUIRED ARGS (positional):
#   1. MANIFEST  — tab-delimited: SAMPLE_NAME<TAB>R1_PATH<TAB>R2_PATH
#                  Lines beginning with '#' are skipped as comments.
#                  Blank lines are skipped.
#                  Absolute paths are recommended; relative paths are resolved
#                  from the current working directory.
#   2. FASTQ_DIR — directory where symlinks will be created
#   3. LIST_DIR  — directory where GENLIST.txt will be written
#
# Usage:
#   bash 00_prep_local_fastq.sh /path/to/manifest.tsv /path/to/fastq_dir /path/to/list_dir
#
# Outputs:
#   {FASTQ_DIR}/{sample}_R1.{ext}   — symlink → source R1 FASTQ
#   {FASTQ_DIR}/{sample}_R2.{ext}   — symlink → source R2 FASTQ
#   {LIST_DIR}/GENLIST.txt          — one sample name per line
#
# After this script completes, pass GENLIST.txt and FASTQ_DIR to step 03:
#   n=$(wc -l < /path/to/list_dir/GENLIST.txt)
#   qsub -J 1-$n \
#        -v "SAMPLE_LIST=/path/to/list_dir/GENLIST.txt,\
#            FQ_DIR=/path/to/fastq_dir,\
#            OUT_DIR=/path/to/bam" \
#        03_trim_align_markdup.pbs

set -euo pipefail

usage() {
    echo "Usage: bash $0 <MANIFEST> <FASTQ_DIR> <LIST_DIR>" >&2
    echo "" >&2
    echo "  MANIFEST  — tab-delimited: SAMPLE_NAME<TAB>R1_PATH<TAB>R2_PATH" >&2
    echo "              Lines beginning with '#' and blank lines are skipped." >&2
    echo "  FASTQ_DIR — directory where symlinks will be created" >&2
    echo "  LIST_DIR  — directory where GENLIST.txt will be written" >&2
    exit 1
}

[[ $# -ne 3 ]] && usage

manifest_raw="$1"
fastq_dir="$2"
list_dir="$3"

[[ ! -f "$manifest_raw" ]] && { echo "ERROR: manifest not found: $manifest_raw" >&2; exit 1; }
manifest=$(realpath "$manifest_raw")

mkdir -p "$fastq_dir" "$list_dir"

genlist="${list_dir}/GENLIST.txt"
> "$genlist"

errors=0
lineno=0
declare -A seen_samples

while IFS=$'\t' read -r sample r1 r2; do
    lineno=$(( lineno + 1 ))

    # Strip carriage returns (Windows line endings)
    sample="${sample//$'\r'/}"
    r1="${r1//$'\r'/}"
    r2="${r2//$'\r'/}"

    # Skip blank lines and comment lines
    [[ -z "$sample" || "$sample" == \#* ]] && continue

    # Validate all three fields are present
    if [[ -z "$r1" || -z "$r2" ]]; then
        echo "ERROR: line $lineno: expected 3 tab-delimited fields (SAMPLE_NAME, R1_PATH, R2_PATH)" >&2
        errors=$(( errors + 1 ))
        continue
    fi

    # Detect duplicate sample names
    if [[ -v seen_samples["$sample"] ]]; then
        echo "ERROR: line $lineno: duplicate sample name '$sample' (first seen at line ${seen_samples[$sample]})" >&2
        errors=$(( errors + 1 ))
        continue
    fi
    seen_samples["$sample"]=$lineno

    # Check that source files exist and resolve to absolute paths
    if [[ ! -f "$r1" ]]; then
        echo "ERROR: line $lineno ($sample): R1 not found: $r1" >&2
        errors=$(( errors + 1 ))
        continue
    fi
    r1=$(realpath "$r1")

    if [[ ! -f "$r2" ]]; then
        echo "ERROR: line $lineno ($sample): R2 not found: $r2" >&2
        errors=$(( errors + 1 ))
        continue
    fi
    r2=$(realpath "$r2")

    # Detect file extension from R1 to preserve it in the symlink name.
    # 03_trim_align_markdup.pbs matches: ${sample}_*1.f*q.gz and ${sample}_*2.f*q.gz
    # Both .fastq.gz and .fq.gz satisfy that glob.
    case "$r1" in
        *.fastq.gz) ext="fastq.gz" ;;
        *.fq.gz)    ext="fq.gz" ;;
        *.fastq|*.fq)
            echo "WARNING: line $lineno ($sample): source file appears uncompressed." >&2
            echo "         Downstream steps expect gzip-compressed input (.fastq.gz or .fq.gz)." >&2
            echo "         Compress with: pigz -p 8 $r1" >&2
            ext="${r1##*.}"
            ;;
        *)
            echo "WARNING: line $lineno ($sample): unrecognised extension for $(basename "$r1") — symlinking as-is." >&2
            ext="${r1##*.}"
            ;;
    esac

    sym_r1="${fastq_dir}/${sample}_R1.${ext}"
    sym_r2="${fastq_dir}/${sample}_R2.${ext}"

    # -s: symbolic, -f: overwrite any existing symlink at that path
    ln -sf "$r1" "$sym_r1"
    ln -sf "$r2" "$sym_r2"

    echo "$sample" >> "$genlist"
    printf "  %-30s %s / %s\n" "$sample" "$(basename "$r1")" "$(basename "$r2")"

done < "$manifest"

if [[ $errors -gt 0 ]]; then
    echo "" >&2
    echo "ERROR: $errors problem(s) found above. Fix the manifest and re-run." >&2
    rm -f "$genlist"
    exit 1
fi

n=$(wc -l < "$genlist")

echo ""
echo "Done: $n sample(s)"
echo "  GENLIST : $genlist"
echo "  Symlinks: $fastq_dir"
echo ""
echo "Next step:"
echo "  qsub -J 1-${n} \\"
echo "       -v \"SAMPLE_LIST=${genlist},FQ_DIR=${fastq_dir},OUT_DIR=/path/to/bam\" \\"
echo "       03_trim_align_markdup.pbs"
