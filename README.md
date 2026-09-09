# Sunflower Variant Calling Pipeline

> **Running SLURM or want automatic dependency tracking and failure recovery?**
> See the [`snakemake` branch](../../tree/snakemake) for a Snakemake version of this pipeline that runs on both PBS/Torque and SLURM out of the box.

A PBS/Torque job scheduler pipeline for population-scale variant calling in sunflower (*Helianthus annuus*). It produces imputed, population-stratified biallelic SNP sets ready for downstream analysis. The pipeline can start from NCBI SRA accessions **or** from paired-end FASTQs that already exist on the HPC.

## Pipeline Overview

Steps 03–06 are parameterized PBS jobs; all inputs are passed via `qsub -v` so no script needs to be edited between runs. Two entry points feed into step 03:

**Path A — starting from NCBI SRA accessions:**
```
00_prep_manifest        — query NCBI → MANIFEST.txt, GENLIST.txt, SRR_LIST.txt
01_download_srr         — prefetch SRA files from NCBI (one job per SRR)
02_fastq_convert_concat — convert SRA → FASTQ and concatenate runs per genotype (one job per genotype)
  ↓
```
**Path B — starting from existing paired-end FASTQs:**
```
00_prep_local_fastq     — create symlinks + GENLIST.txt from a user-supplied FASTQ manifest (shell script, no job submission needed)
  ↓
```
**Shared steps (both paths):**
```
03_trim_align_markdup   — fastp trim → bwa-mem2 align → samtools markdup → BAM (one job per sample)
04_haplotype_caller     — GATK HaplotypeCaller in GVCF mode (one job per sample, 17 chroms in parallel)
05_genotype_gvcfs       — GATK GenomicsDBImport + GenotypeGVCFs (one job per chromosome)
06_imputation           — bcftools concat → BEAGLE 5.5 imputation → biallelic SNP filter (serial)
```

## Requirements

**Cluster:** PBS/Torque job scheduler with `$SCRATCH` available per job.

**Modules** (loaded automatically within each script):
- `sratoolkit/3.2.1`
- `samtools/1.20`
- `gatk` (v4.6+)
- `bcftools/1.20`
- `parallel`

**Additional software** (expected at `/mmfs1/projects/brent.hulke/software/`):
- `fastp`, `bwa-mem2`, `mosdepth`, `pigz`
- `edirect` (for `00_prep_manifest.pbs` NCBI queries)
- BEAGLE 5.5 jar (path overridable via `BEAGLE_JAR=`)

**Reference genome:** Ha412HOv2 with chloroplast and mitochondrial sequences (HA412HOv2_w_CPMT).

## Porting to a New System

All site-specific paths are consolidated in **`site_config.sh`**. Edit the three values there before running the pipeline on a new cluster — nothing else needs to change:

```bash
# site_config.sh
GENOME_DIR=/path/to/HA412HOv2_w_CPMT    # reference genome directory
SOFTWARE_DIR=/path/to/software           # fastp, bwa-mem2/, mosdepth, pigz
BEAGLE_JAR=/path/to/beagle.jar          # BEAGLE 5.5 jar
```

All three can also be overridden per-submission without editing the file:
```bash
qsub -v "GENOME_DIR=/alt/ref,..." 04_haplotype_caller.pbs
```

**PBS scheduler settings** (`-q queue` and `-W group_list=`) are set in each script's `#PBS` header and can similarly be overridden at submission time:
```bash
qsub -q myqueue -W group_list=mygroup 03_trim_align_markdup.pbs
```

## Input Files

| File | Format | Used by |
|------|--------|---------|
| `SRR_LIST.txt` | one SRR/ERR accession per line | `01_download_srr.pbs` |
| `MANIFEST.txt` | tab-delimited `SAMPLE_NAME<TAB>SRR_ID`, one SRR per row | `02_fastq_convert_concat.pbs` |
| `GENLIST.txt` | unique sample names, one per line | steps 02–05 |
| `population_lists/*.txt` | one sample name per line, one file per population | `06_imputation.pbs` |

`MANIFEST.txt` and `GENLIST.txt` are generated automatically by `00_prep_manifest.pbs`, or can be supplied manually. See `MANIFEST.txt` and `GENLIST.txt` in this repository for examples.

## Quick Start

```bash
# 1. Generate MANIFEST and GENLIST from an SRR list (skip if you already have them)
qsub -v "INPUT=/path/to/srr_list.txt,INPUT_TYPE=srr,OUT_DIR=/path/to/lists" \
     00_prep_manifest.pbs

# 2. Download SRA files (set -J to number of SRRs)
N=$(wc -l < /path/to/lists/SRR_LIST.txt)
qsub -J 1-$N \
     -v "SRR_LIST=/path/to/lists/SRR_LIST.txt,SRA_DIR=/path/to/sra" \
     01_download_srr.pbs

# 3. Convert and concatenate FASTQs (set -J to number of genotypes)
N=$(wc -l < /path/to/lists/GENLIST.txt)
qsub -J 1-$N \
     -v "MANIFEST=/path/to/lists/MANIFEST.txt,GENLIST=/path/to/lists/GENLIST.txt,\
         SRA_DIR=/path/to/sra,FASTQ_DIR=/path/to/fastq,\
         DRIVER=$(pwd)/02_fastq_convert_concat_driver.sh" \
     02_fastq_convert_concat.pbs

# 4. Trim, align, mark duplicates
qsub -J 1-$N \
     -v "SAMPLE_LIST=/path/to/lists/GENLIST.txt,FQ_DIR=/path/to/fastq,OUT_DIR=/path/to/bam" \
     03_trim_align_markdup.pbs

# 5. Call variants (GVCF mode)
qsub -J 1-$N \
     -v "SAMPLE_LIST=/path/to/lists/GENLIST.txt,BAM_DIR=/path/to/bam,OUT_DIR=/path/to/gvcf" \
     04_haplotype_caller.pbs

# 6. Joint genotyping (always 17 jobs — one per chromosome)
qsub -J 1-17 \
     -v "SAMPLE_LIST=/path/to/lists/GENLIST.txt,GVCF_DIR=/path/to/gvcf,\
         COHORT=MyProject,OUT_DIR=/path/to/vcf" \
     05_genotype_gvcfs.pbs

# 7. Imputation by population
qsub -v "VCF_DIR=/path/to/vcf,COHORT=MyProject,\
         POP_LIST_DIR=/path/to/population_lists,OUT_DIR=/path/to/imputed" \
     06_imputation.pbs
```

To test a single sample before a full run, override the array range: `qsub -J 1-1 -v ...`

## Quick Start — local FASTQs (Path B)

If paired-end FASTQs already exist on the HPC, skip steps 00–02. Instead:

**1. Create a manifest** (tab-delimited, one sample per line):
```
# SAMPLE_NAME	R1_PATH	R2_PATH
HA89	/path/to/HA89_R1.fastq.gz	/path/to/HA89_R2.fastq.gz
HA412	/path/to/HA412_L001_R1_001.fastq.gz	/path/to/HA412_L001_R2_001.fastq.gz
```
Lines beginning with `#` and blank lines are ignored. Absolute paths are recommended. Source files may follow any naming convention; the script normalises them.

**2. Run the prep script:**
```bash
bash 00_prep_local_fastq.sh /path/to/manifest.tsv /path/to/fastq_dir /path/to/list_dir
```

**3. Continue with steps 03–06 as normal** (the script prints the exact `qsub` command):
```bash
N=$(wc -l < /path/to/list_dir/GENLIST.txt)
qsub -J 1-$N \
     -v "SAMPLE_LIST=/path/to/list_dir/GENLIST.txt,FQ_DIR=/path/to/fastq_dir,OUT_DIR=/path/to/bam" \
     03_trim_align_markdup.pbs
```

## Script Details

### 00_prep_manifest.pbs
Queries NCBI for sample metadata and builds the manifest files needed by downstream scripts. Accepts either a flat list of SRR accessions or a BioProject/study accession. SRRs are queried in batches of 100 to stay within NCBI rate limits.

**Args:** `INPUT`, `INPUT_TYPE` (`srr` or `bioproject`), `OUT_DIR`  
**Optional:** `EMAIL` (NCBI queries; default: `joseph.barham@ndsu.edu`)  
**Outputs:** `MANIFEST.txt`, `GENLIST.txt`, `SRR_LIST.txt`

---

### 00_prep_local_fastq.sh *(alternative to steps 00–02)*
Creates symlinks in `FASTQ_DIR` and writes `GENLIST.txt` from a user-supplied manifest of existing paired-end FASTQs. No job submission is needed — this is a lightweight shell script that runs on the login node. Source files may follow any naming convention; the symlinks are normalised to `{sample}_R1.{ext}` / `{sample}_R2.{ext}` so that `03_trim_align_markdup.pbs` picks them up without modification.

**Args (positional):** `MANIFEST`, `FASTQ_DIR`, `LIST_DIR`  
**Manifest format:** tab-delimited `SAMPLE_NAME<TAB>R1_PATH<TAB>R2_PATH`; `#` comment lines and blank lines are skipped  
**Outputs:** symlinks in `FASTQ_DIR`, `{LIST_DIR}/GENLIST.txt`  
**Note:** Source files must be gzip-compressed (`.fastq.gz` or `.fq.gz`); the script warns if uncompressed files are detected.

---

### 01_download_srr.pbs
Downloads SRA files using `prefetch`. Skips accessions already present in `SRA_DIR`, so the job array can be safely resubmitted after partial failures.

**Args:** `SRR_LIST`, `SRA_DIR`

---

### 02_fastq_convert_concat.pbs + 02_fastq_convert_concat_driver.sh
**Array indexed by genotype** — each job processes all SRRs for one sample, then concatenates them. This is the most storage-sensitive step; several design choices keep peak disk usage low:

- `fasterq-dump` writes uncompressed temp files to `$SCRATCH` (not project quota)
- Each SRR is compressed with `pigz` immediately after dumping; uncompressed files are deleted before the next SRR begins
- Multiple `.fastq.gz` files are concatenated using `cat` (gzip stream concatenation — valid per RFC 1952, compatible with all downstream tools)
- Per-SRR `.fastq.gz` files are deleted once the concatenated output passes size and `gzip -t` validation
- Jobs skip genotypes whose output already exists (safe to resubmit)
- Storage is checked before each `fasterq-dump`; the job exits with a clear message if space is insufficient

**Chunking for storage limits:** if you are near your quota, process genotypes in batches by passing a subset `GENLIST` and adjusting `-J` accordingly. Jobs that find finished output skip automatically, so chunks can overlap.

**Args:** `MANIFEST`, `GENLIST`, `SRA_DIR`, `FASTQ_DIR`, `DRIVER`  
**Optional:** `DELETE_SRA` (default `false`), `MIN_FREE_GB` (default `100`), `INFLATION_FACTOR` (default `8`)  
**Outputs:** `{SAMPLE}_{SRR1}[_{SRR2}...]_1.fastq.gz` and `_2.fastq.gz` per genotype

---

### 03_trim_align_markdup.pbs
Runs fastp → bwa-mem2 → samtools collate/fixmate/sort/markdup in a single pipe to avoid writing intermediate BAMs. Sort temp files go to `$SCRATCH`.

**Args:** `SAMPLE_LIST`, `FQ_DIR`, `OUT_DIR`  
**Outputs per sample:** `{sample}.bam`, `{sample}.bam.bai`, QC dirs `fastp_out/`, `markdup_out/`, `mosdepth_out/`, `flagstat_out/`

---

### 04_haplotype_caller.pbs
Runs GATK HaplotypeCaller in GVCF mode. All 17 chromosomes are processed in parallel within each job using GNU parallel. Repeat regions are excluded via a BED file. Heterozygosity priors are tuned for sunflower (0.01 ± 0.1; Todesco et al. 2020).

**Args:** `SAMPLE_LIST`, `BAM_DIR`, `OUT_DIR`  
**Outputs:** `{OUT_DIR}/{chrom}/{sample}.{chrom}.g.vcf.gz`

---

### 05_genotype_gvcfs.pbs
Runs GenomicsDBImport followed by GenotypeGVCFs, one chromosome per array job. Accepts a single-column `GENLIST` as the sample list and builds the two-column GATK sample map internally. The GenomicsDB workspace is written to `$SCRATCH`.

**Args:** `SAMPLE_LIST`, `GVCF_DIR`, `COHORT`, `OUT_DIR`  
**Outputs:** `{OUT_DIR}/{COHORT}.{chrom}.vcf.gz`

---

### 06_imputation.pbs
Concatenates all 17 chromosome VCFs, subsets by population, runs BEAGLE 5.5 phasing and imputation, then filters to biallelic SNPs. All intermediate files are deleted on completion.

**Args:** `VCF_DIR`, `COHORT`, `POP_LIST_DIR`, `OUT_DIR`  
**Optional:** `BEAGLE_JAR`, `N_THREADS` (default `16`), `MEM_GB` (default `120`)  
**Outputs:** `{COHORT}.{POP_NAME}.imputed.biallelic_snps.vcf.gz` (+ `.tbi`) per population

## License
This work is released into the public domain under the CC0 1.0 Universal (CC0 1.0) Public Domain Dedication. To the extent possible under law, all copyright and related or neighboring rights to this work have been waived. This work is published from the United States. For more information, see <https://creativecommons.org/publicdomain/zero/1.0/>.
