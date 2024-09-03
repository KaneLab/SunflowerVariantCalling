# GATK4.6.0_generateGVCFs_HA412HOv2cpmt.pbs

## Description
This PBS script implements a pipeline for calling variants on multiple sunflower samples using GATK 4.6.0 using the HA412HOv2_w_cpmt reference genome. It performs automated read quality control, alignment, and variant calling, utilizing parallelization for efficiency. The pipeline starts with fastq files and produces gVCF files for each sample.

## Author
Brian Smart and Kyle Keepers

## Date
September 3, 2024

## Prerequisites
- PBS job scheduling system
- GATK 4.6.0
- fastp
- bwa-mem2
- samtools 1.20
- bcftools 1.20
- GNU parallel

## Installation
1. Ensure all required software is installed and accessible in your environment.
2. Clone this repository or download the `generate412GVCFs.pbs` script.

## Usage
Submit the job using qsub with the following command:

```bash
qsub -v FASTQ_DIR="/path/to/fastqs",OUTPUT_DIR="/path/to/output",PL="ILLUMINA" GATK4.6.0_generateGVCFs_HA412HOv2cpmt.pbs
```

### Parameters:
- `FASTQ_DIR`: Path to the directory containing input fastq files
- `OUTPUT_DIR`: Path to the directory where output files will be saved
- `PL`: Sequencing platform (e.g., ILLUMINA)

### Input Requirements:
- Fastq files must be uniquely named and end with `1.f*q.gz` and `2.f*q.gz` for paired-end reads.
- All samples to be processed must be in a single folder specified by `FASTQ_DIR`.

## Pipeline Steps
1. **Quality Control and Alignment:**
   - Performs read quality control using fastp
   - Aligns reads to the reference genome using bwa-mem2
   - Sorts, marks duplicates, and indexes the resulting BAM files

2. **Variant Calling:**
   - Uses GATK HaplotypeCaller to call variants
   - Processes the genome in parallel using regions defined in a BED file
   - Concatenates results into a single gVCF file per sample

## Output

- Sorted and indexed BAM files for each sample
- gVCF files for each sample
- Index files for the gVCFs

## Resource Requirements

- 53 CPUs
- 477GB of memory
- Maximum walltime: 168 hours

## Customization

You may need to modify the following paths in the script:

- `pathToFastp`
- `pathTobwamem2`
- `pathToBED`
- `pathToReference`
- `FileListCopy`

Ensure these paths are correct for your system before running the pipeline.

## Notes

- The script uses a BED file to split the genome into 53 regions for parallel processing.
- It's designed to work with the HA412HOv2_w_CPMT reference genome.

## Support

For issues or questions, please open an issue in this GitHub repository.
