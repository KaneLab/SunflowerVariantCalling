#!/usr/bin/env bash
# site_config.sh — sourced by PBS scripts at runtime.
#
# To port this pipeline to a new system, edit the three values below.
# Everything else in the pipeline derives from these roots.
#
# All three variables support qsub-time override:
#   qsub -v "GENOME_DIR=/my/ref,..." script.pbs
# The values here are used only when the variable is not already set.
#
# GENOME_DIR   — directory containing the Ha412HOv2_w_CPMT reference genome,
#                chromosome list, and filtered-repeat BED file
# SOFTWARE_DIR — directory containing fastp, bwa-mem2/, mosdepth, pigz
# BEAGLE_JAR   — path to the BEAGLE 5.5 jar file

GENOME_DIR="${GENOME_DIR:-/mmfs1/projects/brent.hulke/sunflower_reference_genomes/HA412HOv2_w_CPMT}"
SOFTWARE_DIR="${SOFTWARE_DIR:-/mmfs1/projects/brent.hulke/software}"
BEAGLE_JAR="${BEAGLE_JAR:-/mmfs1/home/brian.smart/projects/software/beagle.27Feb25.75f.jar}"
