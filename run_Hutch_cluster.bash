#!/bin/bash
#
#SBATCH -c 16

snakemake -j 16 --use-conda -s dms-vep-pipeline-3/Snakefile
