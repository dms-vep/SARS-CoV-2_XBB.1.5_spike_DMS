#!/bin/bash
#
#SBATCH -c 16

snakemake -j 16 --rerun-incomplete --use-conda -s dms-vep-pipeline-3/Snakefile
