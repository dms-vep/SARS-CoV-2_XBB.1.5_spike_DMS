#!/bin/bash
#
#SBATCH -c 32

snakemake -j 32 --use-conda  --rerun-incomplete
