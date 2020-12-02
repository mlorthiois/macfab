#!/bin/bash

#SBATCH --job-name=macfab
#SBATCH --output=slurm_macfab.out

. /local/env/envsnakemake-5.20.1.sh

snakemake --use-conda \
    --cores 30 \
    --rerun-incomplete \
    --conda-frontend mamba \
    --cluster "sbatch --cpus-per-task={threads} --mem={resources.ram} --output=slurm_{rule}-%j.out --job-name={rule}" \
    -j 15 \
    --latency-wait 60 \
    --keep-going \
    --restart-times 3
