#!/bin/bash
    	
#SBATCH --job-name=macfab_analysis
#SBATCH --nodelist=cl1n036
#SBATCH --cpus-per-task=30
#SBATCH --mem=70G
#SBATCH --output=slurm_macfab.out

. /local/env/envsnakemake-5.20.1.sh

snakemake --use-conda --cores 30 --rerun-incomplete
