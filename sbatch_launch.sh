#!/bin/bash

source /local/env/envsnakemake-5.20.1.sh

#time snakemake --use-conda -s Snakefile.py --cores 4 --cluster "sbatch --cpus-per-task={threads} --mem={resources.ram}" -j 5
time snakemake --use-conda -s Snakefile.py --cores 4 --cluster "sbatch --cpus-per-task={threads} --mem={resources.ram}" -j 5

# sbatch -J Snakefile_slurm -o Snakefile_slurm."%j".out -e Snakefile_slurm."%j".err --mail-user tderrien@univ-rennes1.fr --mail-type=ALL sbatch_launch.sh
