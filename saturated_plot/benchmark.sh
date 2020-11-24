#!/bin/bash
    	
#SBATCH --job-name=bench_macfab
#SBATCH --cpus-per-task=20
#SBATCH --mem=60G
#SBATCH --output=slurm_benchmacfab.out

echo "Launch benchmark"
numbers="500000"
samples="/groups/dog/mlorthiois/manopore/guppy_4.0.14/Bear/results/guppy/fastq/Bear.fastq.gz
        /home/genouest/cnrs_umr6290/mlorthiois/workspace/manopore/guppy_4.0.14/CML10/results/guppy/fastq/CML10.fastq.gz
        /home/genouest/cnrs_umr6290/mlorthiois/workspace/manopore/guppy_4.0.14/popsi/results/guppy/fastq/popsi.fastq.gz
        /home/genouest/cnrs_umr6290/mlorthiois/workspace/manopore/guppy_4.0.14/twiny/results/guppy/fastq/twiny.fastq.gz"


. /local/env/envsnakemake-5.20.1.sh

for sample in $samples; do
    echo "Run analysis on $sample"
    sample_id=$(basename $sample | cut -f1 -d '.')
    if ! [ -d $sample_id ]; then
        mkdir -p $sample_id
    fi
    cd $sample_id
    
    if [ ! -f "./$sample_id.fastq" ]; then
        echo "Unzip $sample"
        gunzip -c $sample > $sample_id.fastq
    fi
    
    count=$(\cat $sample_id.fastq | wc -l)
    
    for number in $numbers; do
        if (( $number*4 < count)); then
            if ! [ -d $number ]; then
                echo "Download macfab"
                git clone -q git@gitlab.com:bioinfog/macfab.git $number
            fi
            
            cd $number
            
            if [ ! -f "./$number.fastq" ]; then
                echo "Create subsample $number for $sample_id ..."
                head -n $((number*4)) ../$sample_id.fastq > $sample_id_$number.fastq
            fi

            if ! grep -q "$number.fastq" config.yaml; then
                grep -v "fastq_path" config.yaml > temp.yaml
                echo "fastq_path : '$number.fastq'" >> temp.yaml
                mv temp.yaml config.yaml
            fi

            snakemake --use-conda --cores 20 --rerun-incomplete --conda-frontend mamba --conda-prefix /home/genouest/cnrs_umr6290/mlorthiois/workspace/macfab/benchmark_tools/conda

            cd ..
        fi
    done
    cd ..
done

chmod +x parse_isoforms.sh && ./parse_isoforms.sh

conda activate tidyverse

Rscript report.R
