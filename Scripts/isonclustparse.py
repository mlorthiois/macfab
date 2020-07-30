import glob
import sys

## steps :
#1 parse tsv
#2 output fastq

tsv=open(snakemake.input[0]+"/final_cluster_origins.tsv")
fastq=open(snakemake.output[0]+"converted.fastq")

for line in tsv.readlines():
    split=line.split("\t")
    fastq.write("@cluster"+split[0]+"-"+split[1]+"\n")
    fastq.write(split[2]+"\n")
    fastq.write("+"+"\n")
    fastq.write(split[3]+"\n")
