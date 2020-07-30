if (!requireNamespace("devtools", quietly = TRUE))
	install.packages("devtools")
devtools::install_github("GoekeLab/bambu")

if (!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")
BiocManager::install("BSgenome")

library(bambu)
library("BSgenome")

gtf.file <- snakemake@input[["gtf"]]
test.bam <- snakemake@input[["bam"]]
fa.file <- snakemake@input[["fa"]]

bambuAnnotations <- prepareAnnotationsFromGTF(gtf.file)

se <- bambu(reads = test.bam, annotations = bambuAnnotations, genomeSequence = fa.file)
writeBambuOutput(se, path = snakemake@output[["o_dir"]])
file.rename(paste(snakemake@input[["0"]], "transcript_exon.gtf", sep=""), paste(snakemake@input[["0"]], snakemake@output[["o_name"]], sep=""))
