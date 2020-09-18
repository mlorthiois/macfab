library("bambu")
library("BSgenome")

gtf_file <- snakemake@input[["gtf"]]
test_bam <- snakemake@input[["bam"]]
fa_file <- snakemake@input[["fa"]]

bambuAnnotations <- prepareAnnotationsFromGTF(gtf_file)

se <- bambu(reads = test_bam, annotations = bambuAnnotations, genomeSequence = fa_file)
writeBambuOutput(se, path = "./results/bambu/")
file.rename("./results/bambu/transcript_exon.gtf", snakemake@output[["bambu_gtf"]])
