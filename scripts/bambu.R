library("bambu")
library("BSgenome")

gtf_file <- snakemake@input[["gtf"]]
test_bam <- snakemake@input[["bam"]]
fa_file <- snakemake@input[["fa"]]

bambuAnnotations <- prepareAnnotations(gtf_file)

se <- bambu(reads = test_bam, annotations = bambuAnnotations, genome = fa_file)
writeBambuOutput(se, path = "./results/bambu/")
file.rename("./results/bambu/extended_annotations.gtf", snakemake@output[["bambu_gtf"]])
