library("bambu")
library("BSgenome")

gtf.file <- snakemake@input[["gtf"]]
test.bam <- snakemake@input[["bam"]]
fa.file <- snakemake@input[["fa"]]

bambuAnnotations <- prepareAnnotationsFromGTF(gtf.file)

se <- bambu(reads = test.bam, annotations = bambuAnnotations, genomeSequence = fa.file)
writeBambuOutput(se, path = "./")
file.rename("transcript_exon.gtf", snakemake@output[["o_name"]])
