if (!requireNamespace("ggplot2", quietly = TRUE))
	install.packages("ggplot2")

library(ggplot2)

Sensitivity.gffparse <- read.delim(snakemake@input[["Sensitivity"]])
Values.gffparse <- read.delim(snakemake@input[["Values"]])

#ggplot(data=gffdata, aes(x=Type, y=Sensitivity, fill=Source)) + geom_bar(stat="identity", position=position_dodge()) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pdf("Graph.recap.pdf")
ggplot(data=Sensitivity.gffparse, aes(x=Feature, y=Value)) + geom_point(aes(colour=Source), size=2) + facet_grid(Metric ~ Annotation) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + coord_flip() + scale_color_brewer(palette = "Dark2")
ggplot(data=Values.gffparse, aes(x=Feature, y=Value)) + geom_point(aes(colour=Source), size=2) + facet_grid(Annotation ~ Type) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + coord_flip() + scale_color_brewer(palette = "Dark2")
dev.off()
