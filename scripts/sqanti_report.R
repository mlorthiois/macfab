if (!requireNamespace("ggplot2", quietly = TRUE))
	install.packages("ggplot2", repos="cran.univ-lyon1.fr")
library("ggplot2")
if (!requireNamespace("reshape2", quietly = TRUE))
    install.packages("reshape2")
require("reshape2")

summary <- read.csv(snakemake@input[["summary"]], sep="\t")

genes <- summary[, c('sample', 'novel_genes', 'annotated_genes')]
genes_m <- melt(genes)

trans_colnames <- c(colnames(summary)[1], colnames(summary)[4:12])
trans <- summary[, trans_colnames]
trans_m <- melt(trans)

junc_colnames <- c(colnames(summary)[1], colnames(summary)[13:16])
junc <- summary[, junc_colnames]
junc_m <- melt(junc)

pdf(file=snakemake@output[["report"]], width = 6.5, height = 6.5)
par(mar = c(2, 2, 2, 2))

ggplot(genes_m, aes(sample, value, fill = variable)) + 
    geom_bar(stat = "identity") + 
    labs(fill = "Gene classification", x = "Software and Ref Annotation", y = "Number of Genes") +
    scale_fill_discrete(labels = c("Novel Genes", "Annotated Genes")) +
    theme(legend.position="bottom",
        legend.justification = "left",
        plot.title = element_text(hjust = 0.5, face="bold"), 
        plot.margin = margin(1, 1, 0.5, 0.5, "cm"),
        legend.text = element_text(size=8),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
    ggtitle("Distribution of Novel vs Annotated Genes") +
    guides(fill = guide_legend(title.position = "top"))


ggplot(trans_m, aes(sample, value, fill = variable)) + 
    geom_bar(stat = "identity") +
    labs(fill = "Structural classification", x = "Software and Ref Annotation", y = "Number of Isoforms") +
    scale_fill_discrete(labels = c("FSM", "NNC", "Intergenic", "NIC", "Antisense", 
                                   "Genic", "Fusion", "ISM", "Genic Intron")) +
    theme(legend.position="bottom",
        legend.justification = "left",
        plot.title = element_text(hjust = 0.5, face="bold"), 
        plot.margin = margin(1, 1, 0.5, 0.5, "cm"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        legend.text = element_text(size=8),
        legend.text.align = 0) +
    ggtitle("Distribution of structural isoform classification") +
    guides(fill = guide_legend(title.position = "top"))

ggplot(junc_m, aes(sample, value, fill = variable)) + 
    geom_bar(position="dodge2", stat='identity') + 
    labs(fill = "Splice junctions", x = "Software and Ref Annotation", y = "% of junctions") +
    scale_fill_discrete(labels = c("Known canonical", "Novel canonical", "Known non-canonical", "Novel non-canonical")) +
    theme(legend.position="bottom",
        legend.justification = "left",
        plot.title = element_text(hjust = 0.5, face="bold"), 
        plot.margin = margin(1, 1, 0.5, 0.5, "cm"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        legend.text = element_text(size=8),
        legend.text.align = 0) +
    ggtitle("Distribution of splice junctions") +
    guides(fill = guide_legend(title.position = "top"))

dev.off()