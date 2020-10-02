if (!requireNamespace("tidyverse", quietly = TRUE))
	install.packages("tidyverse", repos="cran.univ-lyon1.fr")
library("tidyverse")
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

pdf(file=snakemake@output[["report"]], width = 8, height = 7)

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

junc_m %>% group_by(sample) %>% 
    summarise(total = sum(value)) %>% 
    left_join(junc_m) %>% 
    mutate(percent = round(100*value/total, 1)) %>% 
    ggplot() + 
        aes(x=sample, y = value, fill = variable, label = percent) + 
        geom_bar(stat = "identity", position = "dodge", colour="black") +
        geom_text(aes(label=paste(percent, '%', sep = ""),
                    vjust=-0.5),
                    position = position_dodge(0.9),
                    size=2.8) +
        labs(fill = "Junction classification", x = "Software and Ref Annotation", y = "Number of Junctions") +
        scale_fill_discrete(labels = c("Known Canonical", "Novel Canonical",
                                        "Known non-canonical", "Novel non-canonical")) +
        theme(legend.position="bottom",
            legend.justification = "left",
            plot.title = element_text(hjust = 0.5, face="bold"), 
            plot.margin = margin(1, 1, 0.5, 0.5, "cm"),
            axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
            legend.text = element_text(size=8),
            legend.text.align = 0) +
        ggtitle("Distribution of junction classification") +
        guides(fill = guide_legend(title.position = "top"))

dev.off()