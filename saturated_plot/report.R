library(ggplot2)

data <- read.table("results_bench.csv", header=TRUE, sep=",")

isoforms <- ggplot(data, aes(reads, isoforms, colour = tool)) + 
  facet_wrap(~sample)+
  geom_point() +
  geom_line() +
  # geom_smooth(method="lm",formula=y~log(x), se=FALSE, fullrange=TRUE) +
  scale_x_continuous(name="Number of reads (millions)") +
  scale_y_continuous(name="Number of isoforms", labels=function(x) format(x, big.mark = ",", scientific = FALSE)) +
  theme_bw()

genes <- ggplot(data, aes(reads, genes, colour = tool)) + 
  facet_wrap(~sample)+
  geom_point() +
  geom_line() +
  # geom_smooth(method="lm",formula=y~log(x), se=FALSE, fullrange=TRUE) +
  scale_x_continuous(name="Number of reads (millions)") +
  scale_y_continuous(name="Number of genes", labels=function(x) format(x, big.mark = ",", scientific = FALSE)) +
  theme_bw()

pdf("results.pdf", width=8, height=6)
print(list(isoforms, genes))
dev.off()

print("PDF created")