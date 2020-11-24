library(ggplot2)

data <- read.table("results_bench.csv", header=TRUE, sep=",")

ggplot(data, aes(reads, isoforms, colour = tool)) + 
  facet_wrap(~sample)+
  geom_point() +
  geom_line() +
  # geom_smooth(method="lm",formula=y~log(x), se=FALSE, fullrange=TRUE) +
  scale_x_continuous(name="Number of reads (millions)") +
  scale_y_continuous(name="Number of isoforms", labels=function(x) format(x, big.mark = ",", scientific = FALSE)) +
  theme_bw()
ggsave(file="bench_tools.pdf", width=8, height=6, dpi=300)

ggplot(data, aes(reads, isoforms, colour = sample)) + 
  facet_wrap(~tool)+
  geom_point() +
  geom_line() +
  # geom_smooth(method="lm",formula=y~log(x), se=FALSE, fullrange=TRUE) +
  scale_x_continuous(name="Number of reads (millions)") +
  scale_y_continuous(name="Number of isoforms", labels=function(x) format(x, big.mark = ",", scientific = FALSE)) +
  theme_bw()
ggsave(file="bench_tools_bytools.pdf", width=8, height=6, dpi=300)