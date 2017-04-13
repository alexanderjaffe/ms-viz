setwd("~/Desktop")

library(ggplot2)

samples1 = read.table("sample_count_results.txt", sep="\t", header=T, check.names = F)
samples2 = read.table("sample_count_results_wild.txt", sep="\t", header=T, check.names = F)

samples = rbind(samples1, samples2)

ggplot(data=samples, aes(x=sample_count, y=unique_cmpds, colour=group)) + geom_line(size=2) + geom_point(size=4) + theme_bw() + theme(axis.text=element_text(size=24), axis.title=element_text(size=24), legend.text = element_text(size=18), legend.title = element_text(size=18))
filter1 = read.table("filter_threshold_results_50.txt", sep="\t", header=T, check.names = F)

ggplot(data=filter1, aes(x=threshold, y=unique_cmpds)) + geom_line(size=2) + geom_point(size=4) + theme_bw() + theme(axis.text=element_text(size=24), axis.title=element_text(size=24), legend.text = element_text(size=18), legend.title = element_text(size=18))
