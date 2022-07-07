### Title: Figuring out why abundant windows are all left shifted ### 
### Date: 12/20/21
library(gtools)

#Create dataframe with all 96 props for each window 
i <- unique(allele_counts$MST)[1]
mst_data <- read.table(paste0("/Users/kevinliao/Desktop/Michigan Research/Selection/tajD_MST_corrected/tajD_", i, " (2020_02_29 01_13_00 UTC).txt"), header=TRUE)
mst_data$tajD_sig <- ifelse(mst_data$taj_D < quantile(mst_data$taj_D, probs = 0.05), 1, 0)
top_windows_df <- mst_data[,c(1,2,9)]

for(i in unique(allele_counts$MST)){
  print(i)
  if(i == '.'){next}
  
  mst_data <- read.table(paste0("/Users/kevinliao/Desktop/Michigan Research/Selection/tajD_MST_corrected/tajD_", i, " (2020_02_29 01_13_00 UTC).txt"), header=TRUE)
  mst_data$tajD_sig <- ifelse(mst_data$taj_D < quantile(mst_data$taj_D, probs = 0.05), 1, 0)
  mst_data$prop_weird_sig <- ifelse(mst_data$prop_weird > quantile(mst_data$prop_weird, probs = 0.90), 1, 0)
  keep <- mst_data[, c(1,2,10)]
  colnames(keep) <- c('chr', 'window', paste0("abundant_", i))
  
  top_windows_df <- merge(top_windows_df, keep, by=c('chr','window'), all.x=TRUE)
}

top_windows_df$num_abundant <- rowSums(top_windows_df[,4:99], na.rm=TRUE)

t.test(subset(top_windows_df, top_windows_df$tajD_sig == 1)$num_abundant, subset(top_windows_df, top_windows_df$tajD_sig == 0)$num_abundant)
tapply(top_windows_df$num_abundant, top_windows_df$tajD_sig, mean)

sig_tajD <- subset(top_windows_df, top_windows_df$tajD_sig == 1)
nonsig_tajD <- subset(top_windows_df, top_windows_df$tajD_sig == 0)

plot1 <- ggplot(top_windows_df, aes(x=num_abundant, fill=as.factor(tajD_sig))) + geom_density(alpha = 0.5) +
  xlab("# Abundant Subtypes in 100Kb Windows") + ylab("Density") + 
  ggtitle("Distribution of # of Abundant Subtypes\nBy Significant Tajima's D Status") + 
  theme(plot.title = element_text(hjust = 0.5)) + scale_fill_discrete(label=c('No', "Yes"), name=("D Significant"))
plot1

#Look at which ones are in
num_abundant <- as.data.frame(colSums(sig_tajD[,4:99], na.rm=TRUE))
num_abundant$MST <- substr(rownames(num_abundant), 10, 16)
num_abundant$MT <- substr(num_abundant$MST, 1,3)
tapply(num_abundant$`colSums(sig_tajD[, 4:99], na.rm = TRUE)`, num_abundant$MT, mean)

#Look at correlation between genomewide D and num abundat windows 
merge_abundant <- merge(merged_data, num_abundant, by='MST')
cor.test(merge_abundant$D, merge_abundant$`colSums(sig_tajD[, 4:99], na.rm = TRUE)`)

plot2 <- ggplot(merge_abundant, aes(x=D, y=`colSums(sig_tajD[, 4:99], na.rm = TRUE)`)) + geom_point() +
  xlab("Genome wide Tajima's D ") + ylab("# Windows Abundant in Subtype") + 
  ggtitle("Number of 100Kb Windows Abundant in Subtype\nby Subtype's Genome Wide Tajima's D\n*Significant D Windows Only") + 
  theme(plot.title = element_text(hjust = 0.5)) + geom_smooth(method='lm', se=FALSE)

ggarrange(plot1, plot2, ncol=1,labels=c("a)", "b)"))
