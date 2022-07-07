### Make plot for combining mut rate/gBGC with distribution of variants. Also permuted data too ### 
### Date: 1/13/22
library(gtools)

i <- unique(allele_counts$MST)[1]
mst_data <- read.table(paste0("/Users/kevinliao/Desktop/Michigan Research/Selection/tajD_MST_corrected/tajD_", i, " (2020_02_29 01_13_00 UTC).txt"), header=TRUE)
mst_data$tajD_group <- quantcut(mst_data$taj_D, q=20)
top_windows_df <- mst_data[,c(1,2,9)]

for(i in unique(allele_counts$MST)){
  print(i)
  if(i == '.'){next}
  
  mst_data <- read.table(paste0("/Users/kevinliao/Desktop/Michigan Research/Selection/tajD_MST_corrected/tajD_", i, " (2020_02_29 01_13_00 UTC).txt"), header=TRUE)
  mst_data$tajD_sig <- ifelse(mst_data$taj_D < quantile(mst_data$taj_D, probs = 0.05), 1, 0)
  mst_data$prop_weird_sig <- ifelse(mst_data$prop_weird > quantile(mst_data$prop_weird, probs = 0.90), 1, 0)
  keep <- mst_data[, c(1,2,10)]
  colnames(keep) <- c('chr', 'window', i)
  
  top_windows_df <- merge(top_windows_df, keep, by=c('chr','window'), all.x=TRUE)
}

top_windows_df[is.na(top_windows_df)] <- 0

#Figure out which subtypes fall into following categories: CpG -> TpG, SW low mut, SW high mut, etc
median_mut_rate <- median(merged_data$ERV_rel_rate)
SW_types <- c("C_T", "C_A"); WS_types <- c("A_C", "A_G")
neutral_types <- c("C_G", "A_T")
CpG_TpG <- c("C_T.ACG", "C_T.CCG", "C_T.GCG", "C_T.TCG")

SW_low_mut <- subset(no_CpG, substr(no_CpG$MST, 1, 3) %in% SW_types & no_CpG$ERV_rel_rate <= median_mut_rate)$MST
SW_high_mut <- subset(no_CpG, substr(no_CpG$MST, 1, 3) %in% SW_types & no_CpG$ERV_rel_rate > median_mut_rate)$MST

WS_low_mut <- subset(no_CpG, substr(no_CpG$MST, 1, 3) %in% WS_types & no_CpG$ERV_rel_rate <= median_mut_rate)$MST
WS_high_mut <- subset(no_CpG, substr(no_CpG$MST, 1, 3) %in% WS_types & no_CpG$ERV_rel_rate > median_mut_rate)$MST

neutral_low_mut <- subset(no_CpG, substr(no_CpG$MST, 1, 3) %in% neutral_types & no_CpG$ERV_rel_rate <= median_mut_rate)$MST
neutral_high_mut <- subset(no_CpG, substr(no_CpG$MST, 1, 3) %in% neutral_types & no_CpG$ERV_rel_rate > median_mut_rate)$MST

top_windows_df$CpG_TpG <- rowSums(top_windows_df[,CpG_TpG])
top_windows_df$SW_low_mut <- rowSums(top_windows_df[,SW_low_mut])
top_windows_df$SW_high_mut <- rowSums(top_windows_df[,SW_high_mut])
top_windows_df$WS_low_mut <- rowSums(top_windows_df[,WS_low_mut])
top_windows_df$WS_high_mut <- rowSums(top_windows_df[,WS_high_mut])
top_windows_df$neutral_low_mut <- rowSums(top_windows_df[,neutral_low_mut])
top_windows_df$neutral_high_mut <- rowSums(top_windows_df[,neutral_high_mut])
head(top_windows_df)

CpG_TpG_quantiles <- as.data.frame(tapply(top_windows_df$CpG_TpG, top_windows_df$tajD_group, mean))
SW_low_mut_quantiles <- as.data.frame(tapply(top_windows_df$SW_low_mut, top_windows_df$tajD_group, mean))
SW_high_mut_quantiles <- as.data.frame(tapply(top_windows_df$SW_high_mut, top_windows_df$tajD_group, mean))
WS_low_mut_quantiles <- as.data.frame(tapply(top_windows_df$WS_low_mut, top_windows_df$tajD_group, mean))
WS_high_mut_quantiles <- as.data.frame(tapply(top_windows_df$WS_high_mut, top_windows_df$tajD_group, mean))
neutral_low_mut_quantiles <- as.data.frame(tapply(top_windows_df$neutral_low_mut, top_windows_df$tajD_group, mean))
neutral_high_mut_quantiles <- as.data.frame(tapply(top_windows_df$neutral_high_mut, top_windows_df$tajD_group, mean))

all_categories <- cbind(CpG_TpG_quantiles, SW_low_mut_quantiles, SW_high_mut_quantiles,
                        WS_low_mut_quantiles, WS_high_mut_quantiles, neutral_low_mut_quantiles, neutral_high_mut_quantiles)
all_categories$D_quantile <- as.factor(rownames(all_categories))
all_categories$D_quantile_rank <- rank(all_categories$D_quantile)
colnames(all_categories) <- c("CpG_TpG", "SW_low", "SW_high", "WS_low", "WS_high", "Neutral_low", "Neutral_high", "D_quantiles", "Quantile_rank")

library(tidyr)
library(dplyr)
categories_long <- gather(all_categories, group, num_abundant, CpG_TpG:Neutral_high, factor_key=TRUE)
categories_long$gBGC <- sub("\\_.*", "", categories_long$group)
categories_long$mut_rate <- sub(".*\\_", "", categories_long$group)
categories_long$gBGC <- ifelse(categories_long$gBGC == 'CpG', 'CpG TpG', categories_long$gBGC)
categories_long$mut_rate <- ifelse(categories_long$mut_rate == 'TpG', 'CpG TpG', categories_long$mut_rate)

ggplot(categories_long, aes(x=D_quantiles, y=num_abundant, fill=group, group=group)) + geom_area(position='stack')

ggplot(categories_long, aes(x=D_quantiles, y=num_abundant, fill=group)) + 
  geom_bar(stat='identity') 

ggplot(data=categories_long, aes(x=D_quantiles, y=num_abundant, fill=gBGC, pattern=mut_rate)) +
  geom_bar_pattern(position='stack', stat='identity',
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025) + 
  scale_pattern_manual("Mutation Rate", values = c(high = "stripe", low = "none", `CpG TpG`='circle'),
                       labels=c("CpG TpG", "High", "Low")) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
        fill = guide_legend(override.aes = list(pattern = "none"))) +
  xlab("D Quantiles") + ylab("Average # Abundant Subtypes") +
  ggtitle("Average # Abundant Subtypes by D Quantiles across Subtype Categories") + 
  scale_x_discrete(limits = rev) + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#Make proportions plot suggested by Jean 
all_categories$total_abun <- rowSums(all_categories[,1:7])
all_categories_prop <- all_categories[,1:7]/rowSums(all_categories[,1:7])
all_categories_prop$tajD_group <- rownames(all_categories)
all_categories_prop$D_quantile_rank <- rank(desc(all_categories_prop$tajD_group))

cat_only <- all_categories_prop[,1:7]
(cat_only[nrow(cat_only),] - cat_only[1,]) / cat_only[1,]

prop_long <- gather(all_categories_prop, cat_group, prop_abundant, CpG_TpG:Neutral_high, factor_key=TRUE)

ggplot(data=prop_long, aes(x=reorder(tajD_group, desc(tajD_group)), y=prop_abundant, color=cat_group)) + 
  geom_line(aes(group=cat_group)) +
  geom_point() +
  xlab("D Quantiles") + ylab("Proportion of Abundant Subtypes") + 
  ggtitle("Proportion of Abundant by D Quantiles across Subtype Categories") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_discrete(name = "Category Group")












###################### REPEAT FOR PERMUTED DATA ################################

i <- unique(allele_counts$MST)[1]
mst_data <- read.table(paste0("/Users/kevinliao/Desktop/Michigan Research/Selection/tajD_permuted/tajD_", i, " (2020_02_29 01_13_00 UTC).txt"), header=TRUE)
mst_data$tajD_group <- quantcut(mst_data$taj_D, q=20)
top_windows_df <- mst_data[,c(1,2,9)]

for(i in unique(allele_counts$MST)){
  print(i)
  if(i == '.'){next}
  
  mst_data <- read.table(paste0("/Users/kevinliao/Desktop/Michigan Research/Selection/tajD_permuted/tajD_", i, " (2020_02_29 01_13_00 UTC).txt"), header=TRUE)
  mst_data$tajD_sig <- ifelse(mst_data$taj_D < quantile(mst_data$taj_D, probs = 0.05, na.rm=TRUE), 1, 0)
  mst_data$prop_weird_sig <- ifelse(mst_data$prop_weird > quantile(mst_data$prop_weird, probs = 0.90), 1, 0)
  keep <- mst_data[, c(1,2,10)]
  colnames(keep) <- c('chr', 'window', i)
  
  top_windows_df <- merge(top_windows_df, keep, by=c('chr','window'), all.x=TRUE)
}

top_windows_df[is.na(top_windows_df)] <- 0

#Figure out which subtypes fall into following categories: CpG -> TpG, SW low mut, SW high mut, etc
median_mut_rate <- median(merged_data$ERV_rel_rate)
SW_types <- c("C_T", "C_A"); WS_types <- c("A_C", "A_G")
neutral_types <- c("C_G", "A_T")
CpG_TpG <- c("C_T.ACG", "C_T.CCG", "C_T.GCG", "C_T.TCG")

SW_low_mut <- subset(no_CpG, substr(no_CpG$MST, 1, 3) %in% SW_types & no_CpG$ERV_rel_rate <= median_mut_rate)$MST
SW_high_mut <- subset(no_CpG, substr(no_CpG$MST, 1, 3) %in% SW_types & no_CpG$ERV_rel_rate > median_mut_rate)$MST

WS_low_mut <- subset(no_CpG, substr(no_CpG$MST, 1, 3) %in% WS_types & no_CpG$ERV_rel_rate <= median_mut_rate)$MST
WS_high_mut <- subset(no_CpG, substr(no_CpG$MST, 1, 3) %in% WS_types & no_CpG$ERV_rel_rate > median_mut_rate)$MST

neutral_low_mut <- subset(no_CpG, substr(no_CpG$MST, 1, 3) %in% neutral_types & no_CpG$ERV_rel_rate <= median_mut_rate)$MST
neutral_high_mut <- subset(no_CpG, substr(no_CpG$MST, 1, 3) %in% neutral_types & no_CpG$ERV_rel_rate > median_mut_rate)$MST

top_windows_df$CpG_TpG <- rowSums(top_windows_df[,CpG_TpG])
top_windows_df$SW_low_mut <- rowSums(top_windows_df[,SW_low_mut])
top_windows_df$SW_high_mut <- rowSums(top_windows_df[,SW_high_mut])
top_windows_df$WS_low_mut <- rowSums(top_windows_df[,WS_low_mut])
top_windows_df$WS_high_mut <- rowSums(top_windows_df[,WS_high_mut])
top_windows_df$neutral_low_mut <- rowSums(top_windows_df[,neutral_low_mut])
top_windows_df$neutral_high_mut <- rowSums(top_windows_df[,neutral_high_mut])
head(top_windows_df)

CpG_TpG_quantiles <- as.data.frame(tapply(top_windows_df$CpG_TpG, top_windows_df$tajD_group, mean))
SW_low_mut_quantiles <- as.data.frame(tapply(top_windows_df$SW_low_mut, top_windows_df$tajD_group, mean))
SW_high_mut_quantiles <- as.data.frame(tapply(top_windows_df$SW_high_mut, top_windows_df$tajD_group, mean))
WS_low_mut_quantiles <- as.data.frame(tapply(top_windows_df$WS_low_mut, top_windows_df$tajD_group, mean))
WS_high_mut_quantiles <- as.data.frame(tapply(top_windows_df$WS_high_mut, top_windows_df$tajD_group, mean))
neutral_low_mut_quantiles <- as.data.frame(tapply(top_windows_df$neutral_low_mut, top_windows_df$tajD_group, mean))
neutral_high_mut_quantiles <- as.data.frame(tapply(top_windows_df$neutral_high_mut, top_windows_df$tajD_group, mean))

all_categories <- cbind(CpG_TpG_quantiles, SW_low_mut_quantiles, SW_high_mut_quantiles,
                        WS_low_mut_quantiles, WS_high_mut_quantiles, neutral_low_mut_quantiles, neutral_high_mut_quantiles)
all_categories$D_quantile <- as.factor(rownames(all_categories))
all_categories$D_quantile_rank <- rank(all_categories$D_quantile)
colnames(all_categories) <- c("CpG_TpG", "SW_low", "SW_high", "WS_low", "WS_high", "Neutral_low", "Neutral_high", "D_quantiles", "Quantile_rank")

library(tidyr)
library(dplyr)
categories_long <- gather(all_categories, group, num_abundant, CpG_TpG:Neutral_high, factor_key=TRUE)
categories_long$gBGC <- sub("\\_.*", "", categories_long$group)
categories_long$mut_rate <- sub(".*\\_", "", categories_long$group)
categories_long$gBGC <- ifelse(categories_long$gBGC == 'CpG', 'CpG TpG', categories_long$gBGC)
categories_long$mut_rate <- ifelse(categories_long$mut_rate == 'TpG', 'CpG TpG', categories_long$mut_rate)

#Figure 4 in paper
ggplot(data=categories_long, aes(x=D_quantiles, y=num_abundant, fill=gBGC, pattern=mut_rate)) +
  geom_bar_pattern(position='stack', stat='identity',
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025) + 
  scale_pattern_manual("Mutation Rate", values = c(high = "stripe", low = "none", `CpG TpG`='circle'),
                       labels=c("CpG TpG", "High", "Low")) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  xlab("D Quantiles") + ylab("Average # Abundant Subtypes") +
  ggtitle("Average # Abundant Subtypes by D Quantiles across Subtype Categories\n*Artificial Windows from Permuted Data") + 
  scale_x_discrete(limits = rev) + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#Make proportions plot suggested by Jean 
all_categories$total_abun <- rowSums(all_categories[,1:7])
all_categories_prop <- all_categories[,1:7]/rowSums(all_categories[,1:7])
all_categories_prop$tajD_group <- rownames(all_categories)
all_categories_prop$D_quantile_rank <- rank(desc(all_categories_prop$tajD_group))

cat_only <- all_categories_prop[,1:7]
(cat_only[nrow(cat_only),] - cat_only[1,]) / cat_only[1,]

prop_long <- gather(all_categories_prop, cat_group, prop_abundant, CpG_TpG:Neutral_high, factor_key=TRUE)

ggplot(data=prop_long, aes(x=reorder(tajD_group, desc(tajD_group)), y=prop_abundant, color=cat_group)) + 
  geom_line(aes(group=cat_group)) +
  geom_point() +
  xlab("D Quantiles") + ylab("Proportion of Abundant Subtypes") + 
  ggtitle("Proportion of Abundant by D Quantiles across Subtype Categories\n*Artificial Windows from Permuted Data") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_discrete(name = "Category Group")

all_categories$total_abun[c(1,20)]
(all_categories$total_abun[10] + all_categories$total_abun[11])/2

abundant_normalized <- as.data.frame(sweep(as.matrix(all_categories[,1:7]), 2, c(4,7,21,16,16,25,7), `/`))
