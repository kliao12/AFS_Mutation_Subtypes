##### Title: Script to run the main analysis in paper 
##### Date: 7/6/22

#Import libraries
library(ggplot2)
library(gridExtra)
library(scales)
library(MASS)
library(data.table)
library(gee)
library(ggpubr)

# Use chr22 file to easily get MST vector
allele_counts <- read.table("/Users/kevinliao/Desktop/Michigan Research/Mutation Rates/chr22_fullanno_to_output (2020_02_29 01_13_00 UTC).vcf")
#allele_counts <- read.table("C:\\Users\\Kevin Liao\\Desktop\\Michigan Research\\Mutation Rates\\chr22_fullanno_to_output.vcf")
colnames(allele_counts) <- c("CHR", "REF","ALT", "AA","AC", "AC_correct", "ANNO", "MT", "KMER", "MST")



##### Begin creating dataframe with summary stats for each MST: Tajima D, F*, Ratio of singletons to doubletons, Mutation Rates, etc #####


# Compute Tajimas D using genome wide SFS
taj_d_list = list()

for(i in unique(allele_counts$MST)){
  if(i == '.'){next}
  MST_sfs <-read.table(paste0("/Users/kevinliao/Desktop/Michigan Research/AFS Paper Analysis/Genomewide SFS/",i, ' (2020_02_29 01_13_00 UTC).txt'), header = FALSE)
  colnames(MST_sfs) <- c("AC_correct", "site_num")
  MST_sfs <- subset(MST_sfs, site_num != 0)
  
  #Compute pairwise diff
  MST_sfs$n <- 3556*2
  MST_sfs$pi <- MST_sfs$AC_correct*MST_sfs$site_num*(MST_sfs$n - MST_sfs$site_num)/
    choose(3556*2, 2)
  sum_pi <- sum(MST_sfs$pi) 
  
  #Compute numerator: Pi - S/a1
  S <- sum(MST_sfs$AC_correct)
  a1 <- sum(1/seq(1,3555*2))
  a2 <- sum(1/seq(1,3555*2)^2)
  n <- 7112
  
  numer <- sum_pi - S/a1
  
  #Compute denomenator
  e1 <- (1/a1)*((n+1)/(3*(n-1)) - 1/a1)
  e2 <- 1/(a1^2 + a2)*((2*(n^2 + n + 3)/(9*n*(n-1))) - (n+2)/(n*a1) + (a2/a1^2))
  
  denom <- sqrt(e1*S + e2*S*(S-1))

  taj_d_list[[i]] <- numer/denom
  
}

#Create Taj D dataframe from list
taj_d_df <- as.data.frame(cbind(names(taj_d_list), unlist(taj_d_list, use.names=FALSE)))
colnames(taj_d_df) <- c("MST", "D")
taj_d_df$MST <- as.character(taj_d_df$MST)
taj_d_df$D <- as.numeric(as.vector(taj_d_df$D))


# Compute F* using genome wide SFS
f_list = list()

for(i in unique(allele_counts$MST)){
  if(i == '.'){next}
  
  MST_sfs <-read.table(paste0("/Users/kevinliao/Desktop/Michigan Research/AFS Paper Analysis/Genomewide SFS/",i, ' (2020_02_29 01_13_00 UTC).txt'), header = FALSE)
  colnames(MST_sfs) <- c("AC_correct", "site_num")
  MST_sfs <- subset(MST_sfs, site_num != 0)
  
  #Compute numerator: MPD - singletons
  MST_sfs$n <- 7112
  MST_sfs$pi <- MST_sfs$AC_correct*MST_sfs$site_num*(MST_sfs$n - MST_sfs$site_num)/
    choose(7112, 2)
  sum_pi <- sum(MST_sfs$pi) 
  
  eps_1 <- MST_sfs[1, "AC_correct"]
  eps_n_1 <- MST_sfs[nrow(MST_sfs), "AC_correct"]
  delta_1_n_1 <- 0
  n <- 7112
  eta_1 <- (eps_1 + eps_n_1) / (1 + delta_1_n_1)
  
  f_numer <- sum_pi - ((n-1)/n)*eta_1
  
  #Compute denomenator: Use formula from Fu and Li's paper
  eta <- sum(MST_sfs$AC_correct)
  a_n <- sum(1/seq(1,n-1))
  a_n_1 <- sum(1/seq(1,n))
  b_n <- sum(1/seq(1,n-1)^2)
  c_n <- 2*(n*a_n - 2*(n-1))/((n-1)*(n-2))
  d_n <- c_n + (n-2)/(n-1)^2 + (2/(n-1))*(1.5 - (2*a_n_1 - 3)/(n-2) - 1/n)
  
  v_f <- (d_n + 2*(n^2+n+3)/(9*n*(n-1)) - 2*(4*b_n - 6 + 8/n)/(n-1)) / (a_n^2+b_n)
  u_f <- (n/(n+1) + (n+1)/(3*(n-1)) - 4/(n*(n-1)) + (2*(n+1)/(n-1)^2)*(a_n_1 - 2*n/(n+1)))/a_n - v_f
  
  var_f <- u_f*eta + v_f*eta^2
  
  f_denom <- sqrt(var_f)
  
  f_list[[i]] <- f_numer / f_denom
}

#Create F* dataframe from list
f_stat_df <- as.data.frame(cbind(names(f_list), unlist(f_list, use.names=FALSE)))
colnames(f_stat_df) <- c("MST", "F")
f_stat_df$MST <- as.character(f_stat_df$MST)
f_stat_df$F <- as.numeric(as.vector(f_stat_df$F))


# Compute proportion of singletons to doubletons
prop_list = list()

for(i in unique(allele_counts$MST)){
  if(i == '.'){next}
  
  MST_sfs <-read.table(paste0("/Users/kevinliao/Desktop/Michigan Research/AFS Paper Analysis/Genomewide SFS/",i, ' (2020_02_29 01_13_00 UTC).txt'), header = FALSE)
  colnames(MST_sfs) <- c("AC_correct", "site_num")
  MST_sfs <- subset(MST_sfs, site_num != 0)
  
  prop_list[[i]] <- MST_sfs[1, "AC_correct"]/MST_sfs[2, "AC_correct"]
}

singles_list = list()
doubles_list = list()
triples_list = list() 

for(i in unique(allele_counts$MST)){
  if(i == '.'){next}
  
  MST_sfs <-read.table(paste0("/Users/kevinliao/Desktop/Michigan Research/AFS Paper Analysis/Genomewide SFS/",i, ' (2020_02_29 01_13_00 UTC).txt'), header = FALSE)
  colnames(MST_sfs) <- c("AC_correct", "site_num")
  MST_sfs <- subset(MST_sfs, site_num != 0)
  singles_list[[i]] <- prop.table(MST_sfs$AC_correct)[1]
  doubles_list[[i]] <- prop.table(MST_sfs$AC_correct)[2]
  triples_list[[i]] <- prop.table(MST_sfs$AC_correct)[3]
}

props_df <- as.data.frame(cbind(names(singles_list), unlist(singles_list, use.names=FALSE)))
props_df <- cbind(props_df, unlist(doubles_list, use.names=FALSE), unlist(triples_list, use.names=FALSE))
colnames(props_df) <- c("MST", "singles", "doubles","triples")

#Create prop dataframe from list
S_D_df <- as.data.frame(cbind(names(prop_list), unlist(prop_list, use.names=FALSE)))
colnames(S_D_df) <- c("MST", "S_D_ratio")
S_D_df$MST <- as.character(S_D_df$MST)
S_D_df$S_D_ratio <- as.numeric(as.vector(S_D_df$S_D_ratio))

singles <- as.data.frame(cbind(names(singles_list), unlist(singles_list, use.names=FALSE)))
colnames(singles) <- c("MST", "prop_singles")
singles$MST <- as.character(singles$MST)
singles$prop_singles <- as.numeric(as.vector(singles$prop_singles))

#Import mutation rates
mut_rates <-read.table("/Users/kevinliao/Desktop/Michigan Research/Mutation Rates/Jeds subtype mutation rates (2020_02_29 01_13_00 UTC).csv", sep=',', header = TRUE, stringsAsFactors=FALSE)

mut_rates$MST <- ifelse(substr(mut_rates$Type,1,1) == 'A',
                        paste0(substr(mut_rates$Type,1,1), substr(mut_rates$Type,3,4),'.', substr(mut_rates$Motif,1,3)),
                        paste0(substr(mut_rates$Type,2,3), substr(mut_rates$Type,5,5),'.', substr(mut_rates$Motif,1,3))
)

#Create combined dataset
merged_data <- merge(taj_d_df, f_stat_df, by='MST')
merged_data <- merge(merged_data, S_D_df, by='MST')
merged_data <- merge(merged_data, mut_rates, by='MST')
#remove(allele_counts)

#Subset subtypes to without CpGs
no_CpG <- merged_data[-c(83,87,91,95), ]

CpGs <- c("C_T.ACG","C_T.CCG","C_T.GCG","C_T.TCG")
merged_data$CpG_fill <- ifelse(merged_data$MST %in% CpGs, 2, 1)

merged_data$mt_grouping <- with(merged_data, ifelse(substr(MST,1,3) == 'A_C', 1,0))
merged_data$mt_grouping <- with(merged_data, ifelse(substr(MST,1,3) == 'A_T', 2, merged_data$mt_grouping))
merged_data$mt_grouping <- with(merged_data, ifelse(substr(MST,1,3) == 'C_G', 3, merged_data$mt_grouping))




####### Run analysis. Also create graphs and figures #########


#0.5) Create plot for Tajimas D by MST
taj_d <- ggplot(data=merged_data, aes(x=MST, y=D, fill = as.factor(mt_grouping)) ) +
  geom_bar(stat="identity", show.legend = FALSE) + 
  xlab("Mutation Subtype") +
  ylab("D") +
  theme(axis.text.x=element_text(angle=90)) +
  ggtitle("Tajimas D by Mutation Subtype") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits=c(-2.3, -1.4),oob = rescale_none) +
  scale_fill_manual(values = c("3"="cornflowerblue", "2" = "cornflowerblue",
                               "1" = "cornflowerblue", "0" = "grey70"))
                     
taj_d


# 1) Relationship between mutation rate and ratio of singleton to doubletons by MST
cor.test(no_CpG$ERV_rel_rate, no_CpG$S_D_ratio)
cor.test(merged_data$ERV_rel_rate, merged_data$S_D_ratio)

merged_data$MT <- substr(merged_data$MST, 1, 3)
no_CpG$MT <- substr(no_CpG$MST, 1, 3)

a1 <- ggplot(merged_data, aes(x=ERV_rel_rate, y=S_D_ratio)) +
  geom_point(aes(color=as.factor(CpG_fill))) + 
  xlab("Mutation Rate") + ylab("Singletons / Doubletons") + geom_smooth(method = "lm", se = FALSE) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits=c(3.5,8),oob = rescale_none) +
  scale_color_manual(values = c("2" = "orangered2", "1" = "black"),
                    name="Mutation\nSubtype",
                    breaks = c(1,2),
                    labels =c("Other","CpG TpG")) 
a1

a2 <- ggplot(no_CpG, aes(x=ERV_rel_rate, y=S_D_ratio)) + 
  geom_point(aes(color = "black")) +
  xlab("ERV Mutation Rate") +  ylab("Singletons / Doubletons") + geom_smooth(method = "lm", se = FALSE) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits=c(5,8),oob = rescale_none) + scale_color_manual(values = 1, name = "Mutation\nSubtype", labels = "Non CpG TpG")
a2

grid.arrange(a1, a2)

no_CpG$MT_arrow <- paste0(substr(no_CpG$MT,1,1), "->", substr(no_CpG$MT,3,3))
merged_data$MT_arrow <- paste0(substr(merged_data$MT,1,1), "->", substr(merged_data$MT,3,3))

breaks_fun <- function(x) {
  if (max(x) > 24) {
    seq(0, 120, 24)
  } else {
    seq(0, 8, 1)
  }
}

equal_breaks <- function(n = 4, s = 0.05, ...){
  function(x){
    # rescaling
    d <- s * diff(range(x)) / (1+2*s)
    seq(min(x)+d, max(x)-d, length=n)
  }
}

all_six <- ggplot(data=merged_data, aes(x=ERV_rel_rate, y=S_D_ratio)) + geom_point(aes(color=as.factor(CpG_fill))) + geom_smooth(method = "lm", se=FALSE) + 
  facet_wrap(~MT_arrow, scales = "free_x") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Mutation Rate") + ylab("Singletons / Doubletons") + 
  scale_color_manual(values = c("2" = "orangered2", "1" = "black"),
                     name="Mutation\nSubtype",
                     breaks = c(1,2),
                     labels =c("Other","CpG TpG")) +
  theme(panel.spacing = unit(2, "lines")) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.0001), breaks=equal_breaks(n=4))
all_six
ggarrange(a1, all_six, ncol=1,labels=c("a)", "b)"))

cor.test(subset(merged_data, substr(merged_data$MST, 1, 3) == 'A_C')$ERV_rel_rate, subset(merged_data, substr(merged_data$MST, 1, 3) == 'A_C')$S_D_ratio)
cor.test(subset(merged_data, substr(merged_data$MST, 1, 3) == 'A_G')$ERV_rel_rate, subset(merged_data, substr(merged_data$MST, 1, 3) == 'A_G')$S_D_ratio)
cor.test(subset(merged_data, substr(merged_data$MST, 1, 3) == 'A_T')$ERV_rel_rate, subset(merged_data, substr(merged_data$MST, 1, 3) == 'A_T')$S_D_ratio)
cor.test(subset(merged_data, substr(merged_data$MST, 1, 3) == 'C_A')$ERV_rel_rate, subset(merged_data, substr(merged_data$MST, 1, 3) == 'C_A')$S_D_ratio)
cor.test(subset(merged_data, substr(merged_data$MST, 1, 3) == 'C_G')$ERV_rel_rate, subset(merged_data, substr(merged_data$MST, 1, 3) == 'C_G')$S_D_ratio)
cor.test(subset(merged_data, substr(merged_data$MST, 1, 3) == 'C_T')$ERV_rel_rate, subset(merged_data, substr(merged_data$MST, 1, 3) == 'C_T')$S_D_ratio)


# 2): Biased gene conversion (gBGC) on the AFS. Split into mutation types and WS, SW, neutral 
A_C <- subset(merged_data, substr(MST,1,3) == 'A_C')
non_A_C <- subset(merged_data, substr(MST,1,3) != 'A_C')

A_G <- subset(merged_data, substr(MST,1,3) == 'A_G')
non_A_G <- subset(merged_data, substr(MST,1,3) != 'A_G')

A_T <- subset(merged_data, substr(MST,1,3) == 'A_T')
non_A_T <- subset(merged_data, substr(MST,1,3) != 'A_T')

C_A <- subset(merged_data, substr(MST,1,3) == 'C_A')
non_C_A <- subset(merged_data, substr(MST,1,3) != 'C_A')

C_G <- subset(merged_data, substr(MST,1,3) == 'C_G')
non_C_G <- subset(merged_data, substr(MST,1,3) != 'C_G')

C_T <- subset(merged_data, substr(MST,1,3) == 'C_T')
non_C_T <- subset(merged_data, substr(MST,1,3) != 'C_T')

C_T_noCpG <- subset(no_CpG, substr(MST,1,3) == 'C_T')
non_C_T_noCpG<- subset(no_CpG, substr(MST,1,3) != 'C_T')

wilcox.test(as.numeric(as.vector(A_C$F)), as.numeric(as.vector(non_A_C$F)), paired=FALSE)
wilcox.test(as.numeric(as.vector(A_C$D)), as.numeric(as.vector(non_A_C$D)), paired=FALSE)

#significant
wilcox.test(as.numeric(as.vector(A_G$F)), as.numeric(as.vector(non_A_G$F)), paired=FALSE)
wilcox.test(as.numeric(as.vector(A_G$D)), as.numeric(as.vector(non_A_G$D)), paired=FALSE)

wilcox.test(as.numeric(as.vector(A_T$F)), as.numeric(as.vector(non_A_T$F)), paired=FALSE)
wilcox.test(as.numeric(as.vector(A_T$D)), as.numeric(as.vector(non_A_T$D)), paired=FALSE)

#significant
wilcox.test(as.numeric(as.vector(C_A$F)), as.numeric(as.vector(non_C_A$F)), paired=FALSE)
wilcox.test(as.numeric(as.vector(C_A$D)), as.numeric(as.vector(non_C_A$D)), paired=FALSE)

wilcox.test(as.numeric(as.vector(C_G$F)), as.numeric(as.vector(non_C_G$F)), paired=FALSE)
wilcox.test(as.numeric(as.vector(C_G$D)), as.numeric(as.vector(non_C_G$D)), paired=FALSE)

#significant
wilcox.test(as.numeric(as.vector(C_T$F)), as.numeric(as.vector(non_C_T$F)), paired=FALSE)
#not significant after taking out CpGs
wilcox.test(as.numeric(as.vector(C_T_noCpG$F)), as.numeric(as.vector(non_C_T_noCpG$F)), paired=FALSE)

###Checking by WS SW Neutral bins
#Try with all mutation subtypes
WS <- subset(merged_data, substr(MST,1,3) %in% c("A_C", "A_G"))
non_WS <- subset(merged_data, !(substr(MST,1,3) %in% c("A_C", "A_G")))
wilcox.test(as.numeric(as.vector(WS$D)), as.numeric(as.vector(non_WS$D)), paired=FALSE)
mean(WS$D)

SW <- subset(merged_data, substr(MST,1,3) %in% c("C_A", "C_T"))
non_SW <- subset(merged_data, !(substr(MST,1,3) %in% c("C_A", "C_T")))
wilcox.test(as.numeric(as.vector(SW$D)), as.numeric(as.vector(non_SW$D)), paired=FALSE)
mean(SW$D)

neutral <- subset(merged_data, substr(MST,1,3) %in% c("C_G","A_T"))
non_neutral <- subset(merged_data, !(substr(MST,1,3) %in% c("C_G","A_T")))
wilcox.test(as.numeric(as.vector(neutral$D)), as.numeric(as.vector(non_neutral$D)), paired=FALSE)
mean(neutral$D)

### Remove CpGs like that one tishkoff paper
remove_CpG <- subset(merged_data, substr(merged_data$MST,6,7) != 'CG')

WS <- subset(remove_CpG, substr(MST,1,3) %in% c("A_C", "A_G"))
non_WS <- subset(remove_CpG, !(substr(MST,1,3) %in% c("A_C", "A_G")))
wilcox.test(as.numeric(as.vector(WS$D)), as.numeric(as.vector(non_WS$D)), paired=FALSE)
t.test(WS$D, non_WS$D)

SW <- subset(remove_CpG, substr(MST,1,3) %in% c("C_A", "C_T"))
non_SW <- subset(remove_CpG, !(substr(MST,1,3) %in% c("C_A", "C_T")))
wilcox.test(as.numeric(as.vector(SW$D)), as.numeric(as.vector(non_SW$D)), paired=FALSE)
t.test(SW$D, non_SW$D)

neutral <- subset(remove_CpG, substr(MST,1,3) %in% c("C_G","A_T"))
non_neutral <- subset(remove_CpG, !(substr(MST,1,3) %in% c("C_G","A_T")))
wilcox.test(as.numeric(as.vector(neutral$D)), as.numeric(as.vector(non_neutral$D)), paired=FALSE)
t.test(neutral$D, non_neutral$D)

A_C <- subset(remove_CpG, substr(MST,1,3) == 'A_C')
non_A_C <- subset(remove_CpG, substr(MST,1,3) != 'A_C')

A_G <- subset(remove_CpG, substr(MST,1,3) == 'A_G')
non_A_G <- subset(remove_CpG, substr(MST,1,3) != 'A_G')

A_T <- subset(remove_CpG, substr(MST,1,3) == 'A_T')
non_A_T <- subset(remove_CpG, substr(MST,1,3) != 'A_T')

C_A <- subset(remove_CpG, substr(MST,1,3) == 'C_A')
non_C_A <- subset(remove_CpG, substr(MST,1,3) != 'C_A')

C_G <- subset(remove_CpG, substr(MST,1,3) == 'C_G')
non_C_G <- subset(remove_CpG, substr(MST,1,3) != 'C_G')

C_T <- subset(remove_CpG, substr(MST,1,3) == 'C_T')
non_C_T <- subset(remove_CpG, substr(MST,1,3) != 'C_T')

wilcox.test(as.numeric(as.vector(A_C$D)), as.numeric(as.vector(non_A_C$D)), paired=FALSE)

#significant
wilcox.test(as.numeric(as.vector(A_G$D)), as.numeric(as.vector(non_A_G$D)), paired=FALSE)

wilcox.test(as.numeric(as.vector(A_T$D)), as.numeric(as.vector(non_A_T$D)), paired=FALSE)

#significant
wilcox.test(as.numeric(as.vector(C_A$D)), as.numeric(as.vector(non_C_A$D)), paired=FALSE)

wilcox.test(as.numeric(as.vector(C_G$D)), as.numeric(as.vector(non_C_G$D)), paired=FALSE)

wilcox.test(as.numeric(as.vector(C_T$D)), as.numeric(as.vector(non_C_T$D)), paired=FALSE)


#3) Effect of AFS heterogeneity on Demographic inference
library(ggplot2)
library(gridExtra)
library(scales)
library(MASS)

#Read in exponential growth demographic inference
data <-read.table("/Users/kevinliao/Desktop/Michigan Research/Demographic Inference/combined_exp (2021_08_24 20_19_58 UTC).txt", header = FALSE)
colnames(data) <- c("theta", "LogLik","nu", "T", "MST", "theta_w", "theta_pi")

data1 <- merge(data, merged_data, by="MST")
data1$theta_dadi <- data1$theta/data1$nMotifs
data1$theta_w_norm <- data1$theta_w/data1$nMotifs 
data1$theta_pi_norm <- data1$theta_pi/data1$nMotifs 

#Compare theta estimates from Dadi, watterson, MPD
theta_dadi3 <- ggplot(data=data1, aes(x=MST, y=theta_dadi)) +
  geom_bar(stat="identity") + 
  xlab("Mutation Subtype") +
  ylab("Theta - Dadi") +
  theme(axis.text.x=element_text(angle=90)) +
  theme(plot.title = element_text(hjust = 0.5)) 

theta_dadi3

theta_w3 <- ggplot(data=data1, aes(x=MST, y=theta_w_norm)) +
  geom_bar(stat="identity") + 
  xlab("Mutation Subtype") +
  ylab("Theta - Watterson") +
  theme(axis.text.x=element_text(angle=90)) +
  theme(plot.title = element_text(hjust = 0.5)) 

theta_w3

theta_pi3 <- ggplot(data=data1, aes(x=MST, y=theta_pi_norm)) +
  geom_bar(stat="identity") + 
  xlab("Mutation Subtype") +
  ylab("Theta - MPD") +
  theme(axis.text.x=element_text(angle=90)) +
  #ggtitle("Theta Standardized by Mutation Subtype Genomewide") +
  theme(plot.title = element_text(hjust = 0.5)) 

theta_pi3

grid.arrange(theta_dadi3, theta_w3, theta_pi3, nrow = 3, ncol = 1)

##### Compute absolute mutation rate b/c previous paper used relative rates
lambda <- 60/sum(data1$nMotifs*data1$mut_rate)
data1$abs_mutRate <- data1$mut_rate*lambda

data1$N_ref <- data1$theta/(4*data1$nMotifs*data1$abs_mutRate)
data1$T_gen <- 2*data1$N_ref * data1$T
data1$growth_rate <- (data1$nu)^(1/data1$T_gen) - 1
#data1$growth_rate <- log(data1$nu)/data1$T_gen

#Export to file 
library(openxlsx)
write.xlsx(data1, "Desktop/Michigan Research/AFS Paper Analysis/Figures in Paper/exp_parameters.xlsx")

#Make plots comparing Nu and T
T_plot <- ggplot(data=data1, aes(x=MST, y=T)) +
  geom_bar(stat="identity") + 
  xlab("Mutation Subtype") + 
  ylab("Time since ancestral population started growing") +
  theme(axis.text.x=element_text(angle=90)) +
  #ggtitle("Theta Standardized by Mutation Subtype Genomewide") +
  theme(plot.title = element_text(hjust = 0.5)) 
T_plot

nu_plot <- ggplot(data=data1, aes(x=MST, y=nu)) +
  geom_bar(stat="identity") + 
  xlab("Mutation Subtype") +
  ylab("Current N_eff relative to ancestral N_eff") +
  theme(axis.text.x=element_text(angle=90)) +
  #ggtitle("Theta Standardized by Mutation Subtype Genomewide") +
  theme(plot.title = element_text(hjust = 0.5))  + 
  ylim(0,275)
nu_plot

grid.arrange(T_plot, nu_plot, nrow = 2)

data2 <- merge(data1, singles, by='MST')

ggplot(data=data2, aes(x=prop_singles, y=nu)) +
  geom_point(stat="identity") + 
  xlab("Proportion of Singletons") +
  ylab("Current Pop Size / Ancestral Pop Size") +
  theme(axis.text.x=element_text(angle=90)) +
  ggtitle("Population Growth by Proportion of Singletons for 96 MSTs") +
  theme(plot.title = element_text(hjust = 0.5))  

cor.test(data2$nu, data2$prop_singles)

ggplot(data=subset(data1, mut_rate < 0.06), aes(x=mut_rate, y=nu)) +
  geom_point(stat="identity") + 
  xlab("Mutation Rate") +
  ylab("Present N_eff / Ancestral N_eff") +
  theme(axis.text.x=element_text(angle=90)) +
  #ggtitle("Theta Standardized by Mutation Subtype Genomewide") +
  theme(plot.title = element_text(hjust = 0.5))  

cor.test(subset(data1, mut_rate < 0.06)$mut_rate, subset(data1, mut_rate < 0.06)$nu)


N_ref_plot <- ggplot(data=data1, aes(x=MST, y=N_ref)) +
  geom_bar(stat="identity") + 
  xlab("Mutation Subtype") +
  ylab("Effective Population Size") +
  theme(axis.text.x=element_text(angle=90)) +
  ggtitle("Ancestral Effective Population Size" ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=16,face="bold"),
        plot.title = element_text(size = 18, face = "bold"))
N_ref_plot

Tgen_plot <- ggplot(data=data1, aes(x=MST, y=T_gen)) +
  geom_bar(stat="identity") + 
  xlab("Mutation Subtype") +
  ylab("Time (generations)") +
  theme(axis.text.x=element_text(angle=90)) +
  ggtitle("Time Since Ancestral Population Started Growing") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=16,face="bold"),
        plot.title = element_text(size = 18, face = "bold"))
Tgen_plot

grid.arrange(N_ref_plot, Tgen_plot, nrow = 2)

gr_plot <- ggplot(data=data1, aes(x=MST, y=growth_rate)) +
  geom_bar(stat="identity") + 
  xlab("Mutation Subtype") +
  ylab("Growth Rate per Generation") +
  theme(axis.text.x=element_text(angle=90)) +
  #ggtitle("Theta Standardized by Mutation Subtype Genomewide") +
  theme(plot.title = element_text(hjust = 0.5)) 
gr_plot


#Look at output for bottlneck model
bottleneck <- read.table("/Users/kevinliao/Desktop/Michigan Research/Demographic Inference/bottleneck_7_1_21 (2021_08_24 20_19_58 UTC).txt")
colnames(bottleneck) <- c("theta", "max_ll", "nuB","nuF","TB","TF","MST")
bottleneck$time <- bottleneck$TB + bottleneck$TF

bottleneck1 <- merge(bottleneck, mut_rates, by="MST")
bottleneck2 <- merge(bottleneck1, singles, by='MST')
bottleneck3 <- merge(bottleneck2, taj_d_df, by='MST')

lambda <- 60/sum(bottleneck3$nMotifs*bottleneck3$ERV_rel_rate)
bottleneck3$abs_mutRate <- bottleneck3$ERV_rel_rate*lambda
bottleneck3$N_ref <- bottleneck3$theta/(4*bottleneck3$nMotifs*bottleneck3$abs_mutRate)
bottleneck3$T_gen <- 2*bottleneck3$N_ref * bottleneck3$time

#export for paper
write.xlsx(bottleneck3, "Desktop/Michigan Research/AFS Paper Analysis/Figures in Paper/BN_parameters.xlsx")

time <- ggplot(data=bottleneck, aes(x=MST, y=time)) +
  geom_bar(stat="identity") + 
  xlab("Mutation Subtype") +
  theme(axis.text.x=element_text(angle=90)) +
  theme(plot.title = element_text(hjust = 0.5)) 

nuB <- ggplot(data=bottleneck, aes(x=MST, y=nuB)) +
  geom_bar(stat="identity") + 
  xlab("Mutation Subtype") +
  theme(axis.text.x=element_text(angle=90)) +
  theme(plot.title = element_text(hjust = 0.5)) 

nuF <- ggplot(data=bottleneck, aes(x=MST, y=nuF)) +
  geom_bar(stat="identity") + 
  xlab("Mutation Subtype") +
  theme(axis.text.x=element_text(angle=90)) +
  theme(plot.title = element_text(hjust = 0.5)) 

grid.arrange(time, nuB, nuF)

cor.test(bottleneck3$nuB, bottleneck3$D)
cor.test(subset(bottleneck3, bottleneck3$D < -1.6)$nuB, subset(bottleneck3, bottleneck3$D < -1.6)$D)

bottleneck3$CpG <- as.factor(ifelse(substr(bottleneck3$MST, 1,3) == 'C_T' & substr(bottleneck3$MST,6,7) == 'CG', 2, 1))

cor.test(bottleneck3$nuB, bottleneck3$TF)
ggplot(data=bottleneck3, aes(x=TF, y=nuB)) + geom_point() 

jpeg("C:\\Users\\Kevin Liao\\Desktop\\Michigan Research\\AFS Paper Analysis\\Figures in Paper\\bottleneck.jpg", units = 'in', width = 7.5, height = 5.5, res = 600)
ggplot(data=bottleneck3, aes(x=D, y=nuB)) + geom_point() +
  ggtitle("Population Contraction During Bottleneck by Tajima's D for 96 MSTs") +
  xlab("Tajima's D") + ylab("Population Contraction During Bottleneck") +
  geom_smooth(method='lm') + theme(plot.title = element_text(hjust = 0.5, size = 13))
dev.off()

#4 Effect of AFS heterogeneity on selection. Look at distribution of subtypes and Tajima's D across the genome 
### Look at distribution of subtypes vs Tajima's D
source(Abundant_windows.R)

### Make figure for abundant windows by gBGC, mut rate across Tajima's D quantiles 
source(Make_figure_abundant_categories.R)

### Look at increase in false positives for windows abundant in a given subtype 
single <- "C_G.GCG"
single_MST <- read.table(paste0("/Users/kevinliao/Desktop/Michigan Research/Selection/tajD_MST_corrected/tajD_", single, " (2020_02_29 01_13_00 UTC).txt"), header=TRUE)

fig_3a <- ggplot(single_MST, aes(x=prop_weird, y=taj_D)) +
  geom_point() +
  xlab("Proportion of C_G.GCG") +
  ylab("Tajimas D") +
  ggtitle("Tajimas D by Proportion of C_G.GCG in 100kb Windows") +
  theme(plot.title = element_text(hjust = 0.5)) +
  #annotate("text", x=0.03, y=-0.5, label=" cor = -0.08, p < 2.2e-16 ", size = 4.5) +
  geom_smooth(method = "lm", se = FALSE) 

fig_3a

lower_D <- quantile(single_MST$taj_D, probs = 0.05)
single_MST$class <- "Overall"

cutoff = .95
single_MST$prop_weird_sig <- ifelse(single_MST$prop_weird > quantile(single_MST$prop_weird, probs = cutoff),
                                    1, 0)
high <- subset(single_MST, prop_weird_sig == 1)
high$class <- "High"

df_95 <- rbind(single_MST, high)

fig_3b <- ggplot(df_95, aes(x=taj_D, fill= class)) +
  geom_density(alpha=0.25) +
  scale_x_continuous(limits=c(-3,-0.25)) +
  labs(title = "Distribution of D for Windows Abundant in G[C->G]G", x= "Tajima's D",
       y = "Density") +
  scale_fill_discrete(name="Window Type", labels=c(paste0("High ","G[C->G]G" ), "All Windows")) +
  geom_vline(aes(xintercept = lower_D)) +
  theme(plot.title = element_text(hjust = 0.5))
fig_3b

grid.arrange(fig_3a, fig_3b)

#Using top 10% in figure for paper
cutoff = .90
single_MST$prop_weird_sig <- ifelse(single_MST$prop_weird > quantile(single_MST$prop_weird, probs = cutoff),
                                    1, 0)
high <- subset(single_MST, prop_weird_sig == 1)
high$class <- "High"

df_90 <- rbind(single_MST, high)

ggplot(df_90, aes(x=taj_D, fill= class)) +
  geom_density(alpha=0.25) +
  scale_x_continuous(limits=c(-3,-0.25)) +
  labs(title = "Distribution of D for Windows Abundant in G[C->G]G", x= "Tajima's D",
       y = "Density") +
  scale_fill_discrete(name="Window Type", labels=c(paste0("High ","G[C->G]G" ), "All Windows")) +
  geom_vline(aes(xintercept = lower_D)) +
  theme(plot.title = element_text(hjust = 0.5))

### False positive rates
### Run loop to get all FPR increases
MST_fdr <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(MST_fdr) <- c("MST","corr","corr_p", "low","high","FDR_increase")

for(i in unique(allele_counts$MST)){
  print(i)
  if(i == '.'){next}
  
  mst_data <- read.table(paste0("/Users/kevinliao/Desktop/Michigan Research/Selection/tajD_MST_corrected/tajD_", i, " (2020_02_29 01_13_00 UTC).txt"), header=TRUE)
  #write a loop to get FDR by cutoff for "high" proportion of weird sites
  FDR_df <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(FDR_df) <- c("quantile", "prop_cutoff" ,"FDR")

  mst_data$tajD_sig <- ifelse(mst_data$taj_D < quantile(mst_data$taj_D, probs = 0.05), 1, 0)
  
  for(cutoff in seq(0.05,0.90, by = 0.01)){
    if(cutoff > 0.10 & cutoff < 0.85){next}
    mst_data$prop_weird_sig <- ifelse(mst_data$prop_weird > quantile(mst_data$prop_weird, probs = cutoff),
                                      1, 0)
    
    table_D_prop <- table(mst_data$tajD_sig, mst_data$prop_weird_sig)
    
    FDR <- table_D_prop[2,2]/(table_D_prop[1,2] + table_D_prop[2,2])
    
    FDR_df[nrow(FDR_df)+1, ] <- c(cutoff, quantile(mst_data$prop_weird, probs = cutoff) , FDR) 
  }
  
  fdr_summary <- summary(FDR_df$FDR)
  FDR_increase <- FDR_df$FDR[length(FDR_df$FDR)] - 0.05

  corr <- cor.test(mst_data$prop_weird, mst_data$taj_D)

  MST_fdr[nrow(MST_fdr)+1, ] <- c(i, corr$estimate, corr$p.value, 
                                  as.numeric(fdr_summary[1]), as.numeric(fdr_summary[6]), as.numeric(FDR_increase))
}

MST_fdr$corr <- as.numeric(MST_fdr$corr)
MST_fdr$FDR_increase <- as.numeric(MST_fdr$FDR_increase)
MST_fdr$corr_p <- as.numeric(MST_fdr$corr_p)
MST_fdr$fill_p <- ifelse(MST_fdr$corr_p > 0.05, 1, 0)
MST_fdr$high <- as.numeric(MST_fdr$high)

MST_fdr$CpG <- ifelse(substr(MST_fdr$MST, 1,3) == 'C_T' & substr(MST_fdr$MST,6,7) == 'CG', 1,0)

ggplot(data=MST_fdr, aes(x=reorder(MST, as.numeric(high)), y=as.numeric(high), fill=as.factor(CpG))) +
  geom_bar(stat="identity") + xlab("Mutation Subtype") + ylab("FPR for Windows in Top 10% Abundance") +
  theme(axis.text.x=element_text(angle=90)) +
  ggtitle("FPR Increase by Mutation Subtype") +
  theme(plot.title = element_text(hjust = 0.5, size=22), legend.position="bottom", 
        axis.text=element_text(size=12)) +
  scale_y_continuous(limits=c(0.04, 0.12),oob = rescale_none) +
  guides(fill=guide_legend(title="")) + 
  scale_fill_manual(values=c("#999999", "firebrick4"), labels =c("Non CpG", "CpG")) +
  geom_hline(yintercept = 0.05)








