### Script to run newD_neutral_sim.R in parallel b/c takes too damn long 
args = commandArgs(trailingOnly=TRUE)
subtype = as.character(args[1])
print(subtype)

library(data.table)
library(ggplot2)

# Step 1) Import genomewide sites. Need sites rather than SFS b/c will be bootstrapping 
gw_sites <- fread("/net/wonderland/home/ksliao/zoellner_research/demographic_inf/data/window_VCFs/genWide_window.txt")
segregating <- subset(gw_sites, gw_sites$V3 != 0) #54,721,662 segregating sites

#Function compute traditional D 
compute_oldD <- function(mst_afs){
  n <- 3556*2; S <- sum(mst_afs$AC_correct); a1 <- sum(1/seq(1,n-1))
  a2 <- sum(1/seq(1,3555*2)^2)

  theta_s <- S/a1
  mst_afs$pi <- mst_afs$AC_correct*mst_afs$site_num*(n - mst_afs$site_num)/choose(n, 2)
  sum_pi <- sum(mst_afs$pi)


  D <- sum_pi - theta_s
  e1 <- (1/a1)*((n+1)/(3*(n-1)) - 1/a1)
  e2 <- 1/(a1^2 + a2)*((2*(n^2 + n + 3)/(9*n*(n-1))) - (n+2)/(n*a1) + (a2/a1^2))

  denom <- sqrt(e1*S + e2*S*(S-1))
  return(denom)
}

#Function to compute D-2 statistic
compute_newD <- function(mst_afs){
  n <- 3556*2; S <- sum(mst_afs$AC_correct); a1 <- sum(1/seq(1,n-1))
  singles <- mst_afs$AC_correct[1]; doubles <- mst_afs$AC_correct[2]
  theta_s <- (S-singles-doubles)/(a1 - 3/2)

  #mst_afs$pi <- mst_afs$AC_correct*mst_afs$site_num*(n - mst_afs$site_num)/choose(n, 2)
  mst_afs$pi <- mst_afs$AC_correct*mst_afs$site_num*(n - mst_afs$site_num)
  sum_pi <- sum(mst_afs$pi)

  #theta_w <- ((n-1)/(choose(n,2)*(n-2))) * (sum_pi - singles*(n-1) - doubles*(n-2)*2)
  #theta_w <- ((1-(2/(n*(n-1))))^-1) * 1/choose(n,2) * (sum_pi - singles*(n-1) - doubles*(n-2)*2)
  theta_w <- ((1-(2/n)-(2*(n-2))/(n*(n-1)))^-1) * 1/choose(n,2) * (sum_pi - singles*(n-1) - doubles*(n-2)*2)
  newD_num <- theta_w - theta_s
  return(newD_num)
}

newD_df <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(newD_df) <- c("MST","newD")

#Subset sites to subtype
mst <- subset(segregating, segregating$V4 == subtype)

#Get # of sites to export as well   
num_mst_sites <- nrow(mst)

# Step 3) Make AFS from sites and compute numerator of D statistic without singletons and doubletons 
afs <- as.data.frame(table(mst$V3))
colnames(afs) <- c("site_num", "AC_correct")
afs$site_num <- as.numeric(afs$site_num); afs$AC_correct <- as.numeric(afs$AC_correct)

numerator <- compute_newD(afs)

# Step 4) Use neutral simulations to get variance
sim_data <- fread(paste0("/net/wonderland/home/ksliao/fsc26_linux64/", subtype, "/", subtype, "_DAFpop0.obs"), skip=2)
sim_data2 <- sim_data[,c(2:7113)]


#get variance in denominator to use in next loop
sim_newD <- vector()
for(i in 1:nrow(sim_data2)){
  temp_AFS <- as.data.frame(cbind(seq(1,7112), t(sim_data2[i,])))
  colnames(temp_AFS) <- c("site_num", "AC_correct")
  D_2 <- compute_newD(temp_AFS)
  sim_newD <- c(sim_newD, D_2)
}

mst_sd <- sd(sim_newD)

#Computing new D for 10,000 null simulated AFS. Reuse SD from 10,000 simulations in denominator of each one
sim_newD2 <- vector()
for(i in 1:nrow(sim_data2)){
  temp_AFS <- as.data.frame(cbind(seq(1,7112), t(sim_data2[i,])))
  colnames(temp_AFS) <- c("site_num", "AC_correct")
  D_2 <- compute_newD(temp_AFS)/mst_sd
  sim_newD2 <- c(sim_newD2, D_2)
}

summary(sim_newD2)
table(sim_newD2 > 2)
tab <- table(abs(sim_newD2) > 1.96)
t1E <- tab[2]/(tab[1]+tab[2])

m <- as.data.frame(matrix(0, ncol = 3, nrow = 1))
m[1,] <- c(subtype, num_mst_sites, t1E)

write.table(m, paste0("~/zoellner_research/AFS_project/output/newD_neutral_sim/", subtype,".txt"), col.names=FALSE, row.names=FALSE, quote=FALSE)
