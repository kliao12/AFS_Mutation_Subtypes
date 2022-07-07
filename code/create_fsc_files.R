### Script to simulate neutral SFS For each mutation subtype ###
library(data.table)

#Import mutation subtype AFS to get total # of sites 
afs_files <- list.files("/net/wonderland/home/ksliao/zoellner_research/mst_afs/", pattern = "\\.txt$")

for(file in afs_files){
    data <- fread(paste0("/net/wonderland/home/ksliao/zoellner_research/mst_afs/", file))
    mst <- substr(file, 1,7)
    num_sites <- sum(data$V1)

    filename <- paste0("/net/wonderland/home/ksliao/fsc26_linux64/mst_par_files/", mst, ".par")

    cat("//Number of population samples (demes)\n", file = filename, append=FALSE)
    cat("1\n", file = filename, append=TRUE)
    cat("//Population effective sizes (number of genes)\n", file = filename, append=TRUE)
    cat("20000\n", file = filename, append=TRUE)
    cat("//Sample sizes\n", file = filename, append=TRUE)
    cat("7112\n", file = filename, append=TRUE)
    cat("//Growth rates  : negative growth implies population expansion\n", file = filename, append=TRUE)
    cat("0\n", file = filename, append=TRUE)
    cat("//Number of migration matrices : 0 implies no migration between demes\n", file = filename, append=TRUE)
    cat("0\n", file = filename, append=TRUE)
    cat("//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix\n", file = filename, append=TRUE)
    cat("0  historical event\n", file = filename, append=TRUE)
    cat("//Number of independent loci [chromosome]\n", file = filename, append=TRUE)
    cat("1 0\n", file = filename, append=TRUE)
    cat("//Per chromosome: Number of linkage blocks\n", file = filename, append=TRUE)
    cat("1\n", file = filename, append=TRUE)
    cat("//per Block: data type, num loci, rec. rate and mut rate + optional parameters\n", file = filename, append=TRUE)
    cat(paste0("DNA ", num_sites, " 0 0.00000005"),file = filename, append=TRUE)
}
