####################################################
### FUNCTIONS FOR FERNANDEZ-FOURNIER ET AL. 2021 ###
####################################################

# Title: Do we need to identify adaptive genetic variation when prioritizing populations for conservation?
# Authors: Philippe Fernandez-Fournier1*, Jayme M. M. Lewthwaite1*, Arne Ø. Mooers1
# 1Department of Biological Sciences, Simon Fraser University, 8888 University Drive, Burnaby, British Columbia V5A 1S6, Canada.
# *PFF and JL contributed equally to this work and should be considered joint first author.
# Journal: Conservation Genetics


#### function changing ordered shapley into ordered pops and their ranks ####
shap2ranks_pine <- function(shap){
  rankos <- 1:length(shap)
  names(rankos) <- names(shap)
  rankos <- rankos[order(as.numeric(names(rankos)))]
  names(rankos) <- paste("pop", names(rankos), sep="")
  return(rankos)
}


shap2ranks <- function(shap){
  rankos <- vector()
  for(i in 1:length(shap)){
    rankos[as.numeric(names(shap))[i]] <- i
    names(rankos)[as.numeric(names(shap))[i]] <- paste("pop", names(shap)[i], sep="")
  }
  return(rankos)
}

#### Calculate Shapley index for NJ or NeighborNet object (phangorn package) ####
shapley_index <- function(obj){
  require(phangorn)
  if(length(grep("phylo", class(obj)))==0){
    stop("argument needs to be a Neighbor Joining or NeighborNet object - function NJ() or neighborNet()")
  }
  splits.n <- as.splits(obj)
  splits.m <- as.matrix(splits.n)
  shap <- vector()
  for(i in 1:length(splits.m[1,])){
    sx <- vector() #size of a split set containing the taxon 
    sbarx <- vector() #size of the complementary set that does NOT contain the taxon 
    X <- length(splits.m[1,]) #total number of taxa
    wgt <- attr(splits.n, "weights") #split weights
    wgt[wgt<0] <- 0 #Set negative branch lengths to zero

    for(j in 1:length(splits.m[,1])){
      sx[j] <- length(which(splits.m[j,]==splits.m[j,i]))
      sbarx[j] <- length(which(splits.m[j,]!=splits.m[j,i]))
    }
    shap[i] <- sum(abs(sbarx)/(X*abs(sx))*(wgt)) #Equation from Volkmann et al. 2014, PLOS One
  }
  names(shap) <- colnames(splits.m)
  return(shap[order(shap,decreasing=TRUE)])
}

# Compare distance matrix between a PCA using a set of SNPs and the overall (Input is set of snps)
compare_d_rk <- function(signif.snps){
  
  reso <- list()
  
  snps.gds <- read.gdsn(index.gdsn(YWAR.gds, "snp.rs.id"))
  snps.sel <- read.gdsn(index.gdsn(YWAR.gds, "snp.id"))[snps.gds %in% signif.snps]
  
  pca.adapt <- snpgdsPCA(YWAR.gds, autosome.only=FALSE, bayesian = TRUE, verbose = FALSE, missing.rate = 0.1, 
                         sample.id = read.gdsn(index.gdsn(YWAR.gds, "sample.id"))[miss.rate.samp<0.1],
                         snp.id = snps.sel) # retains all axes
  
  # Choose significant PCs (Tracy-Widom distribution test)
  tw.adapt <- tw(na.exclude(pca.adapt$eigenval), length(na.exclude(pca.adapt$eigenval)), criticalpoint = 0.9793)
  (n.pc.adapt <- tw.adapt$SigntEigenL)
  
  #make data set of centroids
  pca.adapt.cen <- aggregate(pca.adapt$eigenvect, list(pop = pops.sr), mean)
  d.pca.adapt <- stats::dist(pca.adapt.cen[,2:(n.pc.adapt+1)], method = "euclidean")
  
  # NeighborNet and rank population distinctiveness
  nnet.pca.adapt <- neighborNet(d.pca.adapt)
  rk.pca.adapt <- shap2ranks(shapley_index(nnet.pca.adapt))
  
  reso$obs.dist.mantel <- mantel(d.pca.all, d.pca.adapt)$statistic
  
  reso$obs.rank.kendall <- cor(rk.pca.all, rk.pca.adapt, method = "kendall")
  
  return(reso)
}

# create a NeighborNet on SplitsTRee via R
splitstree_pff <- function(dist.file, nexus.out = NULL, splitstree.exe = "E:/Programs/SplitsTree/SplitsTree.exe") {
  # generate an appropriate file name, if none is provided
  if(missing(nexus.out)) {
    nexus.out <- 'splitstree-output.nex'
  }
  # plotting commands to be passed to splitstree
  splitstree_script <- paste0(
    # "begin SplitsTree;\n",
    "LOAD FILE='",  dist.file, "'\n",
    "EXPORT FILE='",  nexus.out, "' REPLACE=yes\n",
    "EXECUTE FILE='",  nexus.out, "'\n",
    "EXPORT FILE='",  nexus.out, "' REPLACE=yes\n",
    "QUIT\n" # , "end;"
  )
  writeLines(splitstree_script, "cmd_file_temp.txt")
  # run splitstree
  system2(splitstree.exe, args = c("-S --commandLineMode", paste0("--commandFile ", file.path(getwd()), "/cmd_file_temp.txt")))   
  invisible(file.remove("cmd_file_temp.txt"))
}


#Function for pairwise Jost D for vcfR object
pairwise_JOST <- function(vcf, pops){
  require(spaa)
  pops <- as.factor(pops)
  diffo <- pairwise_genetic_diff(vcf, pops, method = "jost") #long
  pw_j <- colMeans(diffo[,grep("Dest_Chao",colnames(diffo))], na.rm = TRUE)
  poink <- strsplit(names(pw_j), split = "_")
  pop1<-vector()
  pop2<-vector()
  for(i in 1:length(poink)){
    pop1[i]<- poink[[i]][3]
    pop2[i]<- poink[[i]][4]
  }
  df<- data.frame(distcols = pop1, distrows = pop2, dist = pw_j)
  jostd.m <- list2dist(df)
  jostd.m[which(jostd.m<0)] <- min(jostd.m[jostd.m>0])/100 #change negative value to minimum D divided by 100...
  return(jostd.m)
}


