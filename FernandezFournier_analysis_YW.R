#################################################
### R CODE FOR FERNANDEZ-FOURNIER ET AL. 2021 ###
#################################################

# Title: Do we need to identify adaptive genetic variation when prioritizing populations for conservation?
# Authors: Philippe Fernandez-Fournier1*, Jayme M. M. Lewthwaite1*, Arne Ø. Mooers1
# 1Department of Biological Sciences, Simon Fraser University, 8888 University Drive, Burnaby, British Columbia V5A 1S6, Canada.
# *PFF and JL contributed equally to this work and should be considered joint first author.
# Journal: Conservation Genetics

library(SNPRelate)
library(AssocTests)
library(ggplot2)
library(adegenet)
library(phangorn)
library(vegan)
library(rnaturalearth)
library(png)
library(geosphere)
library(spaa)

setwd("E:/sfuvault/SFU/MooersLab/PopGen/Analysis/Final R code - Cons Bio")

# Load functions
source("popgen_functions.R")

#Yellow Warbler data
YWAR.meta <- read.delim("YWAR_RAD_meta.txt")
pop.coord <- read.csv("pop_coord.csv")
YWAR.gds <- snpgdsOpen("ywar.gds", allow.duplicate = TRUE)
YWAR.gds$root

#Environmental SNPs
cands.all <- read.csv("pvalues_all_BIOCLIM.csv", row.names = 1)

#### IDENTIFY SIGNIFICANT SNPS ####
# p-values dataframe replaced by 0/1 for if p-val is  significant or not
'
signif.matrix <- matrix(data = 0, nrow = dim(cands.all)[1], ncol = (dim(cands.all)[2]-1), 
                        dimnames = list(rownames(cands.all), colnames(cands.all)[-length(colnames(cands.all))]))
for(i in 1:length(signif.matrix[1,])){
  biovar <- cands.all[,i]
  bh.i <- rank(biovar)
  bh.m <- length(biovar)
  bh.Q <- 0.001
  bh.corr <- (bh.i/bh.m)*0.001
  signifs <- which(biovar < bh.corr)
  signif.matrix[signifs,i] <- 1
}

# Make it binary (1 = significant SNP, 0 = not significant)
significant.snp <- vector()
for(i in 1:length(signif.matrix[,1])){
  if(sum(signif.matrix[i,])>0){
    significant.snp[i] <- 1
  } else {
    significant.snp[i] <- 0
  }
}

cands.all$significant.snp <- significant.snp #Already done and saved in cands.all

# How many significant SNPs per bioclim variable?
as.data.frame(colSums(signif.matrix))
'
#### Explore ####

#Proportion of missing SNPs per sample
miss.rate.samp <- apply(read.gdsn(index.gdsn(YWAR.gds, "genotype")), 1, 
                        function(x) length(which(x==3))/length(x)) 

hist(miss.rate.samp, breaks = 50, main = "Proportion of missing SNPs per sample")

# Proportion of missing alleles per SNP
miss.rate.snp <- apply(read.gdsn(index.gdsn(YWAR.gds, "genotype")), 2, 
                       function(x) length(which(x==3))/length(x)) 

hist(miss.rate.snp, breaks = 30,  main = "Proportion of missing alleles per SNP")

# Number of retained samples and SNPs
length(which(miss.rate.samp < 0.1))
length(which(miss.rate.snp  < 0.1))


#### PCA OF ALL SNPS ####

#Select for SNPs with <10% missing rate and Samples with <10% missing SNPs
pca.all <- snpgdsPCA(YWAR.gds, autosome.only = FALSE, bayesian = TRUE, missing.rate = 0.1, 
                        sample.id = read.gdsn(index.gdsn(YWAR.gds, "sample.id"))[miss.rate.samp<0.1]) # retains all axes

# Vector of populations for each individual
pops.sr <- vector()
for(i in 1:length(pca.all$sample.id)){
  pops.sr[i] <- YWAR.meta$Pop[YWAR.meta$Field_Number==pca.all$sample.id[i]]
}

# Choose significant PCs (Tracy-Widom distribution test)
tw.all <- tw(na.exclude(pca.all$eigenval), length(na.exclude(pca.all$eigenval)), criticalpoint = 0.9793)
(n.pc.all <- tw.all$SigntEigenL) # number of significant principal components retained

# PCA plot
pca.all.gg <- as.data.frame(pca.all$eigenvect)
pca.all.gg$group <- as.character(pops.sr)

p1<-ggplot(pca.all.gg,aes(x=V1,y=V2,color=group)) +
  geom_point(size=2) +
  scale_color_manual(values=funky(21)) +
  xlab(paste("PC1 (", round(pca.all$varprop[1]*100, 2), "%)", sep = "")) +
  ylab(paste("PC2 (", round(pca.all$varprop[2]*100, 2), "%)", sep = "")) +
  xlim(c(-0.2, 0.2)) +
  ylim(c(-0.2, 0.2)) +
  theme_bw() +
  theme(legend.position="none")
p1

# Contribution of PCs
pc.percent.all <- as.vector(na.exclude(pca.all$varprop*100))
barplot(pc.percent.all, 
        col = c(rep("grey50", n.pc.all), rep("grey90", length(pc.percent.all)-n.pc.all)),
        xlab = "Principal Components (dark grey = significant)",
        ylab = "Percent contribution")

#make data set of centroids
pca.all.cen <- aggregate(pca.all$eigenvect, list(pop = pops.sr), mean)
d.pca.all <- stats::dist(pca.all.cen[,2:(n.pc.all+1)], method = "euclidean") 

# NeighborNet and rank population distinctiveness
nnet.pca.all <- neighborNet(d.pca.all)
(rk.pca.all <- shap2ranks(shapley_index(nnet.pca.all)))

#Plot NNET and color tips by rank
colfunc <- colorRampPalette(c("red","yellow","green"))
ord.all <- rk.pca.all[rank(as.numeric(nnet.pca.all$tip.label))]
tipcols.all <- colfunc(dim(as.matrix(d.pca.all))[1])[ord.all]

plot(nnet.pca.all, tip.color = tipcols.all, type = '2D')


#### PCA OF ADAPTIVE SNPS ####
signif.snps <- row.names(cands.all)[which(cands.all$significant.snp==1)]
snps.gds <- read.gdsn(index.gdsn(YWAR.gds, "snp.rs.id"))
snps.select <- read.gdsn(index.gdsn(YWAR.gds, "snp.id"))[snps.gds %in% signif.snps]

pca.adapt <- snpgdsPCA(YWAR.gds, autosome.only=FALSE, bayesian = TRUE, verbose = FALSE, missing.rate = 0.1, 
                        sample.id = read.gdsn(index.gdsn(YWAR.gds, "sample.id"))[miss.rate.samp<0.1],
                        snp.id = snps.select) # retains all axes

# Choose significant PCs (Tracy-Widom distribution test)
tw.adapt <- tw(na.exclude(pca.adapt$eigenval), length(na.exclude(pca.adapt$eigenval)), criticalpoint = 0.9793)
(n.pc.adapt <- tw.adapt$SigntEigenL)

## PCA plot
pca.adapt.gg <- as.data.frame(pca.adapt$eigenvect)
pca.adapt.gg$group <- as.character(pops.sr)

p2<-ggplot(pca.adapt.gg,aes(x=V1,y=-V2,color=group)) + # V2 flipped horizontally
  geom_point(size=2) +
  scale_color_manual(values=funky(21)) +
  xlab(paste("PC1 (", round(pca.all$varprop[1]*100, 2), "%)", sep = "")) +
  ylab(paste("PC2 (", round(pca.all$varprop[2]*100, 2), "%)", sep = "")) +
  xlim(c(-0.2, 0.2)) +
  ylim(c(-0.2, 0.2)) +
  theme_bw() +
  theme(legend.position="none")
p2 

# Contribution of PCs
pc.percent.adapt <- as.vector(na.exclude(pca.adapt$varprop*100))
barplot(pc.percent.adapt, 
        col = c(rep("grey50", n.pc.adapt), rep("grey90", length(pc.percent.adapt)-n.pc.adapt)),
        xlab = "Principal Components (dark grey = significant)",
        ylab = "Percent contribution")

#make data set of centroids
pca.adapt.cen <- aggregate(pca.adapt$eigenvect, list(pop = pops.sr), mean)
d.pca.adapt <- stats::dist(pca.adapt.cen[,2:(n.pc.adapt+1)], method = "euclidean") #How many axes do we retain?????

# NeighborNet and rank population distinctiveness
nnet.pca.adapt <- neighborNet(d.pca.adapt)
rk.pca.adapt <- shap2ranks(shapley_index(nnet.pca.adapt))

#Plot NNET and color tips by rank
colfunc <- colorRampPalette(c("red","yellow","green"))
ord.clim <- rk.pca.adapt[rank(as.numeric(nnet.pca.adapt$tip.label))]
tipcols.clim <- colfunc(dim(as.matrix(d.pca.adapt))[1])[ord.clim]

plot(nnet.pca.adapt, tip.color = tipcols.clim, type = '2D')


#### COMPARE ALL vs. ADAPTIVE SNPS ####

# Distance matrices
(obs.dist.mantel <- mantel(d.pca.all, d.pca.adapt)$statistic)

# Ranks
rk.pca.all
rk.pca.adapt
(obs.rank.kendall <- cor(rk.pca.all, rk.pca.adapt, method = "kendall"))

(ranks.table <- data.frame(ranks.all = rk.pca.all,
                          ranks.clim = rk.pca.adapt))

#### COMPARE TO RANDOM SUBSETS ####
# Draw many different sets of (any) SNPs of same size as adaptive SNPs

comp.dist <- vector()
comp.rank <- vector()
(n.snps <- length(which(cands.all$significant.snp==1)))

for(i in 1:100){ # < 5min for 1k iterations
  cat("\r", paste("Neutral SNP iteration number:", i, "   "))
  set.seed(i+12)
  
  ## PCA of SUBSET SNPs ##
  sampled.snps <- sample(rownames(cands.all), n.snps, replace = FALSE) #SAMPLING OF SNPS HERE
  snps.gds <- read.gdsn(index.gdsn(YWAR.gds, "snp.rs.id"))
  snps.select <- read.gdsn(index.gdsn(YWAR.gds, "snp.id"))[snps.gds %in% sampled.snps]
  
  pca.sub <- snpgdsPCA(YWAR.gds, autosome.only=FALSE, bayesian =  TRUE, verbose = FALSE, missing.rate = 0.1, 
                       sample.id = read.gdsn(index.gdsn(YWAR.gds, "sample.id"))[miss.rate.samp<0.1],
                       snp.id = snps.select) # retains all axes
  
  pops.sub <- vector()
  for(j in 1:length(pca.sub$sample.id)){
    pops.sub[j] <- YWAR.meta$Pop[YWAR.meta$Field_Number == pca.sub$sample.id[j]]
  }
  
  ## CHOOSE SIGNIFICANT PRINCIPAL COMPONENTS
  #Tracy-Widom distribution test
  tw.sub <- tw(na.exclude(pca.sub$eigenval), length(na.exclude(pca.sub$eigenval)), criticalpoint = 0.9793)
  (n.pc.sub <- tw.sub$SigntEigenL)
  
  #make data set of centroids
  pca.sub.cen <- aggregate(pca.sub$eigenvect, list(pop = pops.sub), mean)
  d.pca.sub <- stats::dist(pca.sub.cen[,2:(n.pc.sub+1)], method = "euclidean") #How many axes do we retain?????
  
  rk.pca.sub <- shap2ranks(shapley_index(neighborNet(d.pca.sub)))
  
  #nnet.pca.sub <- neighborNet(d.pca.sub)
  #colfunc <- colorRampPalette(c("red","yellow","green"))
  #ord.sub <- rk.pca.sub[rank(as.numeric(nnet.pca.sub$tip.label))]
  #tipcols.sub <- colfunc(dim(as.matrix(d.pca.sub))[1])[ord.sub]
  #plot(nnet.pca.sub, tip.color = tipcols.sub, type = '2D')

  # compare dist matrices
  mtest.i <- mantel(d.pca.all, d.pca.sub)
  comp.dist[i] <- mtest.i$statistic
  
  # compare ranks
  comp.rank[i] <- cor(rk.pca.all, rk.pca.sub, method = "kendall")
  
}

comp.subs.pca <- data.frame(iteration = c(1:length(comp.rank)),
                           compare.dist = comp.dist,
                           compare.ranks = comp.rank)
head(comp.subs.pca)

# Compare distance matrices
(pval <- length(which(comp.subs.pca$compare.dist < obs.dist.mantel))/length(comp.subs.pca$compare.dist))

g1 <- ggplot(comp.subs.pca, aes(compare.dist)) +
  geom_histogram(binwidth = 0.01, color = "grey15", fill = "grey15") +
  geom_vline(xintercept = obs.dist.mantel, color = "red", size = 2) +
  xlab("Mantel correlation to distance matrix using all SNPs")+
  geom_text(x = min(comp.subs.pca$compare.dist)+((range(comp.subs.pca$compare.dist)[2]-range(comp.subs.pca$compare.dist)[1])/20),
            y = 25, label=paste("p =", pval), size = 6, fontface = "italic", color = "grey15")+
  theme_classic() 
g1

# Compare ranks
(pval2 <- length(which(comp.subs.pca$compare.ranks < obs.rank.kendall))/length(comp.subs.pca$compare.ranks))

g2 <- ggplot(comp.subs.pca, aes(compare.ranks)) +
  geom_histogram(binwidth = 0.02, color = "grey15", fill = "grey15") +
  geom_vline(xintercept = obs.rank.kendall, color = "red", size = 2) +
  xlab("Kendall Tau correlation to ranks using all SNPs")+
  geom_text(x = min(comp.subs.pca$compare.ranks)+((range(comp.subs.pca$compare.ranks)[2]-range(comp.subs.pca$compare.ranks)[1])/20),
            y = 25, label=paste("p =", pval2), size = 6, fontface = "italic", color = "grey15")+
  theme_classic()
g2


#### EFFECT OF SIGNIFICANCE FDR CUTOFF ####
# Compare to random subsets while changing number of significant SNPs
# Draw many different sets of (any) SNPs of same size as adaptive SNPs

# Include most significant SNPs with different cutoffs
cutoff.Q = c(0.0001, 0.00001, 0.000001) #vector of multiple cutoffs (default = 0.001)

snp.list <- list()

for(j in 1:length(cutoff.Q)){
  signif.matrix <- matrix(data = 0, nrow = dim(cands.all)[1], ncol = (dim(cands.all)[2]-1), 
                          dimnames = list(rownames(cands.all), colnames(cands.all)[-length(colnames(cands.all))]))
  
  for(i in 1:length(signif.matrix[1,])){
    biovar <- cands.all[,i]
    bh.i <- rank(biovar)
    bh.m <- length(biovar)
    bh.Q <- cutoff.Q[j]
    bh.corr <- (bh.i/bh.m)*bh.Q
    signifs <- which(biovar < bh.corr)
    signif.matrix[signifs,i] <- 1
  }
  
  # Make it binary (1 = significant SNP, 0 = not significant)
  signif.bin <- vector()
  for(i in 1:length(signif.matrix[,1])){
    if(sum(signif.matrix[i,])>0){
      signif.bin[i] <- 1
    } else {
      signif.bin[i] <- 0
    }
  }
  snp.list[[j]] <- row.names(cands.all)[signif.bin==1]
  names(snp.list)[j] <- paste("Q", cutoff.Q[j])
}

snp.list
names(snp.list)
lapply(snp.list, length)


# Algorithm
n.iter = 1000

comp.subs.nsnp <- data.frame(iteration = numeric(), 
                             compare.dist = numeric(), 
                             compare.ranks = numeric(), 
                             number.of.snps = numeric())

for(j in 1:length(snp.list)){ #number of different significant snps you want
  
  comp.dist.nsnp <- vector()
  comp.rank.nsnp <- vector()

  for(i in 1:n.iter){ # < 5min for 1k iterations
    print(paste(names(snp.list[j]), "-> Running iteration number", i))
    
    set.seed(i+32)
    n.snps <- length(snp.list[[j]])
    sampled.snps <- sample(rownames(cands.all), n.snps, replace = FALSE) #SAMPLING OF SNPS HERE
    snps.gds <- read.gdsn(index.gdsn(YWAR.gds, "snp.rs.id"))

    pca.sub <- snpgdsPCA(YWAR.gds, autosome.only=FALSE, bayesian =  TRUE, verbose = FALSE, missing.rate = 0.1, 
                         sample.id = read.gdsn(index.gdsn(YWAR.gds, "sample.id"))[miss.rate.samp<0.1],
                         snp.id = read.gdsn(index.gdsn(YWAR.gds, "snp.id"))[snps.gds %in% sampled.snps]) # retains all axes
    
    pops.sub <- vector()
    for(k in 1:length(pca.sub$sample.id)){
      pops.sub[k] <- YWAR.meta$Pop[YWAR.meta$Field_Number==pca.sub$sample.id[k]]
    }
    
    ## CHOOSE SIGNIFICANT PRINCIPAL COMPONENTS
    #Tracy-Widom distribution test
    tw.sub <- tw(na.exclude(pca.sub$eigenval), length(na.exclude(pca.sub$eigenval)), criticalpoint = 0.9793)
    n.pc.sub <- tw.sub$SigntEigenL
    
    #make data set of centroids
    pca.sub.cen <- aggregate(pca.sub$eigenvect, list(pop = pops.sub), mean)
    d.pca.sub <- stats::dist(pca.sub.cen[,2:(n.pc.sub+1)], method = "euclidean")
    
    rk.pca.sub <- shap2ranks(shapley_index(neighborNet(d.pca.sub)))
    
    # compare dist matrices
    mtest.i <- mantel(d.pca.all, d.pca.sub)
    comp.dist.nsnp[i] <- mtest.i$statistic
    
    # compare ranks
    comp.rank.nsnp[i] <- cor(rk.pca.all, rk.pca.sub, method = "kendall")
    
  }
  comp.subs.nsnp <- rbind(comp.subs.nsnp, 
                     data.frame(iteration = c(1:length(comp.rank.nsnp)),
                                compare.dist = comp.dist.nsnp,
                                compare.ranks = comp.rank.nsnp,
                                number.of.snps = rep(n.snps, length(comp.rank.nsnp))))
  
}

comp.subs.nsnp


# P-values
comp.subs.nsnp$number.of.snps <- factor(comp.subs.nsnp$number.of.snps)

pvals.ranks <- vector()
for(i in 1:length(snp.list)){
  obs.rank.kendall.i <- compare_d_rk(snp.list[[i]])$obs.rank.kendall
  pvals.ranks[i] <- length(which(comp.subs.nsnp$compare.ranks[comp.subs.nsnp$number.of.snps==length(snp.list[[i]])] < obs.rank.kendall.i))/
    length(comp.subs.nsnp$compare.ranks[comp.subs.nsnp$number.of.snps==length(snp.list[[i]])])
  names(pvals.ranks)[i] <- length(snp.list[[i]])
}
#names(pvals.ranks) <- names(snp.list)
pvals.ranks

#### MAPS ####
topo <- ne_download(scale = 10, type = 'GRAY_HR_SR_W', category = 'raster', returnclass = "sf")

colfunc2 <- grey.colors(5, start = 1, end = 0.4)

yw.png <- readPNG("E:/sfuvault/SFU/MooersLab/PopGen/Graphs/yellowwarbler.png")

#Map of Ranks using ALL SNPs
raster::plot(topo, xlim=c(-152,-50), ylim=c(15,75),
             col = colfunc2, xlab = "Longitude", ylab = "Latitude", legend = FALSE)
maps::map(xlim=c(-152,-50), ylim=c(15,75), fill = FALSE, lwd = 1, add = TRUE) 
points(x=pop.coord$lon, y = pop.coord$lat, pch = 20, cex = 3, 
       col = colfunc(length(pop.coord$pop))[rk.pca.all])
text(x = pop.coord$lon, y = pop.coord$lat, labels = pop.coord$pop, cex = 0.5)
text(x = -133, y = 16.8, labels = "All genetic variation", cex = 0.7)
rasterImage(yw.png, -143, 18, -123, 29)


#Map of Ranks using CLIMATE SNPs
raster::plot(topo, xlim=c(-152,-50), ylim=c(15,75),
             col = colfunc2, xlab = "Longitude", ylab = "Latitude", legend = FALSE)
maps::map(xlim=c(-152,-50), ylim=c(15,75), fill = FALSE, lwd = 1, add = TRUE) 
points(x=pop.coord$lon, y = pop.coord$lat, pch = 20, cex = 3, 
       col = colfunc(length(pop.coord$pop))[rk.pca.adapt])
text(x = pop.coord$lon, y = pop.coord$lat, labels = pop.coord$pop, cex = 0.5)
text(x = -133, y = 16.8, labels = "Adaptive genetic variation", cex = 0.7)
rasterImage(yw.png, -143, 18, -123, 29)




#### APPENDIX ####

#### Overall Fst ####

# Overall Fst
fst.all <- snpgdsFst(YWAR.gds, as.factor(pops.sr), method = "W&C84", autosome.only = FALSE, missing.rate = 0.1, 
                     sample.id = read.gdsn(index.gdsn(YWAR.gds, "sample.id"))[miss.rate.samp<0.1]) 
fst.all$Fst
fst.all$MeanFst


#### Pairwise Fst ####
library(BEDASSLE)

signif.snps <- row.names(cands.all)[which(cands.all$significant.snp==1)]
snps.gds <- read.gdsn(index.gdsn(YWAR.gds, "snp.rs.id"))
snps.select <- read.gdsn(index.gdsn(YWAR.gds, "snp.id"))[snps.gds %in% signif.snps]

YWAR.gds$root
geno <- read.gdsn(index.gdsn(YWAR.gds, "genotype"))
colnames(geno) <- read.gdsn(index.gdsn(YWAR.gds, "snp.rs.id"))
rownames(geno) <- read.gdsn(index.gdsn(YWAR.gds, "sample.id"))
geno[geno == 3] <- NA # 3 means missing data

pop.fst <- vector()
for(i in 1:length(rownames(geno))){
  pop.fst[i] <- YWAR.meta$Pop[YWAR.meta$Field_Number==rownames(geno)[i]]
}
pop.fst <- factor(pop.fst)

# allele count
geno.pop <- matrix(0L, nrow = length(levels(pop.fst)), ncol = dim(geno)[2])
for(i in 1:length(levels(pop.fst))){
  geno.pop[i,] <- apply(geno[which(pop.fst == levels(pop.fst)[i]),], 2, function(x) sum(x, na.rm=TRUE))
}
colnames(geno.pop) <- read.gdsn(index.gdsn(YWAR.gds, "snp.rs.id"))
rownames(geno.pop) <- as.character(levels(pop.fst))

# chromosome sample size
chrom.pop <- matrix(0L, nrow = length(levels(pop.fst)), ncol = dim(geno)[2])
for(i in 1:length(levels(pop.fst))){
  chrom.pop[i,] <- apply(geno[which(pop.fst == levels(pop.fst)[i]),], 2, function(x) length(which(!is.na(x)))*2)
}
dimnames(chrom.pop) <- dimnames(geno.pop)

# distances
d.Fst.all <- as.dist(calculate.all.pairwise.Fst(geno.pop, chrom.pop))

d.Fst.adapt <- as.dist(calculate.all.pairwise.Fst(geno.pop[,snps.gds %in% signif.snps], 
                                                  chrom.pop[,snps.gds %in% signif.snps]))
# ranks
rk.Fst.all   <- shap2ranks(shapley_index(neighborNet(d.Fst.all)))
rk.Fst.adapt <- shap2ranks(shapley_index(neighborNet(d.Fst.adapt)))

(obs.dist.mantel.fst <- mantel(d.Fst.all, d.Fst.adapt)$statistic)

(obs.rank.kendall.fst <- cor(rk.Fst.all, rk.Fst.adapt, method = "kendall"))

# NeighborNet and rank population distinctiveness
nnet.Fst.all <- neighborNet(d.Fst.all)
colfunc <- colorRampPalette(c("red","yellow","green"))
ord.Fst <- rk.Fst.all[rank(as.numeric(nnet.Fst.all$tip.label))]
tipcols.Fst <- colfunc(dim(as.matrix(d.Fst.all))[1])[ord.Fst]
plot(nnet.Fst.all, tip.color = tipcols.Fst, type = '2D')

nnet.Fst.adapt <- neighborNet(d.Fst.adapt)
colfunc <- colorRampPalette(c("red","yellow","green"))
ord.Fst <- rk.Fst.adapt[rank(as.numeric(nnet.Fst.adapt$tip.label))]
tipcols.Fst <- colfunc(dim(as.matrix(d.Fst.adapt))[1])[ord.Fst]
plot(nnet.Fst.adapt, tip.color = tipcols.Fst, type = '2D')


# Compare to neutral subsets
comp.dist.fst <- vector()
comp.rank.fst <- vector()

for(i in 1:1000){ 
  cat("\r", paste("Neutral SNP iteration number:", i, "   "))
  
  ## Pairwise Fst of SUBSET SNPs ##
  n.snps <- length(signif.snps)
  set.seed(i + 3434)
  sampled.snps <- sample(read.gdsn(index.gdsn(YWAR.gds, "snp.rs.id")), n.snps, replace = FALSE) #SAMPLING OF SNPS HERE
 
  d.Fst.sub <- as.dist(calculate.all.pairwise.Fst(geno.pop[,snps.gds %in% sampled.snps], 
                                                    chrom.pop[,snps.gds %in% sampled.snps]))
  rk.Fst.sub <- shap2ranks(shapley_index(neighborNet(d.Fst.sub)))
  
  # compare dist matrices
  mtest.i <- mantel(d.Fst.all, d.Fst.sub)
  comp.dist.fst[i] <- mtest.i$statistic
  
  # compare ranks
  comp.rank.fst[i] <- cor(rk.Fst.all, rk.Fst.sub, method = "kendall")
  
}

comp.subs.fst <- data.frame(iteration = c(1:length(comp.rank.fst)),
                            compare.dist = comp.dist.fst,
                            compare.ranks = comp.rank.fst)

#write.csv(comp.subs.fst, "compsubs_fst.csv", row.names = FALSE)
#comp.subs.fst <- read.csv("compsubs_fst.csv")

(pval <- length(which(comp.subs.fst$compare.dist < obs.dist.mantel.fst))/length(comp.subs.fst$compare.dist))

g.d.fst <- ggplot(comp.subs.fst, aes(compare.dist)) +
  geom_histogram(binwidth = 0.01, color = "grey15", fill = "grey15") +
  geom_vline(xintercept = obs.dist.mantel.fst, color = "red", size = 2) +
  xlab("Mantel correlation to distance matrix using all SNPs")+
  geom_text(x = min(comp.subs.fst$compare.dist)+((range(comp.subs.fst$compare.dist)[2]-range(comp.subs.fst$compare.dist)[1])/20),
            y = 10, label=paste("p =", pval), size = 6, fontface = "italic", color = "grey15")+
  theme_classic()
g.d.fst

(pval2 <- length(which(comp.subs.fst$compare.ranks < obs.rank.kendall.fst))/length(comp.subs.fst$compare.ranks))

g.r.fst <- ggplot(comp.subs.fst, aes(compare.ranks)) +
  geom_histogram(binwidth = 0.02, color = "grey15", fill = "grey15") +
  geom_vline(xintercept = obs.rank.kendall.fst, color = "red", size = 2) +
  xlab("Kendall Tau correlation to ranks using all SNPs")+
  geom_text(x = min(comp.subs.fst$compare.ranks)+((range(comp.subs.fst$compare.ranks)[2]-range(comp.subs.fst$compare.ranks)[1])/20),
            y = 10, label=paste("p =", pval2), size = 6, fontface = "italic", color = "grey15")+
  theme_classic() 
g.r.fst


#### Pairwise Jost (VCFR) ####
library(vcfR)

YWAR.vcf <- read.vcfR("YWAR_final.vcf")

miss.snp.g  <- apply(YWAR.vcf@gt, 1, function(x) length(grep("\\.", x))/length(x))
miss.samp.g <- apply(YWAR.vcf@gt, 2, function(x) length(grep("\\.", x))/length(x))
YWAR.vcf2 <- YWAR.vcf[which(miss.snp.g<0.1), which(miss.samp.g<0.1)]

pops.vcf <- vector()
for(i in 2:length(colnames(YWAR.vcf@gt))){
  pops.vcf <- append(pops.vcf,YWAR.meta$Pop[YWAR.meta$Field_Number==colnames(YWAR.vcf@gt)[i]])
}
pops.vcf <- factor(pops.vcf)

d.jost.all <- pairwise_JOST(YWAR.vcf2, pops.vcf) #long
d.jost.adapt <- pairwise_JOST(YWAR.vcf2[which(YWAR.vcf2@fix[,3] %in% signif.snps),], pops.vcf)

(obs.dist.mantel.jost <- mantel(d.jost.all, d.jost.adapt)$statistic)
  
nnet.jost.all <- neighborNet(d.jost.all)
rk.jost.all <- shap2ranks(shapley_index(nnet.jost.all))

nnet.jost.adapt <- neighborNet(d.jost.adapt)
rk.jost.adapt <- shap2ranks(shapley_index(nnet.jost.adapt))

(obs.rank.kendall.jost <- cor(rk.jost.all, rk.jost.adapt, method = "kendall"))

plot(nnet.jost.all)
plot(nnet.jost.adapt)

# Compare to neutral subsets

comp.dist.jost <- vector()
comp.rank.jost <- vector()

for(i in 1:1000){ 
  cat("\r", paste0(i/1000*100, "%  "))
  
  ## Pairwise Fst of SUBSET SNPs ##
  n.snps <- length(signif.snps)
  #set.seed(3434)
  sampled.snps <- sample(row.names(cands.all), n.snps, replace = FALSE) #SAMPLING OF SNPS HERE
  d.jost.sub <- pairwise_JOST(YWAR.vcf2[which(YWAR.vcf2@fix[,3] %in% sampled.snps),], pops.vcf)
  rk.jost.sub <- shap2ranks(shapley_index(neighborNet(d.jost.sub)))
  
  # compare dist matrices
  mtest.i <- mantel(d.jost.all, d.jost.sub)
  comp.dist.jost[i] <- mtest.i$statistic
  
  # compare ranks
  comp.rank.jost[i] <- cor(rk.jost.all, rk.jost.sub, method = "kendall")
  
}

comp.subs.jost <- data.frame(iteration = c(1:length(comp.rank.jost)),
                            compare.dist = comp.dist.jost,
                            compare.ranks = comp.rank.jost)

write.csv(comp.subs.jost, "compsubs_jost_vcfR.csv")

# Check significance
(pval <- length(which(comp.subs.jost$compare.dist < obs.dist.mantel.jost))/length(comp.subs.jost$compare.dist))

g.d.jost <- ggplot(comp.subs.jost, aes(compare.dist)) +
  geom_histogram(binwidth = 0.01, color = "grey15", fill = "grey15") +
  geom_vline(xintercept = obs.dist.mantel.jost, color = "red", size = 2) +
  xlab("Mantel correlation to distance matrix using all SNPs")+
  geom_text(x = min(comp.subs.jost$compare.dist)+((range(comp.subs.jost$compare.dist)[2]-range(comp.subs.jost$compare.dist)[1])/20),
            y = 10, label=paste("p =", pval), size = 6, fontface = "italic", color = "grey15")+
  theme_classic()
g.d.jost

(pval2 <- length(which(comp.subs.jost$compare.ranks < obs.rank.kendall.jost))/length(comp.subs.jost$compare.ranks))

g.r.jost <- ggplot(comp.subs.jost, aes(compare.ranks)) +
  geom_histogram(binwidth = 0.01, color = "grey15", fill = "grey15") +
  geom_vline(xintercept = obs.rank.kendall.jost, color = "red", size = 2) +
  xlab("Kendall Tau correlation to ranks using all SNPs")+
  geom_text(x = min(comp.subs.jost$compare.ranks)+((range(comp.subs.jost$compare.ranks)[2]-range(comp.subs.jost$compare.ranks)[1])/20),
            y = 10, label=paste("p =", pval2), size = 6, fontface = "italic", color = "grey15")+
  theme_classic() 
g.r.jost


#### END ####

snpgdsClose(YWAR.gds)


