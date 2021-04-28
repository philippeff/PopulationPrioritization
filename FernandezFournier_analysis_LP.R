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
library(RSplitsTree)
library(phangorn)
library(vegan)
library(svMisc)
library(maps)
library(raster)
library(hierfstat)
library(phangorn)

setwd("E:/sfuvault/SFU/MooersLab/PopGen/Analysis/Final R code - Cons Bio")

# Load functions
source("FernandezFournier_functions.R")

#PHENOTYPE SNP data
cands.pheno <- read.csv("pvalues_phenoSNPs.csv")

#Mahony 2019 populations metadata
pop.coord <- read.csv("pop_coord_LP.csv")

### GDS of Mahony 2019 genomicdata must be loaded every time ###
pine.gds <- snpgdsOpen("pine.gds", allow.duplicate = TRUE) 

#filter populations (selected at random to represent the largest area of distribution possible)
pop.select <- c(61, 171, 162, 26, 143, 245, 33, 14, 7, 208, 233, 277, 59, 55, 118, 248, 242, 70, 216, 279, 82, 30)
pops.gen <- read.gdsn(index.gdsn(pine.gds, "sample.annot"))$pop.group


#### Explore ####

miss.rate.samp <- apply(read.gdsn(index.gdsn(pine.gds, "genotype")), 1, 
                        function(x) length(which(x==3))/length(x)) #Proportion of missing SNPs per sample
miss.rate.snp <- apply(read.gdsn(index.gdsn(pine.gds, "genotype")), 2, 
                       function(x) length(which(x==3))/length(x)) # Proportion of missing alleles per SNP

hist(miss.rate.samp, breaks = 50, main = "Proportion of missing SNPs per sample")
hist(miss.rate.snp, breaks = 30,  main = "Proportion of missing alleles per SNP")

length(which(miss.rate.samp < 0.1))
length(which(miss.rate.snp  < 0.1))


#### PCA OF ALL SNPS ####
#Select for SNPs with <10% missing rate and Samples with <10% missing SNPs
pca.all <- snpgdsPCA(pine.gds, autosome.only = FALSE, bayesian = TRUE, missing.rate = 0.1, 
                     sample.id = read.gdsn(index.gdsn(pine.gds, "sample.id"))[pops.gen %in% pop.select & miss.rate.samp<0.1]) # retains all axes

# Vector of populations for each individual
pops.all <- unlist(lapply(strsplit(pca.all$sample.id, "[A-Z]"), function(x) x[1]))

# Choose significant PCs (Tracy-Widom distribution test)
# significance level of 0.05 (criticalpoint = 0.9793). Default critical point is set to significance level of 0.01 (criticalpoint = 2.0234)
tw.all <- tw(na.exclude(pca.all$eigenval), length(na.exclude(pca.all$eigenval)), criticalpoint = 0.9793)
(n.pc.all <- tw.all$SigntEigenL)

# Contribution to variance for each PC
pc.percent.all <- as.vector(na.exclude(pca.all$varprop*100))
barplot(pc.percent.all, 
        col = c(rep("grey50", n.pc.all), rep("grey90", length(pc.percent.all)-n.pc.all)),
        xlab = "Principal Components (dark grey = significant)",
        ylab = "Percent contribution")

# PCA plot
pca.all.gg <- as.data.frame(pca.all$eigenvect)
pca.all.gg$group <- as.character(pops.all)

p1<-ggplot(pca.all.gg,aes(x=V1,y=V2,color=group)) +
  geom_point(size=2) +
  scale_color_manual(values=funky(length(unique(pops.all)))) +
  xlab(paste("PC1 (", round(pc.percent.all[1], 2), "%)", sep = "")) +
  ylab(paste("PC2 (", round(pc.percent.all[2], 2), "%)", sep = "")) +
  xlim(c(-0.35, 0.2)) +
  ylim(c(-0.35, 0.15)) +
  theme_bw() +
  theme(legend.position="none") 
p1

#make dataset of PCA centroids and pairwise distance object
pca.all.cen <- aggregate(pca.all$eigenvect, list(pop = pops.all), mean)
pca.all.cen <- pca.all.cen[order(as.numeric(pca.all.cen$pop)),]
rownames(pca.all.cen) <- pca.all.cen$pop

d.pca.all <- stats::dist(pca.all.cen[,2:(n.pc.all+1)], method = "euclidean") 


#### NEIGHBORNET AND RANKS OF ALL SNPS ####
# write dist as Phylip distance file for SplitsTree
writeDist(d.pca.all, "d_pca_all_subpop.dist") 

# Open dist in SplitsTree and save as NEX file
splitstree_pff(dist.file = "d_pca_all_subpop.dist", nexus.out = "nnet_all_subpop.nex") #Function to compute neighborNet from dist file in SplitsTree

# NeighborNet
nnet.pca.all <- read.nexus.networx("nnet_all_subpop.nex", splits = TRUE)

# Ranks
rk.pca.all <- shap2ranks_pine(shapley_index(nnet.pca.all))

#Plot NNET and color tips by rank
colfunc<-colorRampPalette(c("red","yellow","green"))
ord.all <- rk.pca.all[rank(as.numeric(nnet.pca.all$tip.label))]
tipcols.all <- colfunc(dim(as.matrix(d.pca.all))[1])[ord.all]

plot(nnet.pca.all,tip.color=tipcols.all, type = '2D')



#### PCA OF PHENOTYPE SNPS ####
pheno.snps <- gsub("_", "-", cands.pheno$SNP[which(cands.pheno$is.signif==1)])

pca.pheno <- snpgdsPCA(pine.gds, autosome.only=FALSE, bayesian = TRUE, verbose = FALSE, missing.rate = 0.1, 
                       sample.id = read.gdsn(index.gdsn(pine.gds, "sample.id"))[pops.gen %in% pop.select & miss.rate.samp<0.1],
                       snp.id = pheno.snps) # retains all axes

# Choose significant PCs (Tracy-Widom distribution test)
# significance level of 0.05 (criticalpoint = 0.9793). Default critical point is set to significance level of 0.01 (criticalpoint = 2.0234)
tw.pheno <- tw(na.exclude(pca.pheno$eigenval), length(na.exclude(pca.pheno$eigenval)), criticalpoint = 0.9793)
(n.pc.pheno <- tw.pheno$SigntEigenL)

# Vector of populations for each individual
pops.pheno <- unlist(lapply(strsplit(pca.pheno$sample.id, "[A-Z]"), function(x) x[1]))

# Contribution of PCs
pc.percent.pheno <- as.vector(na.exclude(pca.pheno$varprop*100))
barplot(pc.percent.pheno, 
        col = c(rep("grey50", n.pc.pheno), rep("grey90", length(pc.percent.pheno)-n.pc.pheno)),
        xlab = "Principal Components (dark grey = significant)",
        ylab = "Percent contribution")

## PCA plot
pca.pheno.gg <- as.data.frame(pca.pheno$eigenvect)
pca.pheno.gg$group <- as.character(pops.pheno)

pheno1 <-ggplot(pca.pheno.gg,aes(x=-V1,y=V2,color=group)) +
  geom_point(size=2) +
  scale_color_manual(values=funky(length(unique(pops.pheno)))) +
  xlab(paste("PC1 (", round(pc.percent.pheno[1], 2), "%)", sep = "")) +
  ylab(paste("PC2 (", round(pc.percent.pheno[2], 2), "%)", sep = "")) +
  xlim(c(-0.35, 0.2)) +
  ylim(c(-0.35, 0.15)) +
  theme_bw() +
  theme(legend.position="none") 
pheno1


#make data set of centroids
pca.pheno.cen <- aggregate(pca.pheno$eigenvect, list(pop = pops.pheno), mean)
pca.pheno.cen <- pca.pheno.cen[order(as.numeric(pca.pheno.cen$pop)),]
rownames(pca.pheno.cen) <- pca.pheno.cen$pop
d.pca.pheno <- stats::dist(pca.pheno.cen[,2:(n.pc.pheno+1)], method = "euclidean") #How many axes do we retain?????


#### NEIGHBORNET AND RANKS OF PHENOTYPE SNPS ####

# write dist as Phylip distance file for SplitsTree
writeDist(d.pca.pheno, "d_pca_pheno_subpop.dist") 

# Open dist in SplitsTree and save as NEX file
splitstree_pff(dist.file = "d_pca_pheno_subpop.dist", nexus.out = "nnet_pheno_subpop.nex") #Function to compute neighborNet from dist file in SplitsTree (I'm a genius)

nnet.pca.pheno <- read.nexus.networx("nnet_pheno_subpop.nex", splits = TRUE)

rk.pca.pheno <- shap2ranks_pine(shapley_index(nnet.pca.pheno))

#Plot NNET and color tips by rank
colfunc <- colorRampPalette(c("red","yellow","green"))
ord.pheno <- rk.pca.pheno[rank(as.numeric(nnet.pca.pheno$tip.label))]
tipcols.pheno <- colfunc(dim(as.matrix(d.pca.pheno))[1])[ord.pheno]

plot(nnet.pca.pheno,tip.color=tipcols.pheno, type = '2D')


#### BOOTSTRAP TEST FOR PHENOTYPE SNPS ####
# Are the pairwise population distance matrix and population ranks MORE DIFFERENT from all SNPs than random subsets?

# Distance matrices
(obs.dist.mantel.pheno <- mantel(d.pca.all, d.pca.pheno)$statistic)

# Ranks
rk.pca.all
rk.pca.pheno
(obs.rank.kendall.pheno <- cor(rk.pca.all, rk.pca.pheno, method = "kendall"))

(ranks.table <- data.frame(pop = names(rk.pca.all),
                          ranks.all = rk.pca.all,
                          ranks.pheno = rk.pca.pheno))

# Compare to neutral subsets
# Draw many different sets of (any) SNPs of same size as adaptive SNPs

(n.phenosnps <- length(cands.pheno$SNP[which(cands.pheno$is.signif==1)]))

comp.subs.pca.pheno <- data.frame(iteration = c(),
                            compare.dist = c(),
                            compare.ranks = c(),
                            n.PC = c())

for(i in 1:10){ 
  print(paste("Iteration number", i))
  set.seed(i+12)
  
  ## PCA of SUBSET SNPs ##
  sampled.snps <- sample(gsub("_", "-", cands.pheno$SNP), n.phenosnps, replace = FALSE) #SAMPLING OF SNPS HERE
  
  pca.sub <- snpgdsPCA(pine.gds, autosome.only=FALSE, bayesian =  TRUE, verbose = FALSE, missing.rate = 0.1, 
                       sample.id = read.gdsn(index.gdsn(pine.gds, "sample.id"))[pops.gen %in% pop.select & miss.rate.samp<0.1],
                       snp.id = sampled.snps) # retains all axes
  
  # Vector of populations for each individual
  pops.sub <- unlist(lapply(strsplit(pca.sub$sample.id, "[A-Z]"), function(x) x[1]))
  
  ## CHOOSE SIGNIFICANT PRINCIPAL COMPONENTS
  #Tracy-Widom distribution test
  # significance level of 0.05 (criticalpoint = 0.9793). Default critical point is set to significance level of 0.01 (criticalpoint = 2.0234)
  tw.sub <- tw(na.exclude(pca.sub$eigenval), length(na.exclude(pca.sub$eigenval)), criticalpoint = 0.9793)
  (n.pc.sub <- tw.sub$SigntEigenL)
  
  #make data set of centroids
  pca.sub.cen <- aggregate(pca.sub$eigenvect, list(pop = pops.sub), mean)
  pca.sub.cen <- pca.sub.cen[order(as.numeric(pca.sub.cen$pop)),]
  rownames(pca.sub.cen) <- pca.sub.cen$pop
  d.pca.sub <- stats::dist(pca.sub.cen[,2:(n.pc.sub+1)], method = "euclidean") 

  #write dist as Phylip distance file for SplitsTree
  writeDist(d.pca.sub, "temp_d_pca_sub.dist") 
  
  #Open dist in SplitsTree and save as NEX file
  splitstree_pff(dist.file = "temp_d_pca_sub.dist", nexus.out = "temp_nnet_sub.nex") #Function to compute neighborNet from dist file in SplitsTree (Im a genius)
  
  nnet.pca.sub <- read.nexus.networx("temp_nnet_sub.nex", splits = TRUE)
  
  rk.pca.sub <- shap2ranks_pine(shapley_index(nnet.pca.sub))
  
  # compare dist matrices
  mtest.i <- mantel(d.pca.all, d.pca.sub)
  comp.dist <- mtest.i$statistic
  
  # compare ranks
  comp.rank <- cor(rk.pca.all, rk.pca.sub, method = "kendall")
  
  comp.subs.pca.pheno <- rbind(comp.subs.pca.pheno, data.frame(iteration = i,
                                                   compare.dist = comp.dist,
                                                   compare.ranks = comp.rank,
                                                   n.PC = n.pc.sub))
  
  write.csv(comp.subs.pca.pheno, "compsubs_pca_mahony_subpop_pheno.csv", row.names = FALSE)

  invisible(file.remove("temp_d_pca_sub.dist"))
  invisible(file.remove("temp_nnet_sub.nex"))
}


# Results and plots
(pheno.pval <- length(which(comp.subs.pca.pheno$compare.dist < obs.dist.mantel.pheno))/length(comp.subs.pca.pheno$compare.dist))

pheno.g1 <- ggplot(comp.subs.pca.pheno, aes(compare.dist)) +
  geom_histogram(binwidth = 0.01, color = "grey25", fill = "grey25") +
  geom_vline(xintercept = obs.dist.mantel.pheno, color = "red", size = 2) +
  xlab("Phenotype SNPs vs. all SNPs (mantel correlation of dist matrices)")+
  geom_text(x = 0.79, y = 25, label=paste("p =", pheno.pval), size = 6, fontface = "italic")+
  theme_classic()
pheno.g1

(pheno.pval2 <- length(which(comp.subs.pca.pheno$compare.ranks < obs.rank.kendall.pheno))/length(comp.subs.pca.pheno$compare.ranks))

pheno.g2 <- ggplot(comp.subs.pca.pheno, aes(compare.ranks)) +
  geom_histogram(binwidth = 0.02, color = "grey25", fill = "grey25") +
  geom_vline(xintercept = obs.rank.kendall.pheno, color = "red", size = 2) +
  xlab("Phenotype SNPs vs. all SNPs (Kendall tau of ranks)")+
  geom_text(x = min(comp.subs.pca.pheno$compare.ranks)+((range(comp.subs.pca.pheno$compare.ranks)[2]-range(comp.subs.pca.pheno$compare.ranks)[1])/20),
            y = 20, label=paste("p =", pheno.pval2), size = 6, fontface = "italic")+
  theme_classic() 
pheno.g2


#### MAPS ####
library(sf)
library(rnaturalearth)
library(png)
library(berryFunctions)

pop.coordsub <- pop.coord[pop.coord$id1 %in% pop.select,]

#Jitter
pop.coordsub$lat[pop.coordsub$id1 == "7"] <- pop.coordsub$lat[pop.coordsub$id1 == "7"] + 0.2
pop.coordsub$lon[pop.coordsub$id1 == "7"] <- pop.coordsub$lon[pop.coordsub$id1 == "7"] + 0.2

topo <- ne_download(scale = 10, type = 'GRAY_HR_SR_W', category = 'raster')

colfunc<-colorRampPalette(c("red","yellow","green"))
colfunc2 <- grey.colors(5, start = 1, end = 0.4)
colfunc2[1:2] <- "#FFFFFF"

pine.png <- readPNG("E:/sfuvault/SFU/MooersLab/PopGen/Graphs/pine.png")
op <- par(no.readonly=TRUE) # original parameters

#Map of Ranks using ALL SNPs
raster::plot(topo, xlim=c(-133.8,-109), ylim=c(48,60),
             col = colfunc2, xlab = "Longitude", ylab = "Latitude", legend = FALSE)
maps::map(xlim=c(-133.8,-109), ylim=c(48,60), fill = FALSE, lwd = 1, add = TRUE) 
points(x=pop.coordsub$lon, y = pop.coordsub$lat, pch = 20, cex = 3, 
       col = colfunc(length(pop.coordsub$id1))[rk.pca.all])
text(x = pop.coordsub$lon, y = pop.coordsub$lat, labels = pop.coordsub$id1, cex = 0.5)
text(x = -129.3, y = 48.3, labels = "All genetic variation", cex = 0.7)
rasterImage(pine.png, -131.6, 48.6, -128.6, 52)
# Inset
par(usr=c(-400, -48, -150, 75)) # c(x1, x2, y1, y2) giving the extremes of the user coordinates of the plotting region.
rect(xleft =-170, xright = -48, ybottom = 11, ytop = 75, col = "white")
maps::map(xlim=c(-150,-48), ylim=c(15,75), fill = TRUE, col = "grey95", add = TRUE)
rect(xleft =-135, xright = -115, ybottom = 48, ytop = 60)
par(op)


#Map of Ranks using PHENOTYPE SNPs
raster::plot(topo, xlim=c(-133.8,-109), ylim=c(48,60),
             col = colfunc2, xlab = "Longitude", ylab = "Latitude", legend = FALSE)
maps::map(xlim=c(-133.8,-109), ylim=c(48,60), fill = FALSE, lwd = 1, add = TRUE) 
points(x=pop.coordsub$lon, y = pop.coordsub$lat, pch = 20, cex = 3, 
       col = colfunc(length(pop.coordsub$id1))[rk.pca.pheno])
text(x = pop.coordsub$lon, y = pop.coordsub$lat, labels = pop.coordsub$id1, cex = 0.5)
text(x = -129.3, y = 48.3, labels = "Adaptive genetic variation", cex = 0.7)
rasterImage(pine.png, -131.6, 48.6, -128.6, 52)

# Inset
par(usr=c(-400, -48, -150, 75)) # c(x1, x2, y1, y2) giving the extremes of the user coordinates of the plotting region.
rect(xleft =-170, xright = -48, ybottom = 11, ytop = 75, col = "white")
maps::map(xlim=c(-150,-48), ylim=c(15,75), fill = TRUE, col = "grey95", add = TRUE)
rect(xleft =-135, xright = -115, ybottom = 48, ytop = 60)
par(op)



#### APPENDIX ####

#### OVERALL Fst ####

fst.all <- snpgdsFst(pine.gds, as.factor(pops.all), method = "W&C84", autosome.only = FALSE, missing.rate = 0.1, 
                     sample.id = read.gdsn(index.gdsn(pine.gds, "sample.id"))[miss.rate.samp<0.1 & pops.gen %in% pop.select]) 
fst.all$Fst
fst.all$MeanFst

#### FST PHENOTYPE SNPS VS. ALL SNPS ####
pine.gds$root
pheno.snps <- gsub("_", "-", cands.pheno$SNP[which(cands.pheno$is.signif==1)])
snps.gds <- read.gdsn(index.gdsn(pine.gds, "snp.id"))
snps.select <- snps.gds[snps.gds %in% pheno.snps]

geno <- read.gdsn(index.gdsn(pine.gds, "genotype"))
colnames(geno) <- read.gdsn(index.gdsn(pine.gds, "snp.id"))
rownames(geno) <- read.gdsn(index.gdsn(pine.gds, "sample.id"))
geno[geno == 3] <- NA # 3 means missing data

pop.fst <- unlist(lapply(strsplit(rownames(geno), "[A-Z]"), function(x) x[1]))
geno <- geno[pop.fst %in% pop.select,]
pop.fst <- unlist(lapply(strsplit(rownames(geno), "[A-Z]"), function(x) x[1]))
pop.fst <- factor(pop.fst)

# allele count
geno.pop <- matrix(0L, nrow = length(levels(pop.fst)), ncol = dim(geno)[2])
for(i in 1:length(levels(pop.fst))){
  geno.pop[i,] <- apply(geno[which(pop.fst == levels(pop.fst)[i]),], 2, function(x) sum(x, na.rm=TRUE))
}
colnames(geno.pop) <- colnames(geno)
rownames(geno.pop) <- as.character(levels(pop.fst))

# chromosome sample size
chrom.pop <- matrix(0L, nrow = length(levels(pop.fst)), ncol = dim(geno)[2])
for(i in 1:length(levels(pop.fst))){
  chrom.pop[i,] <- apply(geno[which(pop.fst == levels(pop.fst)[i]),], 2, function(x) length(which(!is.na(x)))*2)
}
dimnames(chrom.pop) <- dimnames(geno.pop)

# distances
d.Fst.all <- as.dist(calculate.all.pairwise.Fst(geno.pop, chrom.pop))

d.Fst.pheno <- as.dist(calculate.all.pairwise.Fst(geno.pop[,snps.gds %in% pheno.snps], 
                                                  chrom.pop[,snps.gds %in% pheno.snps]))


writeDist(d.Fst.all, "d_fst_all.dist") 
writeDist(d.Fst.pheno, "d_fst_pheno.dist") 

d.Fst.all <- readDist("d_fst_all.dist") 
d.Fst.pheno <- readDist("d_fst_pheno.dist") 

splitstree_pff(dist.file = "d_fst_all.dist", nexus.out = "nnet_fst_all.nex") 
splitstree_pff(dist.file = "d_fst_pheno.dist", nexus.out = "nnet_fst_pheno.nex") 

nnet.fst.all   <- read.nexus.networx("nnet_fst_all.nex", splits = TRUE)
nnet.fst.pheno <- read.nexus.networx("nnet_fst_pheno.nex", splits = TRUE)

plot(nnet.fst.all, type = "2D")
plot(nnet.fst.pheno, type = "2D")

rk.fst.all   <- shap2ranks_pine(shapley_index(nnet.fst.all))
rk.fst.pheno <- shap2ranks_pine(shapley_index(nnet.fst.pheno))

(obs.dist.mantel.fst.pheno <- mantel(d.Fst.all, d.Fst.pheno)$statistic)

(obs.rank.kendall.fst.pheno <- cor(rk.fst.all, rk.fst.pheno, method = "kendall"))

# Compare to neutral subsets
n.snps <- length(cands.pheno$SNP[which(cands.pheno$is.signif==1)])

comp.subs.fst.pheno <- data.frame(iteration = c(),
                                  compare.dist = c(),
                                  compare.ranks = c())

for(i in 1:1000){ 
  print(paste("Iteration number", i))
  set.seed(i+42)
  
  sampled.snps <- sample(gsub("_", "-", cands.pheno$SNP), n.snps, replace = FALSE) #SAMPLING OF SNPS HERE
  d.Fst.sub <- as.dist(calculate.all.pairwise.Fst(geno.pop[,snps.gds %in% sampled.snps], 
                                                    chrom.pop[,snps.gds %in% sampled.snps]))
  writeDist(d.Fst.sub, "d_fst_sub.dist") 
  splitstree_pff(dist.file = "d_fst_sub.dist", nexus.out = "nnet_fst_sub.nex") 
  nnet.fst.sub <- read.nexus.networx("nnet_fst_sub.nex", splits = TRUE)
  rk.fst.sub <- shap2ranks_pine(shapley_index(nnet.fst.sub))
  
  # compare dist matrices
  mtest.i <- mantel(d.Fst.all, d.Fst.sub)
  comp.dist <- mtest.i$statistic
  
  # compare ranks
  comp.rank <- cor(rk.fst.all, rk.fst.sub, method = "kendall")
  comp.subs.fst.pheno <- rbind(comp.subs.fst.pheno, data.frame(iteration = i,
                                                               compare.dist = comp.dist,
                                                               compare.ranks = comp.rank))
  write.csv(comp.subs.fst.pheno, "compsubs_fst_mahony_subpop_pheno.csv", row.names = FALSE)
  
  invisible(file.remove("d_fst_sub.dist"))
  invisible(file.remove("nnet_fst_sub.nex"))
}


# Open the comparisons dataset
comp.subs.fst.pheno <- read.csv("compsubs_fst_mahony_subpop_pheno.csv")

(pheno.fst.pval <- length(which(comp.subs.fst.pheno$compare.dist < obs.dist.mantel.fst.pheno))/length(comp.subs.fst.pheno$compare.dist))

g.d.fst.pheno <- ggplot(comp.subs.fst.pheno, aes(compare.dist)) +
  geom_histogram(binwidth = 0.01, color = "grey15", fill = "grey15") +
  geom_vline(xintercept = obs.dist.mantel.fst.pheno, color = "red", size = 2) +
  xlab("Mantel correlation to distance matrix using all SNPs")+
  geom_text(x = min(comp.subs.fst.pheno$compare.dist)+((range(comp.subs.fst.pheno$compare.dist)[2]-range(comp.subs.fst.pheno$compare.dist)[1])/20),
            y = 10, label=paste("p =", pheno.fst.pval), size = 6, fontface = "italic", color = "grey15")+
  theme_classic()
g.d.fst.pheno

(pheno.fst.pval2 <- length(which(comp.subs.fst.pheno$compare.ranks < obs.rank.kendall.fst.pheno))/length(comp.subs.fst.pheno$compare.ranks))

g.r.fst.pheno <- ggplot(comp.subs.fst.pheno, aes(compare.ranks)) +
  geom_histogram(binwidth = 0.02, color = "grey15", fill = "grey15") +
  geom_vline(xintercept = obs.rank.kendall.fst.pheno, color = "red", size = 2) +
  xlab("Kendall Tau correlation to ranks using all SNPs")+
  geom_text(x = min(comp.subs.fst.pheno$compare.ranks)+((range(comp.subs.fst.pheno$compare.ranks)[2]-range(comp.subs.fst.pheno$compare.ranks)[1])/20),
            y = 10, label=paste("p =", pheno.fst.pval2), size = 6, fontface = "italic", color = "grey15")+
  theme_classic() 
g.r.fst.pheno




