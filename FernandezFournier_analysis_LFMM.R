############Install the BiocMananager package in order to install the LEA package############
install.packages("BiocManager", repos = "http://cran.rstudio.com/", lib = "/home")
library(BiocManager)

############Install LEA package############
BiocManager::install("LEA")
library("LEA")

######Convert the .vcf file into a .geno file, and then to a .lfmm file to use in LFMM######
vcf2geno("YWAR_final.vcf","YWAR_final.geno")
geno2lfmm("YWAR_final.geno", "YWAR_final.lfmm")

######Download necessary packages to process WorldClim data######
install.packages("sp",repos="http://cran.stat.sfu.ca/")
install.packages("raster",repos="http://cran.stat.sfu.ca/")
install.packages("rgdal",repos="http://cran.stat.sfu.ca/")
library(sp)
library(raster)
library(rgdal)

######Get WorldClim data for each of the 19 variables at a resolution of 10 (minutes of a degree)######
all_clim <- getData("worldclim",var="bio",res=10)
names(all_clim) <- c("Annual_Mean_Temp","Mean_Diurnal_Range","Isothermality", "Temp_Seasonality", "Max_Temp_Warmest_Month", "Min_Temp_Coldest_Month", "Temp_Annual_Range", "Mean_Temp_Wettest_Quarter", "Mean_Temp_Driest_Quarter", "Mean_Temp_Warmest_Quarter","Mean_Temp_Coldest_Quarter", "Annual_Prec", "Prec_Wettest_Month", "Prec_Driest_Month", "Prec_ Seasonality", "Prec_Wettest_Quarter", "Prec_Driest_Quarter", "Prec_Warmest_Quarter", "Prec_Coldest_Quarter")




######Upload location data and extract environmental variables at each of these locations######
points <- read.csv(file="YWAR_RAD_meta_filter.csv", header=TRUE, sep=",")
coordinates(points)<-~Long+Lat
values <- extract(all_clim,points)


######Save these extracted values as a .env file to run in LFMM######
write.env(values,"gradients.env")



############Commence LFMM############

project = lfmm("YWAR_final.lfmm",
               "gradients.env",
               K = 2,
               repetitions = 5,
               it=10000,
               project = "new")



######When it is completed, load the project######
project = load.lfmmProject("YWAR_final_gradients.lfmmProject")



######To get p-values, extract the z-scores and combine them using the median value######
zs.table_1 = z.scores(project,d=1, K=2)
zs.table_2 = z.scores(project,d=2, K=2)
zs.table_3 = z.scores(project,d=3, K=2)
zs.table_4 = z.scores(project,d=4, K=2)
zs.table_5 = z.scores(project,d=5, K=2)
zs.table_6 = z.scores(project,d=6, K=2)
zs.table_7 = z.scores(project,d=7, K=2)
zs.table_8 = z.scores(project,d=8, K=2)
zs.table_9 = z.scores(project,d=9, K=2)
zs.table_10 = z.scores(project,d=10, K=2)
zs.table_11 = z.scores(project,d=11, K=2)
zs.table_12 = z.scores(project,d=12, K=2)
zs.table_13 = z.scores(project,d=13, K=2)
zs.table_14 = z.scores(project2,d=14, K=2)
zs.table_15 = z.scores(project2,d=15, K=2)
zs.table_16 = z.scores(project2,d=16, K=2)
zs.table_17 = z.scores(project3,d=17, K=2)
zs.table_18 = z.scores(project3,d=18, K=2)
zs.table_19 = z.scores(project3,d=19, K=2)

zs.median_1 = apply(zs.table_1, MARGIN = 1, median)
zs.median_2 = apply(zs.table_2, MARGIN = 1, median)
zs.median_3 = apply(zs.table_3, MARGIN = 1, median)
zs.median_4 = apply(zs.table_4, MARGIN = 1, median)
zs.median_5 = apply(zs.table_5, MARGIN = 1, median)
zs.median_6 = apply(zs.table_6, MARGIN = 1, median)
zs.median_7 = apply(zs.table_7, MARGIN = 1, median)
zs.median_8 = apply(zs.table_8, MARGIN = 1, median)
zs.median_9 = apply(zs.table_9, MARGIN = 1, median)
zs.median_10 = apply(zs.table_10, MARGIN = 1, median)
zs.median_11 = apply(zs.table_11, MARGIN = 1, median)
zs.median_12 = apply(zs.table_12, MARGIN = 1, median)
zs.median_13 = apply(zs.table_13, MARGIN = 1, median)
zs.median_14 = apply(zs.table_14, MARGIN = 1, median)
zs.median_15 = apply(zs.table_15, MARGIN = 1, median)
zs.median_16 = apply(zs.table_16, MARGIN = 1, median)
zs.median_17 = apply(zs.table_17, MARGIN = 1, median)
zs.median_18 = apply(zs.table_18, MARGIN = 1, median)
zs.median_19 = apply(zs.table_19, MARGIN = 1, median)

#Lambda is your Genomic Inflation Factor (GIF); need to correct your p-values for this
#GIF: artificial differences in allele frequencies due to population stratification, cryptic relatedness and genotyping errors will affect all SNPs and so the test statistics will be inflated across the whole genome
#0·456 corresponds to the median of the chi-square distribution. P-values are correctly calibrated when the inflation factor is close to one
lambda_1 = median(zs.median_1^2)/0.456
lambda_2 = median(zs.median_2^2)/0.456
lambda_3 = median(zs.median_3^2)/0.456
lambda_4 = median(zs.median_4^2)/0.456
lambda_5 = median(zs.median_5^2)/0.456
lambda_6 = median(zs.median_6^2)/0.456
lambda_7 = median(zs.median_7^2)/0.456
lambda_8 = median(zs.median_8^2)/0.456
lambda_9 = median(zs.median_9^2)/0.456
lambda_10 = median(zs.median_10^2)/0.456
lambda_11 = median(zs.median_11^2)/0.456
lambda_12 = median(zs.median_12^2)/0.456
lambda_13 = median(zs.median_13^2)/0.456
lambda_14 = median(zs.median_14^2)/0.456
lambda_15 = median(zs.median_15^2)/0.456
lambda_16 = median(zs.median_16^2)/0.456
lambda_17 = median(zs.median_17^2)/0.456
lambda_18 = median(zs.median_18^2)/0.456
lambda_19 = median(zs.median_19^2)/0.456

adj.p.values_1 = pchisq(zs.median_1^2/lambda_1, df = 1, lower = FALSE)
adj.p.values_2 = pchisq(zs.median_2^2/lambda_2, df = 1, lower = FALSE)
adj.p.values_3 = pchisq(zs.median_3^2/lambda_3, df = 1, lower = FALSE)
adj.p.values_4 = pchisq(zs.median_4^2/lambda_4, df = 1, lower = FALSE)
adj.p.values_5 = pchisq(zs.median_5^2/lambda_5, df = 1, lower = FALSE)
adj.p.values_6 = pchisq(zs.median_6^2/lambda_6, df = 1, lower = FALSE)
adj.p.values_7 = pchisq(zs.median_7^2/lambda_7, df = 1, lower = FALSE)
adj.p.values_8 = pchisq(zs.median_8^2/lambda_8, df = 1, lower = FALSE)
adj.p.values_9 = pchisq(zs.median_9^2/lambda_9, df = 1, lower = FALSE)
adj.p.values_10 = pchisq(zs.median_10^2/lambda_10, df = 1, lower = FALSE)
adj.p.values_11 = pchisq(zs.median_11^2/lambda_11, df = 1, lower = FALSE)
adj.p.values_12 = pchisq(zs.median_12^2/lambda_12, df = 1, lower = FALSE)
adj.p.values_13 = pchisq(zs.median_13^2/lambda_13, df = 1, lower = FALSE)
adj.p.values_14 = pchisq(zs.median_14^2/lambda_14, df = 1, lower = FALSE)
adj.p.values_15 = pchisq(zs.median_15^2/lambda_15, df = 1, lower = FALSE)
adj.p.values_16 = pchisq(zs.median_16^2/lambda_16, df = 1, lower = FALSE)
adj.p.values_17 = pchisq(zs.median_17^2/lambda_17, df = 1, lower = FALSE)
adj.p.values_18 = pchisq(zs.median_18^2/lambda_18, df = 1, lower = FALSE)
adj.p.values_19 = pchisq(zs.median_19^2/lambda_19, df = 1, lower = FALSE)

write.csv(adj.p.values_1,"adj.p.values_1.csv")
write.csv(adj.p.values_2,"adj.p.values_2.csv")
write.csv(adj.p.values_3,"adj.p.values_3.csv")
write.csv(adj.p.values_4,"adj.p.values_4.csv")
write.csv(adj.p.values_5,"adj.p.values_5.csv")
write.csv(adj.p.values_6,"adj.p.values_6.csv")
write.csv(adj.p.values_7,"adj.p.values_7.csv")
write.csv(adj.p.values_8,"adj.p.values_8.csv")
write.csv(adj.p.values_9,"adj.p.values_9.csv")
write.csv(adj.p.values_10,"adj.p.values_10.csv")
write.csv(adj.p.values_11,"adj.p.values_11.csv")
write.csv(adj.p.values_12,"adj.p.values_12.csv")
write.csv(adj.p.values_13,"adj.p.values_13.csv")
write.csv(adj.p.values_14,"adj.p.values_14.csv")
write.csv(adj.p.values_15,"adj.p.values_15.csv")
write.csv(adj.p.values_16,"adj.p.values_16.csv")
write.csv(adj.p.values_17,"adj.p.values_17.csv")
write.csv(adj.p.values_18,"adj.p.values_18.csv")
write.csv(adj.p.values_19,"adj.p.values_19.csv")


#######Next, get the list of candidate loci (significantly associated with environment) from across all 104385 sites#######

######First, a Benjamini-Hochberg test to correct for false positives#######
#A q-value (the false discovery rate) of 0.1% (0.001) means that 0.1% of significant results will be false positives######

q = 0.001

L = length(adj.p.values_1)
w = which(sort(adj.p.values_1) < q * (1:L) / L)
candidates_1 = order(adj.p.values_1)[w]
write.csv(candidates_1,"candidates_1.csv")

L = length(adj.p.values_2)
w = which(sort(adj.p.values_2) < q * (1:L) / L)
candidates_2 = order(adj.p.values_2)[w]
write.csv(candidates_2,"candidates_2.csv")

L = length(adj.p.values_3)
w = which(sort(adj.p.values_3) < q * (1:L) / L)
candidates_3 = order(adj.p.values_3)[w]
write.csv(candidates_3,"candidates_3.csv")


L = length(adj.p.values_4)
w = which(sort(adj.p.values_4) < q * (1:L) / L)
candidates_4 = order(adj.p.values_4)[w]
write.csv(candidates_4,"candidates_4.csv")

L = length(adj.p.values_5)
w = which(sort(adj.p.values_5) < q * (1:L) / L)
candidates_5 = order(adj.p.values_5)[w]
write.csv(candidates_5,"candidates_5.csv")

L = length(adj.p.values_6)
w = which(sort(adj.p.values_6) < q * (1:L) / L)
candidates_6 = order(adj.p.values_6)[w]
write.csv(candidates_6,"candidates_6.csv")

L = length(adj.p.values_7)
w = which(sort(adj.p.values_7) < q * (1:L) / L)
candidates_7 = order(adj.p.values_7)[w]
write.csv(candidates_7,"candidates_7.csv")

L = length(adj.p.values_8)
w = which(sort(adj.p.values_8) < q * (1:L) / L)
candidates_8 = order(adj.p.values_8)[w]
write.csv(candidates_8,"candidates_8.csv")

L = length(adj.p.values_9)
w = which(sort(adj.p.values_9) < q * (1:L) / L)
candidates_9 = order(adj.p.values_9)[w]
write.csv(candidates_9,"candidates_9.csv")

L = length(adj.p.values_10)
w = which(sort(adj.p.values_10) < q * (1:L) / L)
candidates_10 = order(adj.p.values_10)[w]
write.csv(candidates_10,"candidates_10.csv")

L = length(adj.p.values_11)
w = which(sort(adj.p.values_11) < q * (1:L) / L)
candidates_11 = order(adj.p.values_11)[w]
write.csv(candidates_11,"candidates_11.csv")

L = length(adj.p.values_12)
w = which(sort(adj.p.values_12) < q * (1:L) / L)
candidates_12 = order(adj.p.values_12)[w]
write.csv(candidates_12,"candidates_12.csv")

L = length(adj.p.values_13)
w = which(sort(adj.p.values_13) < q * (1:L) / L)
candidates_13 = order(adj.p.values_13)[w]
write.csv(candidates_13,"candidates_13.csv")

L = length(adj.p.values_14)
w = which(sort(adj.p.values_14) < q * (1:L) / L)
candidates_14 = order(adj.p.values_14)[w]
write.csv(candidates_14,"candidates_14.csv")
L = length(p_14$pvalues)
w = which(sort(p_14$pvalues) < q * (1:L) / L)
candidates_14_v2 = order(p_14$pvalues)[w]

L = length(adj.p.values_15)
w = which(sort(adj.p.values_15) < q * (1:L) / L)
candidates_15 = order(adj.p.values_15)[w]
write.csv(candidates_15,"candidates_15.csv")

L = length(adj.p.values_16)
w = which(sort(adj.p.values_16) < q * (1:L) / L)
candidates_16 = order(adj.p.values_16)[w]
write.csv(candidates_16,"candidates_16.csv")

L = length(adj.p.values_17)
w = which(sort(adj.p.values_17) < q * (1:L) / L)
candidates_17 = order(adj.p.values_17)[w]
write.csv(candidates_17,"candidates_17.csv")

L = length(adj.p.values_18)
w = which(sort(adj.p.values_18) < q * (1:L) / L)
candidates_18 = order(adj.p.values_18)[w]
write.csv(candidates_18,"candidates_18.csv")

L = length(adj.p.values_19)
w = which(sort(adj.p.values_19) < q * (1:L) / L)
candidates_19 = order(adj.p.values_19)[w]
write.csv(candidates_19,"candidates_19.csv")
