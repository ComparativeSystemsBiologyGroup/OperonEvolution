'
5.	Average of log2 ratios (mean log2 ori/ter), distance data log transformed (pre)
Distance	Correlation	P-value	Rsq
1/ First Qu. of log dist	0.1693613	0.02265	0.02826573
1/Third Qu. of log dist	-0.1030724	0.1673	
1/Mean of log dist	0.1937886	0.008951	
1/Median of log dist	0.2135094	0.003903	0.04074252
'

## data frame with nGenomes, 1stQu, 3rdQu, mean, median, per species
df <- read.csv("distances_per_species_norm_log.csv", header=F, sep=",")
colnames(df) <- c("Species", "Genomes", "FirstQu", "ThirdQu", "Mean", "Median")
rownames(df) <- df$Species

## Subset by presence in ORI data frame
load("D:/Dottorato/Operons/ORITER_DATA.RData")
ori = data.frame(ORI)
m <- merge(ori, df, by = 0, all=T)
m <- na.omit(m)

## correlation Mean log2 ori/ter -- 1/First Qu.
cor.test(m$mean.log2.ori.ter, (1/m$FirstQu), method = "pearson")

## correlation Mean log2 ori/ter -- 1/Third Qu.
cor.test(m$mean.log2.ori.ter, (1/m$ThirdQu), method = "pearson")

## correlation Mean log2 ori/ter -- 1/Mean
cor.test(m$mean.log2.ori.ter, (1/m$Mean), method = "pearson")

## correlation Mean log2 ori/ter -- 1/Median
cor.test(m$mean.log2.ori.ter, (1/m$Median), method = "pearson")

###############################################################
'
8.	2^log(ori/ter), distance data log transformed (pre)
Distance	Correlation	P-value	Rsq
1/First Qu. of log dist	0.1681241	0.02368	0.02720598
1/Third Qu. of log dist	-0.1002381	0.1794	
1/Mean of log dist	0.190173	0.01034	
1/Median of log dist	0.2018478	0.006433	0.04885576
'


## data frame with nGenomes, 1stQu, 3rdQu, mean, median, per species
df <- read.csv("distances_per_species_norm_log.csv", header=F, sep=",")
colnames(df) <- c("Species", "Genomes", "FirstQu", "ThirdQu", "Mean", "Median")
rownames(df) <- df$Species

## Subset by presence in ORI data frame
load("D:/Dottorato/Operons/ORITER_DATA.RData")
ori = data.frame(ORI)
m <- merge(ori, df, by = 0, all=T)
m <- na.omit(m)

## correlation Mean 2^ log2 ori/ter -- 1/First Qu.
cor.test(2^log(m$mean.log2.ori.ter), (1/m$FirstQu), method = "pearson")

## correlation Mean 2^ log2 ori/ter -- 1/Third Qu.
cor.test(2^log(m$mean.log2.ori.ter), (1/m$ThirdQu), method = "pearson")

## correlation Mean 2^ log2 ori/ter -- 1/Mean
cor.test(2^log(m$mean.log2.ori.ter), (1/m$Mean), method = "pearson")

## correlation Mean 2^ log2 ori/ter -- 1/Median
cor.test(2^log(m$mean.log2.ori.ter), (1/m$Median), method = "pearson")

