source("https://bioconductor.org/biocLite.R")
biocLite("FlowSorted.CordBlood.450k")
biocLite("minfi")
biocLite("SummarizedExperiment")
biocLite("minfiData")
library(minfi)
library(minfiData)

baseDir<- "/Users/bix/Desktop/raw_0021/Allison.Appleton.0036/1520-1521 Appleton EPIC DNA methylation Data Package/ALL_207_IDATS/"
list.files(baseDir)
list.files(file.path(baseDir, "200498360121"))
targets <- read.metharray.sheet(baseDir)
sub(baseDir, "", targets$Basename)
RGSet <- read.metharray.exp(targets = targets, force = T)
pd <- pData(RGSet)
pd[,1:4] # phenodata

save(RGSet, file = "Allison.Appleton.0036.RGChannelSet_IDAT_Read_in.RData")


#### Get raw beta values:
Mset = preprocessRaw(RGSet)
RatioSet = ratioConvert(Mset, what = "both", keepCN = TRUE)
RawBetas = getBeta(RatioSet)
qc = getQC(Mset)
png(file="Sample_QC_Plot.png",width=1400,height=700,pointsize=12)
plotQC(qc)
dev.off()

save(pd, Mset, file = "Allison.Appleton0036.phenoData.mappedSet.RData")












## Functional Normalization:
Mset.norm = preprocessFunnorm(RGSet, sex = NULL, bgCorr = T, dyeCorr = T, nPCs = 2, verbose = TRUE)
FunNormBetas = getBeta(Mset.norm)

## Check distributions before and after normalization
library(RColorBrewer)	
png(file="BetaValue_Distributions_BeforeQC.png",width=1400,height=700,pointsize=12)
par(mfrow=c(1,2))
densityPlot(RawBetas, sampGroups = pd$Slide, main = "Raw Betas", xlab = "Beta")
densityPlot(FunNormBetas, sampGroups = pd$Slide, main = "FunNorm Adjusted Betas", xlab = "Beta")
dev.off()

## Generate Qc report:
qcReport(RGSet, sampNames = pd$Sample_Name, sampGroups = pd$Sample_Plate, pdf = "qcReport.pdf")









## Determine detection p-values within minfi
## Probes with a detection p-value abova signifigance (0.01) should not be trusted

detect.p = detectionP(RGSet, type = "m+u")
binning <- detect.p
binning[binning < 0.01] = 0
binning[binning > 0.01] = 1

prop.table = NULL
for (i in 1:length(colnames(detect.p))){
  table = table(binning[,i])
  prop.table.new = prop.table(table)
  prop.table = rbind(prop.table, prop.table.new)
}
rownames(prop.table) = colnames(binning)
prop.table = as.data.frame(prop.table)
prop.table[order(prop.table[,1]),]
colnames(prop.table) = c("Pass", "Fail")

## Review prop.table to determine which samples should be dropped
list <- subset(prop.table, Fail > 0.01) ## One sample at ~2.2% of probes fail...
fail = rownames(list)
RGSet = RGSet[,-c(which(rownames(pd) %in% fail))] ## Remove samples with poor detection, then re-run FunNorm:

## Estimate Cell Composition
cc = estimateCellCounts(RGSet, compositeCellType = "CordBlood")
#save(prop.table,detect.p,pd,FunNormBetas.sub,Mset.fn.sub, file = "Allison.Appleton.1_BasicQC_and_FunctionalNormalization.RData")




## Re-run FunNorm:
Mset.norm = preprocessFunnorm(RGSet, sex = NULL, bgCorr = T, dyeCorr = T, nPCs = 2, verbose = TRUE)
FunNormBetas = getBeta(Mset.norm)
Mset = preprocessRaw(RGSet)
RatioSet = ratioConvert(Mset, what = "both", keepCN = TRUE)
RawBetas = getBeta(RatioSet)

## Identify and drop probes with at least one probe with poor detection:
detect.p = detectionP(RGSet, type = "m+u")  ##### Starting with 866,836 Probes
## failed <- detect.p >0.01 ## how many failed probes
## colMeans(failed) ## Fraction of failed positions per sample
## sum(rowMeans(failed)>0.5)
## Create a matrix eliminating probes with a p-value greater than 0.01.
## You can then use this matrix after FunNorm to eliminate bad probes.
detect.p[detect.p > 0.01] = NA
filtered.p = na.omit(detect.p)

## Use the above matrix to exclude bad probes from the RawBeta and FunNormBeta data:
intersect = intersect(rownames(Mset.norm), rownames(filtered.p))
length(intersect)
Mset.fn.sub = Mset.norm[intersect,]
FunNormBetas.sub = getBeta(Mset.fn.sub)
RawBetas.sub = RawBetas[intersect,]

## Check density plots after excluding the poorly-detected probes:
png(file="BetaValue_Distributions_AfterQC.png",width=1400,height=700,pointsize=12)
par(mfrow=c(1,2))
densityPlot(RawBetas.sub, sampGroups = pd$Slide, main = "Raw Beta", xlab = "Beta")
densityPlot(FunNormBetas.sub, sampGroups = pd$Slide, main = "Normalized Beta", xlab = "Beta")
dev.off()

## The 'FunNormBetas.sub' object includes data for samples and probes with good detection p-values, and that have undergone functional normalization:
save(prop.table,detect.p,pd,FunNormBetas.sub,Mset.fn.sub, file = "Allison.Appleton.1_BasicQC_and_FunctionalNormalization.RData")

