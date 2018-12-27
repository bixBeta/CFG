##############################################################################################
#### (3) Use ComBat to remove batch effect ###################################################
##############################################################################################

## Evaluate beta value distributions in adjusted data:
png(file="BetaValue_Distribution_After_FunNorm_BMIQ.png",width=700,height=700,pointsize=22)
densityPlot(FunNorm.BMIQ, sampGroups = pd$Sample_Plate, main = "ComBat/FunNorm Betas", xlab = "Beta")
dev.off()

## Read in a function that will calculate the variance of each row
rowVars <- function(x, na.rm=FALSE, dims=1, unbiased=TRUE,
                    SumSquares=FALSE, twopass=FALSE) {
  if (SumSquares) return(rowSums(x^2, na.rm, dims))
  N <- rowSums(!is.na(x), FALSE, dims)
  Nm1 <- if (unbiased) N-1 else N
  if (twopass) {x <- if (dims==0) x - mean(x, na.rm=na.rm) else
    sweep(x, 1:dims, rowMeans(x,na.rm,dims))}
  (rowSums(x^2, na.rm, dims) - rowSums(x, na.rm, dims)^2/N) / Nm1
}

## ComBat assumes normally distributed data, so convert to m-values:
mval <- apply(FunNorm.BMIQ, 2, function(x) log2((x)/(1-x)))
save(mval,pd, file = "Allison.Appleton.0036_3_(mvals)ComBat_Adjusted.RData")

## Calculate the variance of each probe and remove any with a variance of 0 prior to Combat.
vars = as.matrix(rowVars(mval))

## Replace all probes with no variance with NA and remove them from the FunNorm set
vars[vars == 0] = NA
vars = na.omit(vars)
intersect = intersect(rownames(vars), rownames(mval))
print(length(intersect))
fn.sub = FunNorm.BMIQ[intersect,]

## Recalculate m-values from FunNorm set and get sample sheet data parsed into "pheno" object.
## Create a model matrix with the desired aspect to adjust for.
mval <- apply(fn.sub, 2, function(x) log2((x)/(1-x)))
pheno <- pd[colnames(mval),]
table(ifelse(rownames(pheno) == colnames(mval),"Match","Off"))

## Group small batches with others:
pheno$Batch <- NA
pheno$Batch <- pheno$Sample_Plate
batch <- pheno$Batch
mod <- model.matrix(~1, data=pheno)

###########################################################################
## Check variation in array data associated with batch (ie. Slide/plate/box)
###########################################################################

## Run a principle component analysis to determine if there are any remaining
## batch effects following data normalization.

PCobj = prcomp(t(mval), retx = T, center = T, scale. = T)

# Extract the PCs from the PCobj object

PCs = PCobj$x

# Extract the proportion of variability and cummulative proportion of 
# varibility explained by the top R PCs.

R = 5
propvar = summary(PCobj)$importance["Proportion of Variance", 1:R]
cummvar = summary(PCobj)$importance["Cumulative Proportion", 1:R]

## Generate plots of the resulting PCs

# Plot of the proportion of variability explained by the top R PCs
# Plot of the cummulative proportion of variability explained by the top R PCs

png(file="Variation_Explained_by_PCs.png",width=1400,height=700,pointsize=12)
par(mfrow=c(1,2))	
par(mar = c(5,5,4,2))
barplot(propvar*100, xlab = paste("Top", R, "PCs", sep = " "), ylab = "Variation Explained (%)", cex.axis = 1.5, cex.lab = 1.8, cex.names = 1.5)
par(mar = c(5,5,4,2))
barplot(cummvar*100, xlab = paste("Top", R, "PCs", sep = " "), ylab = "Cummulative Variation Explained (%)",cex.axis = 1.5, cex.lab = 1.8, cex.names = 1.5)
abline(a = 100*cummvar[R], b = 0, lwd = 3, col = "red", lty = "dashed")
dev.off()

# Plot of PCX and PCY; by slide number

PCs = PCobj$x
PCs =PCs[,1:5]
Prin.comp<-merge(PCs,pheno, by = "row.names",all=T) 

png(file="PC_Variation_by_Batch.png",width=900,height=900,pointsize=12)
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Batch), xlab = "PC1", ylab = "PC2")
legend("topright", legend=levels(as.factor(Prin.comp$Batch)), text.col=seq_along(levels(as.factor(Prin.comp$Batch))))
dev.off()





install.packages("devtools")
library(devtools)
install_github("Bioconductor-mirror/sva")









## Run ComBat to remove slide # batch effects
library(sva)
library(combat)
combat.adj = ComBat(mval,batch = batch, mod = mod)

## Check to see if batch effect was succesfully removed

PCobj = prcomp(t(combat.adj), retx = T, center = T, scale. = T)
PCs = PCobj$x
R = 5
propvar = summary(PCobj)$importance["Proportion of Variance", 1:R]
cummvar = summary(PCobj)$importance["Cumulative Proportion", 1:R]

png(file="Variation_Explained_by_PCs_AfterComBat.png",width=1400,height=700,pointsize=12)
par(mfrow=c(1,2))	
par(mar = c(5,5,4,2))
barplot(propvar*100, xlab = paste("Top", R, "PCs", sep = " "), ylab = "Variation Explained (%)", cex.axis = 1.5, cex.lab = 1.8, cex.names = 1.5)
par(mar = c(5,5,4,2))
barplot(cummvar*100, xlab = paste("Top", R, "PCs", sep = " "), ylab = "Cummulative Variation Explained (%)",cex.axis = 1.5, cex.lab = 1.8, cex.names = 1.5)
abline(a = 100*cummvar[R], b = 0, lwd = 3, col = "red", lty = "dashed")
dev.off()

PCs = PCobj$x
PCs =PCs[,1:5]
Prin.comp<-merge(PCs,pheno, by = "row.names",all=T) 

#### Check whether batches are still distinguished by first and second PC:
png(file="PC_Variation_by_Batch_AfterComBat.png",width=900,height=900,pointsize=12)
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Batch), xlab = "PC1", ylab = "PC2")
legend("topright", legend=levels(as.factor(Prin.comp$Batch)), text.col=seq_along(levels(as.factor(Prin.comp$Batch))))
dev.off()

###########################################################################
## Convert the adjusted values back into betas

expit2 = function(x) 2^x/(1+2^x)
fnbetas.adj = expit2(combat.adj)

## Save normalized and batch-adjusted beta values
save(fnbetas.adj,pd, file = "Allison.Appleton.0036_3_ComBat_Adjusted.RData")

table(ifelse(rownames(pd) == colnames(fnbetas.adj),"Match","Off"))

## Check how effective the normalization was with a density plot
png(file="BetaValue_Distribution_After_FunNorm_BMIQ_ComBat.png",width=700,height=700,pointsize=22)
densityPlot(fnbetas.adj, sampGroups = pd$Sample_Plate, main = "ComBat/FunNorm Betas", xlab = "Beta")
dev.off() 




















##############################################################################################
#### (Pre) Load relevant packages and functions ##############################################
# ##############################################################################################
# 
# #### Create Functions:
# `%notin%` <- function(x,y) !(x %in% y) 
# 
# ## Load Data
# load("3_ComBat_Adjusted.RData")
# length(match(rownames(pd),colnames(fnbetas.adj)))
# colnames(fnbetas.adj) = substr(as.character(pd$Sample_Name),1,4)
# 
# ## Load Annotation File (450K):
# load("C://Users//Todd//Desktop//450K Annotations//450K_annotation_file.RData")
# rownames(annot) = as.character(annot$Name)
# 
# ## Load Libraries:
library(ggplot2)
# library(minfi)

##############################################################################################
#### (1) Check data for obvious problems #####################################################
##############################################################################################

## Were any missing values generated during processing?
table(is.na(fnbetas.adj))  #All should be FALSE

annot = read.csv("rmDUP_MethylationEPIC_v-1-0_B2.csv",header=T)
rownames(annot) = as.character(annot$Name)
## Are distributions of different probe types aligned (Compare random set of samples)?
ran.sel = round(runif(3,1,length(colnames(fnbetas.adj))),digits=0)
X2 = as.data.frame(merge(fnbetas.adj[,ran.sel],annot[,c("Infinium_Design_Type","Name")],by="row.names"))
X2$Infinium_Design_Type = as.character(X2$Infinium_Design_Type)
colnames(X2) = c("CGid","Sample1","Sample2","Sample3","ProbeType","Name")
## Generate density plots by probe type:
png(file="ProbeTypes_Final(S1).png",width=700,height=700,pointsize=12)
ggplot(X2, aes(x=Sample1, fill=ProbeType)) + geom_density(alpha=.3) + theme_bw()
dev.off()

png(file="ProbeTypes_Final(S2).png",width=700,height=700,pointsize=12)
ggplot(X2, aes(x=Sample2, fill=ProbeType)) + geom_density(alpha=.3) + theme_bw()
dev.off()
png(file="ProbeTypes_Final(S3).png",width=700,height=700,pointsize=12)
ggplot(X2, aes(x=Sample3, fill=ProbeType)) + geom_density(alpha=.3) + theme_bw()
dev.off()

## Are distributions of different probe types aligned (Compare medians of all CGs)?
X3 = as.data.frame(merge(data.frame(CGid = rownames(fnbetas.adj),SampleMeans = rowMeans(fnbetas.adj)),annot[,c("Infinium_Design_Type","Name")],by="row.names"))
X3 = X3[,-1]
colnames(X3) = c("CGid","SampleMeans","ProbeType","Name")
png(file="ProbeTypes_Final(MeanAcrossSamples).png",width=700,height=700,pointsize=12)
ggplot(X3, aes(x=SampleMeans, fill=ProbeType)) + geom_density(alpha=.3) + theme_bw()
dev.off()

##############################################################################################
#### (2) Check for correlation between reps #####################################################################
##############################################################################################

# ## Create dataset that only includes replicates:
# X3 = as.data.frame(fnbetas.adj[,which(colnames(fnbetas.adj) %in% c("200516380035_R05C01","200516380035_R06C01","200516380035_R07C01"))])
# colnames(X3) = c("1198A","1206A","1468A","1198B","1206B","1468B")
# 
# ## Evaluate overall consistency across replicates:
# png(file="C://Users//Todd//Desktop//EPIC Data//NHBCS//Replicates(All CpGs).png",width=2100,height=700,pointsize=22)
# par(mfrow=c(1,3))
# cp = cor.test(X3[,"1198A"],X3[,"1198B"],method="spearman")
# plot(X3[,"1198A"],X3[,"1198B"],xlab="S1-A",ylab="S1-B",main=paste("Replicates (S1)","rho =",round(cp$estimate,digits=2),"P-value",ifelse(cp$p.value < 0.0001,"< 0.0001",paste("=",round(cp$p.value,digits=4))))  )
# cp = cor.test(X3[,"1206A"],X3[,"1206B"],method="spearman")
# plot(X3[,"1206A"],X3[,"1206B"],xlab="S2-A",ylab="S2-B",main=paste("Replicates (S2)","rho =",round(cp$estimate,digits=2),"P-value",ifelse(cp$p.value < 0.0001,"< 0.0001",paste("=",round(cp$p.value,digits=4))))  )
# cp = cor.test(X3[,"1468A"],X3[,"1468B"],method="spearman")
# plot(X3[,"1468A"],X3[,"1468B"],xlab="S3-A",ylab="S3-B",main=paste("Replicates (S3)","rho =",round(cp$estimate,digits=2),"P-value",ifelse(cp$p.value < 0.0001,"< 0.0001",paste("=",round(cp$p.value,digits=4))))  )
# dev.off()
# 
# ## Create dataset that only includes replicates with low methylation levels (0-40%):
# X3 = data.frame(fnbetas.adj[,which(colnames(fnbetas.adj) %in% c("1198","1206","1468"))],CGmax = apply(fnbetas.adj,1,function(qt) max(qt)))
# colnames(X3) = c("1198A","1206A","1468A","1198B","1206B","1468B","CGmax")
# X3 = X3[which(X3$CGmax < 0.4),]
# 
# ## Evaluate overall consistency across replicates:
# png(file="C://Users//Todd//Desktop//EPIC Data//NHBCS//Replicates(Below 0.4 CpGs).png",width=2100,height=700,pointsize=22)
# par(mfrow=c(1,3))
# cp = cor.test(X3[,"1198A"],X3[,"1198B"],method="spearman")
# plot(X3[,"1198A"],X3[,"1198B"],xlab="S1-A",ylab="S1-B",main=paste("Replicates (S1)","rho =",round(cp$estimate,digits=2),"P-value",ifelse(cp$p.value < 0.0001,"< 0.0001",paste("=",round(cp$p.value,digits=4))))  )
# cp = cor.test(X3[,"1206A"],X3[,"1206B"],method="spearman")
# plot(X3[,"1206A"],X3[,"1206B"],xlab="S2-A",ylab="S2-B",main=paste("Replicates (S2)","rho =",round(cp$estimate,digits=2),"P-value",ifelse(cp$p.value < 0.0001,"< 0.0001",paste("=",round(cp$p.value,digits=4))))  )
# cp = cor.test(X3[,"1468A"],X3[,"1468B"],method="spearman")
# plot(X3[,"1468A"],X3[,"1468B"],xlab="S3-A",ylab="S3-B",main=paste("Replicates (S3)","rho =",round(cp$estimate,digits=2),"P-value",ifelse(cp$p.value < 0.0001,"< 0.0001",paste("=",round(cp$p.value,digits=4))))  )
# dev.off()
# 
# ## Create dataset that only includes replicates with low methylation levels (0-4%):
# X3 = data.frame(fnbetas.adj[,which(colnames(fnbetas.adj) %in% c("1198","1206","1468"))],CGmax = apply(fnbetas.adj,1,function(qt) max(qt)))
# colnames(X3) = c("1198A","1206A","1468A","1198B","1206B","1468B","CGmax")
# X3 = X3[which(X3$CGmax < 0.04),]
# 
# ## Evaluate overall consistency across replicates:
# png(file="C://Users//Todd//Desktop//EPIC Data//NHBCS//Replicates(LowMethy CpGs).png",width=2100,height=700,pointsize=22)
# par(mfrow=c(1,3))
# cp = cor.test(X3[,"1198A"],X3[,"1198B"],method="spearman")
# plot(X3[,"1198A"],X3[,"1198B"],xlab="S1-A",ylab="S1-B",main=paste("Replicates (S1)","rho =",round(cp$estimate,digits=2),"P-value",ifelse(cp$p.value < 0.0001,"< 0.0001",paste("=",round(cp$p.value,digits=4))))  )
# cp = cor.test(X3[,"1206A"],X3[,"1206B"],method="spearman")
# plot(X3[,"1206A"],X3[,"1206B"],xlab="S2-A",ylab="S2-B",main=paste("Replicates (S2)","rho =",round(cp$estimate,digits=2),"P-value",ifelse(cp$p.value < 0.0001,"< 0.0001",paste("=",round(cp$p.value,digits=4))))  )
# cp = cor.test(X3[,"1468A"],X3[,"1468B"],method="spearman")
# plot(X3[,"1468A"],X3[,"1468B"],xlab="S3-A",ylab="S3-B",main=paste("Replicates (S3)","rho =",round(cp$estimate,digits=2),"P-value",ifelse(cp$p.value < 0.0001,"< 0.0001",paste("=",round(cp$p.value,digits=4))))  )
# dev.off()

##############################################################################################
#### (3) Save Final Data #####################################################################
##############################################################################################

# ## Drop duplicates then save Final Data:
# Epic.NHBCS = fnbetas.adj[,c(-283:-285)]
# pd = pd[c(-283:-285),]
# save(Epic.NHBCS,pd, file = "C://Users//Todd//Desktop//EPIC Data//NHBCS//4_Final_NHBCS_EPIC_Data.RData")

# Epic.NHBCS = fnbetas.adj[,c(-283:-285)]
# pd = pd[c(-283:-285),]
# save(Epic.NHBCS,pd, file = "C://Users//Todd//Desktop//EPIC Data//NHBCS//4_Final_NHBCS_EPIC_Data.RData")


Epic.NHBCS = fnbetas.adj
save(Epic.NHBCS,pd, file = "Allison.Appleton.0036_4_Final_NHBCS_EPIC_Data.RData")


