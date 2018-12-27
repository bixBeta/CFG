install.packages("RPMM")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("BiocManager")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("minfi", version = "3.8")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("minfiData", version = "3.8")


## Load Libraries
library(RPMM)
library(ggplot2)

## Load Data from previous step:
load("Allison.Appleton.1_BasicQC_and_FunctionalNormalization.RData")

## Load Annotation File (450K):
#load("450K_annotation_file.RData")
#rownames(annot) = as.character(annot$Name)
## Load Annotation File (EPIC):
annot = read.csv("MethylationEPIC_v-1-0_B2.csv",header=T)
## Drop duplicated rows
#index <- which(duplicated(annot$IlmnID))
#annot = annot[-index, ]
rownames(annot) = as.character(annot$IlmnID)
annot$Type = NULL
annot$Type = as.character(annot$Infinium_Design_Type)

##############################################################################################
#### (2) Use BMIQ to standardize across probe types ##########################################
##############################################################################################

## Check if Type I and Type II probe distributions are different:
## Join array data for 3 random samples with probe types:
X2 = as.data.frame(merge(FunNormBetas.sub[,c(19,50,190)],annot[,c("Type","Name")],by="row.names"))
X2$Type = as.character(X2$Type)
colnames(X2) = c("CGid","Sample1","Sample2","Sample3","ProbeType","Name")
## Generate density plots by probe type:
library(ggplot2)
png(file="ProbeTypes_After_FunNorm(S1).png",width=700,height=700,pointsize=12)
ggplot(X2, aes(x=Sample1, fill=ProbeType)) + geom_density(alpha=.3) + theme_bw()
dev.off()
png(file="ProbeTypes_After_FunNorm(S2).png",width=700,height=700,pointsize=12)
ggplot(X2, aes(x=Sample2, fill=ProbeType)) + geom_density(alpha=.3) + theme_bw()
dev.off()
png(file="ProbeTypes_After_FunNorm(S3).png",width=700,height=700,pointsize=12)
ggplot(X2, aes(x=Sample3, fill=ProbeType)) + geom_density(alpha=.3) + theme_bw()
dev.off()

#### Type I and II probes show different distributions, implement BMIQ:

## Read in BMIQ function:
BMIQ <- function(beta.v,design.v,nL=3,doH=TRUE,nfit=50000,th1.v=c(0.2,0.75),th2.v=NULL,niter=5,tol=0.001,plots=TRUE,sampleID=1){
  
  type1.idx <- which(design.v==1);
  type2.idx <- which(design.v==2);
  
  beta1.v <- beta.v[type1.idx];
  beta2.v <- beta.v[type2.idx];
  
  ### check if there are exact 0's or 1's. If so, regularise using minimum positive and maximum below 1 values.
  if(min(beta1.v)==0){
    beta1.v[beta1.v==0] <- min(setdiff(beta1.v,0));
  }
  if(min(beta2.v)==0){
    beta2.v[beta2.v==0] <- min(setdiff(beta2.v,0));
  }
  if(max(beta1.v)==1){
    beta1.v[beta1.v==1] <- max(setdiff(beta1.v,1));
  }
  if(max(beta2.v)==1){
    beta2.v[beta2.v==1] <- max(setdiff(beta2.v,1));
  }
  
  ### estimate initial weight matrix from type1 distribution
  w0.m <- matrix(0,nrow=length(beta1.v),ncol=nL);
  w0.m[which(beta1.v <= th1.v[1]),1] <- 1;
  w0.m[intersect(which(beta1.v > th1.v[1]),which(beta1.v <= th1.v[2])),2] <- 1;
  w0.m[which(beta1.v > th1.v[2]),3] <- 1;
  
  ### fit type1
  print("Fitting EM beta mixture to type1 probes");
  rand.idx <- sample(1:length(beta1.v),nfit,replace=FALSE)
  em1.o <- blc(matrix(beta1.v[rand.idx],ncol=1),w=w0.m[rand.idx,],maxiter=niter,tol=tol);
  subsetclass1.v <- apply(em1.o$w,1,which.max);
  subsetth1.v <- c(mean(c(max(beta1.v[rand.idx[subsetclass1.v==1]]),min(beta1.v[rand.idx[subsetclass1.v==2]]))),mean(c(max(beta1.v[rand.idx[subsetclass1.v==2]]),min(beta1.v[rand.idx[subsetclass1.v==3]]))));
  class1.v <- rep(2,length(beta1.v));
  class1.v[which(beta1.v < subsetth1.v[1])] <- 1;
  class1.v[which(beta1.v > subsetth1.v[2])] <- 3;
  nth1.v <- subsetth1.v;
  print("Done");
  
  ### generate plot from estimated mixture
  if(plots){
    print("Check");
    tmpL.v <- as.vector(rmultinom(1:nL,length(beta1.v),prob=em1.o$eta));
    tmpB.v <- vector();
    for(l in 1:nL){
      tmpB.v <- c(tmpB.v,rbeta(tmpL.v[l],em1.o$a[l,1],em1.o$b[l,1]));
    }
    
    pdf(paste("Type1fit-",sampleID,".pdf",sep=""),width=6,height=4);
    plot(density(beta1.v));
    d.o <- density(tmpB.v);
    points(d.o$x,d.o$y,col="green",type="l")
    legend(x=0.5,y=3,legend=c("obs","fit"),fill=c("black","green"),bty="n");
    dev.off();
  }
  
  ### Estimate Modes 
  d1U.o <- density(beta1.v[class1.v==1])
  d1M.o <- density(beta1.v[class1.v==3])
  mod1U <- d1U.o$x[which.max(d1U.o$y)]
  mod1M <- d1M.o$x[which.max(d1M.o$y)]
  d2U.o <- density(beta2.v[which(beta2.v<0.4)]);
  d2M.o <- density(beta2.v[which(beta2.v>0.6)]);
  mod2U <- d2U.o$x[which.max(d2U.o$y)]
  mod2M <- d2M.o$x[which.max(d2M.o$y)]
  
  ### now deal with type2 fit
  th2.v <- vector();
  th2.v[1] <- nth1.v[1] + (mod2U-mod1U);
  th2.v[2] <- nth1.v[2] + (mod2M-mod1M);
  
  ### estimate initial weight matrix 
  w0.m <- matrix(0,nrow=length(beta2.v),ncol=nL);
  w0.m[which(beta2.v <= th2.v[1]),1] <- 1;
  w0.m[intersect(which(beta2.v > th2.v[1]),which(beta2.v <= th2.v[2])),2] <- 1;
  w0.m[which(beta2.v > th2.v[2]),3] <- 1;
  
  print("Fitting EM beta mixture to type2 probes");
  rand.idx <- sample(1:length(beta1.v),nfit,replace=FALSE)
  em2.o <- blc(matrix(beta2.v[rand.idx],ncol=1),w=w0.m[rand.idx,],maxiter=niter,tol=tol);
  print("Done");
  
  ### for type II probes assign to state (unmethylated, hemi or full methylation)
  subsetclass2.v <- apply(em2.o$w,1,which.max);
  subsetth2.v <- c(mean(max(beta2.v[rand.idx[subsetclass2.v==1]]),min(beta2.v[rand.idx[subsetclass2.v==2]])),mean(max(beta2.v[rand.idx[subsetclass2.v==2]]),min(beta2.v[rand.idx[subsetclass2.v==3]])));
  class2.v <- rep(2,length(beta2.v));
  class2.v[which(beta2.v < subsetth2.v[1])] <- 1;
  class2.v[which(beta2.v > subsetth2.v[2])] <- 3;
  
  ### generate plot
  if(plots){
    tmpL.v <- as.vector(rmultinom(1:nL,length(beta2.v),prob=em2.o$eta));
    tmpB.v <- vector();
    for(lt in 1:nL){
      tmpB.v <- c(tmpB.v,rbeta(tmpL.v[lt],em2.o$a[lt,1],em2.o$b[lt,1]));
    }
    pdf(paste("Type2fit-",sampleID,".pdf",sep=""),width=6,height=4);
    plot(density(beta2.v));
    d.o <- density(tmpB.v);
    points(d.o$x,d.o$y,col="green",type="l")
    legend(x=0.5,y=3,legend=c("obs","fit"),fill=c("black","green"),bty="n");
    dev.off();
  }
  
  classAV1.v <- vector();classAV2.v <- vector();
  for(l in 1:nL){
    classAV1.v[l] <-  em1.o$mu[l,1];
    classAV2.v[l] <-  em2.o$mu[l,1];
  }
  
  ### start normalising type2 probes
  print("Start normalising type 2 probes");
  nbeta2.v <- beta2.v;
  ### select U probes
  lt <- 1;
  selU.idx <- which(class2.v==lt);
  selUR.idx <- selU.idx[which(beta2.v[selU.idx] > classAV2.v[lt])];
  selUL.idx <- selU.idx[which(beta2.v[selU.idx] < classAV2.v[lt])];
  ### find prob according to typeII distribution
  p.v <- pbeta(beta2.v[selUR.idx],em2.o$a[lt,1],em2.o$b[lt,1],lower.tail=FALSE);
  ### find corresponding quantile in type I distribution
  q.v <- qbeta(p.v,em1.o$a[lt,1],em1.o$b[lt,1],lower.tail=FALSE);
  nbeta2.v[selUR.idx] <- q.v;
  p.v <- pbeta(beta2.v[selUL.idx],em2.o$a[lt,1],em2.o$b[lt,1],lower.tail=TRUE);
  ### find corresponding quantile in type I distribution
  q.v <- qbeta(p.v,em1.o$a[lt,1],em1.o$b[lt,1],lower.tail=TRUE);
  nbeta2.v[selUL.idx] <- q.v;
  
  ### select M probes
  lt <- 3;
  selM.idx <- which(class2.v==lt);
  selMR.idx <- selM.idx[which(beta2.v[selM.idx] > classAV2.v[lt])];
  selML.idx <- selM.idx[which(beta2.v[selM.idx] < classAV2.v[lt])];
  ### find prob according to typeII distribution
  p.v <- pbeta(beta2.v[selMR.idx],em2.o$a[lt,1],em2.o$b[lt,1],lower.tail=FALSE);
  ### find corresponding quantile in type I distribution
  q.v <- qbeta(p.v,em1.o$a[lt,1],em1.o$b[lt,1],lower.tail=FALSE);
  nbeta2.v[selMR.idx] <- q.v;
  
  if(doH){ ### if TRUE also correct type2 hemimethylated probes
    ### select H probes and include ML probes (left ML tail is not well described by a beta-distribution).
    lt <- 2;
    selH.idx <- c(which(class2.v==lt),selML.idx);
    minH <- min(beta2.v[selH.idx])
    maxH <- max(beta2.v[selH.idx])
    deltaH <- maxH - minH;
    #### need to do some patching
    deltaUH <- -max(beta2.v[selU.idx]) + min(beta2.v[selH.idx])
    deltaHM <- -max(beta2.v[selH.idx]) + min(beta2.v[selMR.idx])
    
    ## new maximum of H probes should be
    nmaxH <- min(nbeta2.v[selMR.idx]) - deltaHM;
    ## new minimum of H probes should be
    nminH <- max(nbeta2.v[selU.idx]) + deltaUH;
    ndeltaH <- nmaxH - nminH;
    
    ### perform conformal transformation (shift+dilation)
    ## new_beta_H(i) = a + hf*(beta_H(i)-minH);
    hf <- ndeltaH/deltaH ;
    ### fix lower point first
    nbeta2.v[selH.idx] <- nminH + hf*(beta2.v[selH.idx]-minH);
    
  }
  
  pnbeta.v <- beta.v;
  pnbeta.v[type1.idx] <- beta1.v;
  pnbeta.v[type2.idx] <- nbeta2.v;
  
  ### generate final plot to check normalisation
  if(plots){
    print("Generating final plot");
    d1.o <- density(beta1.v);
    d2.o <- density(beta2.v);
    d2n.o <- density(nbeta2.v);
    ymax <- max(d2.o$y,d1.o$y,d2n.o$y);
    pdf(paste("CheckBMIQ-",sampleID,".pdf",sep=""),width=6,height=4)
    plot(density(beta2.v),type="l",ylim=c(0,ymax),xlim=c(0,1));
    points(d1.o$x,d1.o$y,col="red",type="l");
    points(d2n.o$x,d2n.o$y,col="blue",type="l");
    legend(x=0.5,y=ymax,legend=c("type1","type2","type2-BMIQ"),bty="n",fill=c("red","black","blue"));
    dev.off();
  }
  
  print(paste("Finished for sample ",sampleID,sep=""));
  
  return(list(nbeta=pnbeta.v,class1=class1.v,class2=class2.v,av1=classAV1.v,av2=classAV2.v,hf=hf,th1=nth1.v,th2=th2.v));
  
}

CheckBMIQ <- function(beta.v,design.v,pnbeta.v){### pnbeta is BMIQ normalised profile
  
  type1.idx <- which(design.v==1);
  type2.idx <- which(design.v==2);
  
  beta1.v <- beta.v[type1.idx];
  beta2.v <- beta.v[type2.idx];
  pnbeta2.v <- pnbeta.v[type2.idx];
  
}

## Run BMIQ function:
set.seed(123432)
pt = as.data.frame(annot[rownames(FunNormBetas.sub),c("Name","Type")])
pt$ProbeType = ifelse(pt$Type %in% "I",1,2)
FunNorm.BMIQ = apply(FunNormBetas.sub[,1:length(colnames(FunNormBetas.sub))],2,function(a) BMIQ(a,pt$ProbeType,plots=FALSE)$nbeta )


## The 'FunNorm.BMIQ' object includes beta values that have undergone BMIQ probe-type adjustment:
save(pd,FunNorm.BMIQ, file = "Allison.Appleton.0036_2_BMIQ_Adjsuted.RData")

#### Check whether BMIQ corrects probe-type bias:
## Join array data with probe types:
X2 = as.data.frame(merge(FunNorm.BMIQ[,c(9,20,90)],annot[,c("Type","Name")],by="row.names"))
X2$Type = as.character(X2$Type)
colnames(X2) = c("CGid","Sample1","Sample2","Sample3","ProbeType","Name")
## Generate density plots by probe type:
png(file="ProbeTypes_After_FunNorm_BMIQ(S1).png",width=700,height=700,pointsize=12)
ggplot(X2, aes(x=Sample1, fill=ProbeType)) + geom_density(alpha=.3) + theme_bw()
dev.off()
png(file="ProbeTypes_After_FunNorm_BMIQ(S2).png",width=700,height=700,pointsize=12)
ggplot(X2, aes(x=Sample2, fill=ProbeType)) + geom_density(alpha=.3) + theme_bw()
dev.off()
png(file="ProbeTypes_After_FunNorm_BMIQ(S3).png",width=700,height=700,pointsize=12)
ggplot(X2, aes(x=Sample3, fill=ProbeType)) + geom_density(alpha=.3) + theme_bw()
dev.off()

#### Check Density Plots after BMIQ:
png(file="BetaValue_Distributions_AfterQC_FunNorm_BMIQ.png",width=700,height=700,pointsize=12)
densityPlot(FunNorm.BMIQ, sampGroups = pd$Slide, main = "Normalized Beta", xlab = "Beta")
dev.off()
