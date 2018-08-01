# import count tables ($ paste A_D_1.raw.counts A_D_2.raw.counts A_D_3.raw.counts A_N_1.raw.counts A_N_2.raw.counts A_N_3.raw.counts p_D_1.raw.counts p_D_2.raw.counts p_D_3.raw.counts p_N_1.raw.counts p_N_2.raw.counts p_N_3.raw.counts WT_D_1.raw.counts WT_D_2.raw.counts WT_D_3.raw.counts WT_N_1.raw.counts WT_N_2.raw.counts WT_N_3.raw.counts > Experiment.Name.raw.COUNTS.txt)


outputPrefix <- "Experiment.Name"

## -------------------------------------------------------------------------------------------------------------------
# sampleTable prep  <-  phenoTable

sampleFiles <- c("A_D_1.raw.counts", 
                 "A_D_2.raw.counts",
                 "A_D_3.raw.counts",
                 "A_N_1.raw.counts", 
                 "A_N_2.raw.counts", 
                 "A_N_3.raw.counts",
                 "p_D_1.raw.counts", 
                 "p_D_2.raw.counts", 
                 "p_D_3.raw.counts",
                 "p_N_1.raw.counts", 
                 "p_N_2.raw.counts", 
                 "p_N_3.raw.counts", 
                 "WT_D_1.raw.counts", 
                 "WT_D_2.raw.counts", 
                 "WT_D_3.raw.counts", 
                 "WT_N_1.raw.counts", 
                 "WT_N_2.raw.counts", 
                 "WT_N_3.raw.counts")


sampleNames     <- c("AD1",
                     "AD2",
                     "AD3",
                     "AN1",
                     "AN2",
                     "AN3",
                     "PD1",
                     "PD2",
                     "PD3",
                     "PN1",
                     "PN2",
                     "PN3",
                     "WD1",
                     "WD2",
                     "WD3",
                     "WN1",
                     "WN2",
                     "WN3")

sampleCondition <- c("AD",
                     "AD",
                     "AD",
                     "ANU",
                     "ANU",
                     "ANU",
                     "PD",
                     "PD",
                     "PD",
                     "PN",
                     "PN",
                     "PN",
                     "WD",
                     "WD",
                     "WD",
                     "WN",
                     "WN",
                     "WN")

sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)

## -------------------------------------------------------------------------------------------------------------------
# adding headers 
# import Experiment.Name.0026.raw.COUNTS as rawCounts first

n <- (0:18)*2
n[1] <- 1
countTable <- rawCounts[,n]
colnames(countTable) <- c("gene",
                          "AD1",
                          "AD2",
                          "AD3",
                          "AN1",
                          "AN2",
                          "AN3",
                          "PD1",
                          "PD2",
                          "PD3",
                          "PN1",
                          "PN2",
                          "PN3",
                          "WD1",
                          "WD2",
                          "WD3",
                          "WN1",
                          "WN2",
                          "WN3")
row.names(countTable) <- countTable$gene
countTable<- countTable[,2:19]
sampleTable$sampleName == colnames(countTable)


# FOR HTSEQ  -- CHECK 
#grep("__",x = row.names(countTable))
#countTable[grep("__",x = row.names(countTable)),] # rm if TRUE

#countTable <- countTable[1:26364,]
#countTable[grep("__",x = row.names(countTable)),] # verify

## -------------------------------------------------------------------------------------------------------------------
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = countTable,
                              colData = sampleTable,
                              design = ~ condition)

treatments = c("AD", "ANU", "PD", "PN", "WD", "WN") # Treatments of interest
dds$condition <- factor(colData(dds)$condition,
                        levels = treatments)

## -------------------------------------------------------------------------------------------------------------------
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)
res
resultsNames(dds)
mcols(res, use.names=TRUE) # brief description of res headers

## -------------------------------------------------------------------------------------------------------------------
WD.WN <- results(dds, contrast=c("condition","WD","WN"))
pd.pn <- results(dds, contrast=c("condition","PD","PN"))
Ad.An <- results(dds, contrast=c("condition","AD","ANU"))
WD.pd <- results(dds, contrast=c("condition","WD","PD"))
WD.Ad <- results(dds, contrast=c("condition","WD","AD"))
pd.Ad <- results(dds, contrast=c("condition","PD","AD"))
WN.pn <- results(dds, contrast=c("condition","WN","PN"))
WN.An <- results(dds, contrast=c("condition","WN","ANU"))
pn.An <- results(dds, contrast=c("condition","PN","ANU"))

## -------------------------------------------------------------------------------------------------------------------
# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

# save normalized values
write.table(as.data.frame(assay(rld)),file = paste0(outputPrefix, ".rlog-transformed-counts.txt"), sep = '\t')
write.table(as.data.frame(assay(vsd)),file = paste0(outputPrefix, ".vst-transformed-counts.txt"), sep = '\t')

# saving deseq results as .csv
write.csv(WD.WN, file = paste0(outputPrefix, ".WD.WN.results.csv"))
write.csv(pd.pn, file = paste0(outputPrefix, ".pd.pn.results.csv"))
write.csv(Ad.An, file = paste0(outputPrefix, ".Ad.An.results.csv"))
write.csv(WD.pd, file = paste0(outputPrefix, ".WD.pd.results.csv"))
write.csv(WD.Ad, file = paste0(outputPrefix, ".WD.Ad.results.csv"))
write.csv(pd.Ad, file = paste0(outputPrefix, ".pd.Ad.results.csv"))
write.csv(WN.pn, file = paste0(outputPrefix, ".WN.pn.results.csv"))
write.csv(WN.An, file = paste0(outputPrefix, ".WN.An.results.csv"))
write.csv(pn.An, file = paste0(outputPrefix, ".pn.An.results.csv"))


save(countTable, sampleTable, dds, rawCounts, rld, vsd, WD.WN, pd.pn, Ad.An, WD.pd, 
     WD.Ad, pd.Ad, WN.pn, WN.An, pn.An, file = "Experiment.Name.0026.Rdata")
