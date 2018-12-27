library("IlluminaHumanMethylationEPICanno.ilm10b3.hg19")
data(Other)

genes <- c("NR3C1",
           "HSD11B2",
           "FKBP5",
           "SLC6A4",
           "HTR2A",
           "NET",
           "CRH",
           "BDNF",
           "LEP",
           "LEPR",
           "LEPROTL",
           "H19",
           "IGF1",
           "IGF2")

##################################################################################

grepM = function(pattern, varname, data, exact = T, identifier = ";") 
	{
		if(exact) {
			x = data[,varname]
			tmp = grep(pattern, x)
			x.sub = data[tmp,]
			s = apply(as.matrix(x.sub[,varname]), 1, function(x) strsplit(x, identifier))
			ind = lapply(s, function(x) sum(pattern == x[[1]])>=1)
			ind.1 = unlist(ind)
			tmp.1 = tmp[ind.1]
			tmp.1
			
		}
		else { 
			tmp
		}
	}

# Illustration of the function 
#indices = grepM(pattern = "BCL2", varname = "UCSC_RefGene_Name", data = Other, 
#                exact = T, identifier = ";")
#Other[indices,]


##Loop to pull all indices from genes object for Pot Genes (from genes)

		indices=NULL	
	
		for(i in 1:14){
			index <- grepM(pattern = genes[i], varname = "UCSC_RefGene_Name", data=Other, exact = T, identifier = ";")
			index = as.data.frame(index)
			indices <- c(indices,index[,1])
		}  

annot.sub <- Other[unique(indices),]

### rownames(annot.sub) = list of cpgs
cpGs.genes <- annot.sub[3]

write.csv(cpGs.genes, "0036.cpGs.genes.csv")

load("Allison.Appleton.0036_4_Final_NHBCS_EPIC_Data.RData")
subset <- Epic.NHBCS[rownames(Epic.NHBCS) %in% rownames(annot.sub), ]
write.csv(subset, "Allison.Appleton.0036.geneSubset.csv")


