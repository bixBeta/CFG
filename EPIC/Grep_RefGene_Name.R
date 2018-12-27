library("IlluminaHumanMethylationEPICanno.ilm10b3.hg19")
data(Other)

genes <- c("CNR1", "CNR2", "HDAC3", "FAAH", "MGLL", "CHRNB3", "CHRNA10",
	"CHRNA1", "CHRNA2", "CHRNA3", "CHRNA5","CHRNA6","CHRNA3","CHRNE",
	"CHRNA4","CHRNA7","CHRNA9","CHRNB1","CHRNB4","CHRNB2",
	"CHRND","CHRNG", "NR3C1", "HSD11B2")

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
	
		for(i in 1:24){
			index <- grepM(pattern = genes[i], varname = "UCSC_RefGene_Name", data=Other, exact = T, identifier = ";")
			index = as.data.frame(index)
			indices <- c(indices,index[,1])
		}

annot.sub <- Other[unique(indices),]

### rownames(annot.sub) = list of cpgs

