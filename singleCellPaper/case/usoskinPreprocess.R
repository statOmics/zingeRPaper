library(openxlsx)
library(Biobase)
data <- openxlsx::read.xlsx(xlsxFile="/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/usoskin/Usoskin et al. External resources Table 1.xlsx", sheet=2, rows=12:25345, cols=10:873, colNames=FALSE)
#save(data,file="usoskinData.RData")
#library(tools)
#resaveRdaFiles("usoskinData.RData")

## pheno Data
pData <- openxlsx::read.xlsx(xlsxFile="/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/usoskin/Usoskin et al. External resources Table 1.xlsx", sheet=2, rows=1:10, cols=9:873, colNames=FALSE)
pData <- t(pData)
colnames(pData) <- pData[1,]
pData=pData[-1,]
rownames(pData) <- colnames(data)
pData <- AnnotatedDataFrame(as.data.frame(pData))

## feature Data
fData <- openxlsx::read.xlsx(xlsxFile="/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/usoskin/Usoskin et al. External resources Table 1.xlsx", sheet=2, rows=12:25345, cols=1, colNames=FALSE)
fData <- AnnotatedDataFrame(fData)


eset <- ExpressionSet(assayData=as.matrix(data), phenoData=pData, featureData=fData)

## only use single-cell samples
eset <- eset[,pData(eset)$Content=="cell"]

## only use neuronal cells
sum(pData(eset)[,7] == "NF" | pData(eset)[,7] == "NP" | pData(eset)[,7] == "PEP" | pData(eset)[,7] == "TH") # the 622 cells used in the paper
eset= eset[,(pData(eset)[,7] == "NF" | pData(eset)[,7] == "NP" | pData(eset)[,7] == "PEP" | pData(eset)[,7] == "TH")]

## keep in cpm for MAST analysis
esetCpm=eset
save(esetCpm,file="~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/case/esetUsoskinCpm.RData")


## filter genes with no expression, convert back to counts
exprs(eset) = round(sweep(exprs(eset),2,STATS=as.numeric(as.character(pData(eset)$Reads)),FUN="*")/1e6)
eset=eset[rowSums(exprs(eset))>0,]
save(eset,file="~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/case/esetUsoskin.RData")




