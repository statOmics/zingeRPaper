library(Biobase)
load("/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/usoskin/esetUsoskin.RData")
## only get cells from picking session Cold
eset=eset[,pData(eset)[,2]=="Cold"]
## celltype TH
eset=eset[,pData(eset)[,7]=="TH"]
cellType=factor(pData(eset)[,7], levels="TH")
## remove one sample from different lib
eset=eset[,!pData(eset)[,4]=="L281"]
l226Cells=which(pData(eset)[,4]=="L226")
l227Cells=which(pData(eset)[,4]=="L227")
l228Cells=which(pData(eset)[,4]=="L228")


#three libraries. 16v16 for L226, 11v12 for L227, 18v17 for L228, so 45v45.
subsetMatrixUsoskin <- matrix(NA,nrow=30,ncol=45+45)
set.seed(99)
for(i in 1:30){
    l226Condition1 <- sample(l226Cells,16)
    l227Condition1 <- sample(l227Cells,11)
    l228Condition1 <- sample(l228Cells,18)
    l226Condition2 <- l226Cells[!l226Cells%in%l226Condition1]
    l227Condition2 <- l227Cells[!l227Cells%in%l227Condition1]
    l228Condition2 <- l228Cells[!l228Cells%in%l228Condition1]
  subsetMatrixUsoskin[i,] <- c(l226Condition1,l227Condition1,l228Condition1,
			       l226Condition2,l227Condition2,l228Condition2)
}

write.table(subsetMatrixUsoskin,file="~/Dropbox/PhD/Research/singleCell/usoskin/FPR/subsetMatrixUsoskinFPR.txt",row.names=FALSE,col.names = FALSE)



