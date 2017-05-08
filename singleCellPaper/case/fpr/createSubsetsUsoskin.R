library(Biobase)
load("/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/usoskin/esetUsoskin.RData")


### TH celltype, picking session Cold.
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


##### NP celltype, 15 samples from each picking session. ignore the libraries.
library(Biobase)
load("/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/usoskin/esetUsoskin.RData")
## celltype TH
eset=eset[,pData(eset)[,7]=="NP"]
coldCells <- which(pData(eset)[,2]=="Cold")
rt1Cells <- which(pData(eset)[,2]=="RT-1")
rt2Cells <- which(pData(eset)[,2]=="RT-2")


subsetMatrixUsoskin <- matrix(NA,nrow=30,ncol=45+45)
set.seed(99)
for(i in 1:30){
    coldCondition1 <- sample(coldCells,15)
    rt1Condition1 <- sample(rt1Cells,15)
    rt2Condition1 <- sample(rt2Cells,15)
    coldCondition2 <- sample(coldCells[!coldCells%in%coldCondition1],15)
    rt1Condition2 <- sample(rt1Cells[!rt1Cells%in%rt1Condition1],15)
    rt2Condition2 <- sample(rt2Cells[!rt2Cells%in%rt2Condition1],15)
  subsetMatrixUsoskin[i,] <- c(coldCondition1,rt1Condition1,rt2Condition1,
			       coldCondition2,rt1Condition2,rt2Condition2)
}
write.table(subsetMatrixUsoskin,file="~/Dropbox/PhD/Research/singleCell/usoskin/FPR/subsetMatrixUsoskinFPR_cellTypeNP.txt",row.names=FALSE,col.names = FALSE)


##### two celltypes (NP and TH), three picking sessions in every condition.
coldNPCells <- which(pData(eset)[,7]=="NP" & pData(eset)[,2]=="Cold")
coldTHCells <- which(pData(eset)[,7]=="TH" & pData(eset)[,2]=="Cold")
rt1NPCells <- which(pData(eset)[,7]=="NP" & pData(eset)[,2]=="RT-1")
rt1THCells <- which(pData(eset)[,7]=="TH" & pData(eset)[,2]=="RT-1")
rt2NPCells <- which(pData(eset)[,7]=="NP" & pData(eset)[,2]=="RT-2")
rt2THCells <- which(pData(eset)[,7]=="TH" & pData(eset)[,2]=="RT-2")


subsetMatrixUsoskinTHNP <- matrix(NA,nrow=30,ncol=60+60)
set.seed(99)
for(i in 1:30){
    set.seed(i)
    coldNPCondition1 <- sample(coldNPCells,10)
    coldNPCondition2 <- sample(coldNPCells[!coldNPCells %in% coldNPCondition1],10)
    rt1NPCondition1 <- sample(rt1NPCells,10)
    rt1NPCondition2 <- sample(rt1NPCells[!rt1NPCells %in% rt1NPCondition1],10)
    rt2NPCondition1 <- sample(rt2NPCells,10)
    rt2NPCondition2 <- sample(rt2NPCells[!rt2NPCells %in% rt2NPCondition1],10)

    coldTHCondition1 <- sample(coldTHCells,10)
    coldTHCondition2 <- sample(coldTHCells[!coldTHCells %in% coldTHCondition1], 10)
    rt1THCondition1 <- sample(rt1THCells,10)
    rt1THCondition2 <- sample(rt1THCells[!rt1THCells %in% rt1THCondition1], 10)
    rt2THCondition1 <- sample(rt2THCells,10)
    rt2THCondition2 <- sample(rt2THCells[!rt2THCells %in% rt2THCondition1], 10)

  subsetMatrixUsoskinTHNP[i,] <- c(coldNPCondition1, rt1NPCondition1, rt2NPCondition1, coldTHCondition1, rt1THCondition1, rt2THCondition1,
    coldNPCondition2, rt1NPCondition2, rt2NPCondition2, coldTHCondition2, rt1THCondition2, rt2THCondition2)
}
write.table(subsetMatrixUsoskinTHNP,file="~/Dropbox/PhD/Research/singleCell/usoskin/FPR/subsetMatrixUsoskinFPR_cellTypesTHNP.txt",row.names=FALSE,col.names = FALSE)

### highest variability: random pick 15 cells from each picking session across all celltypes in each mock condition.
library(Biobase)
load("/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/usoskin/esetUsoskin.RData")
coldCells <- which(pData(eset)[,2]=="Cold")
rt1Cells <- which(pData(eset)[,2]=="RT-1")
rt2Cells <- which(pData(eset)[,2]=="RT-2")
subsetMatrixUsoskinRandomCellType <- matrix(NA,nrow=30,ncol=45+45)

for(i in 1:30){
    set.seed(i)
    coldCondition1 <- sample(coldCells,15)
    rt1Condition1 <- sample(rt1Cells,15)
    rt2Condition1 <- sample(rt2Cells,15)
    coldCondition2 <- sample(coldCells[!coldCells%in%coldCondition1],15)
    rt1Condition2 <- sample(rt1Cells[!rt1Cells%in%rt1Condition1],15)
    rt2Condition2 <- sample(rt2Cells[!rt2Cells%in%rt2Condition1],15)
  subsetMatrixUsoskinRandomCellType[i,] <- c(coldCondition1,rt1Condition1,rt2Condition1,
                   coldCondition2,rt1Condition2,rt2Condition2)
}
write.table(subsetMatrixUsoskinRandomCellType,file="~/Dropbox/PhD/Research/singleCell/usoskin/FPR/subsetMatrixUsoskinFPR_randomCellTypes.txt",row.names=FALSE,col.names = FALSE)


