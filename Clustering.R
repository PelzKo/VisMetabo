library(subspace)
library(shiny)
library(ape)

source("C:/Users/Konstantin/Desktop/Uni/6Semester/BachelorArbeit/Code/ReadingData.R")
options(java.parameters = "-Xmx16000m")

metabs <- readFile(dftest)
metabTest <- metabs[[1]]

metabTest <- metabTest[complete.cases(metabTest), ]
#zVar <- (myVar - mean(myVar)) / sd(myVar)
metabTest <- scale(metabTest)

clusterInfo <- SubClu(metabTest)
length(clusterInfo)

clusters <- rep("grey",length(metabTest[,1]))
colors <- c("red","blue","green","purple")
for (i in seq_len(length(clusterInfo))){
  metabTest[clusterInfo[[i]][[2]],clusterInfo[[i]][[1]]]<-10#metabTest[clusterInfo[[i]][[2]],clusterInfo[[i]][[1]]]*100
  #clusters[clusterInfo[[i]][[2]]]=colors[[i]]
}
#clusters[metabs[[2]]$GENDER=="male"]<-rgb(20, 0, 0, maxColorValue=255, alpha=255)
#clusters[metabs[[2]]$GENDER=="female"]<-rgb(100, 0, 0, maxColorValue=255, alpha=255)

#phenotype <- as.integer(metabs[[2]]$GENDER=="male")
phenotype <- metabs[[2]]$AGE
clusters <- rgb(range01(phenotype)*255, 0, 0, maxColorValue=255, alpha=255)


test <- dist(metabTest)
res <- pcoa(test)

plot(res$vectors,col=clusters,bg = clusters,pch = 21)


range01 <- function(x){(x-min(x))/(max(x)-min(x))}
