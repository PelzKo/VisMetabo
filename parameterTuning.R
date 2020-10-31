library(subspace)
library(kohonen)
library(RANN)
library(rCOSA)
source("testdata.R")
source("doc.R")

#setwd("C:\\Users\\Konstantin\\Desktop\\Uni\\6Semester\\BachelorArbeit\\paramTuning")

calcW <- function(data,c=1){
  neighbors <- nn2(data,k=2)
  pairs <- neighbors$nn.idx[,2]
  counter <- 1
  wLocals <- lapply(pairs, function(x){
    one <- data[counter,]
    counter <- counter + 1
    two <- data[x,]
    wLocal <- sum(abs(one-two))/length(one)
    wLocal
  })
  w <- (sum(unlist(wLocals))/length(wLocals))*c
  w
}

args <- commandArgs(trailingOnly = TRUE)
print(args)

method <- args[1]
counterSmallBig <- args[2]

N <- 1000  # total number of rows to preallocate--possibly an overestimate
numberOfTimes <- 1:3 #number of times to repeat algorithms which differ each iteration
counter <- 1 #current iteration

result <- data.frame(name=rep(NA, N), dimNames=rep(NA, N),  # as many cols as you need
                 diab1=rep(NA, N),diab2=rep(NA, N),diab3=rep(NA, N),diab4=rep(NA, N),
                 diab5=rep(NA, N),numberOfIndividuals=rep(NA, N),meanDist=rep(NA, N)) 
logs <- list()


if (counterSmallBig=="small"){
  small <- getData()
  data <- small[[1]]
  smallPheno <- small[[2]]
  pheno <- smallPheno[["T2D"]]
} else {
  big <- getData(FALSE)
  data <- big[[1]]
  bigPheno <- big[[2]]
  genderTable <- bigPheno["GENDER"]
  pheno <- match(genderTable$GENDER,c("male","female"))-1 #0=male,1=female
}
mean <- mean(pheno)
name <- ""
ptm <- proc.time()

args <- list("","","0.2","0.5","1")
switch(method, 
       som={
         run <- tryCatch({
           gridSize <- as.numeric(args[3])
           rlen <- as.numeric(args[4])
           alphaSplit <- strsplit(as.character(args[5]),",")[[1]]
           alpha <- c(as.numeric(alphaSplit[[1]]),as.numeric(alphaSplit[[2]]))
           alphaAsString <- paste(alpha[[1]],alpha[[2]],sep = ";")
           
           name <- paste(method,counterSmallBig,gridSize,rlen,alphaAsString,sep = "_")

           somClust <- som(data.matrix(data), grid = somgrid(xdim = sqrt(gridSize),ydim = sqrt(gridSize), topo = "hexagonal"),rlen = rlen, alpha = alpha)
           somClustTwoLayers <- xyf(data.matrix(data),data.matrix(data), grid = somgrid(xdim = sqrt(gridSize),ydim = sqrt(gridSize), topo = "hexagonal"),rlen = rlen, alpha = alpha)
           
           clusterCounts <- sort(table(somClust$unit.classif),decreasing = TRUE)
           clusterCounts <- clusterCounts[clusterCounts>4]
           clusterCountsTwoLayers <- sort(table(somClustTwoLayers$unit.classif),decreasing = TRUE)
           clusterCountsTwoLayers <- clusterCountsTwoLayers[clusterCountsTwoLayers>4]
           
           nums <- clusterCounts[1:5]
           nums[is.na(nums)]<-0
           numsTwoLayers <- clusterCountsTwoLayers[1:5]
           numsTwoLayers[is.na(numsTwoLayers)]<-0
           
           top5 <- nums[1:5]
           top5TwoLayers <- numsTwoLayers[1:5]
           
           top5Metabs <- list()
           top5MetabsTwoLayers <- list()
           
           top5Diab <- rep(0.5,5)
           top5DiabTwoLayers <- rep(0.5,5)
           
           ids <- 1:nrow(data)
           minicount<-1
           for (i in names(clusterCounts)[1:5]){
             currentCluster <- ids[somClust$unit.classif==as.numeric(i)]
             top5Diab[minicount] <- mean(pheno[currentCluster])
             minicount<-minicount+1
           }
           top5Diab[is.na(top5Diab)]<-0.5
           minicount<-1
           for (i in names(clusterCountsTwoLayers)[1:5]){
             currentCluster <- ids[somClustTwoLayers$unit.classif==as.numeric(i)]
             top5DiabTwoLayers[minicount] <- mean(pheno[currentCluster])
             minicount<-minicount+1
           }
           top5DiabTwoLayers[is.na(top5DiabTwoLayers)]<-0.5
           
           result[counter, ] <- list(paste(method,counterSmallBig,"1Layer",gridSize,rlen,alphaAsString,sep = "_"),"NO INFO" ,top5Diab[[1]],top5Diab[[2]],top5Diab[[3]],top5Diab[[4]],top5Diab[[5]],paste0(nums,collapse = "_"),mean(abs(mean-top5Diab)))
           counter <- counter + 1
           result[counter, ] <- list(paste(method,counterSmallBig,"2Layer",gridSize,rlen,alphaAsString,sep = "_"),"NO INFO" ,top5DiabTwoLayers[[1]],top5DiabTwoLayers[[2]],top5DiabTwoLayers[[3]],top5DiabTwoLayers[[4]],top5DiabTwoLayers[[5]],paste0(numsTwoLayers,collapse = "_"),mean(abs(mean-top5DiabTwoLayers)))
           counter <- counter + 1
           return(NA)
         },
         error = function(cond){
           log <- c(paste(method,counterSmallBig,gridSize,rlen,alphaAsString,sep = "_"),cond$message,toString(cond$call))
           return(log)
         })
         if (run[[2]]!="keine Funktion abzubrechen, springe zum Top Level"&run[[2]]!="no function to return from, jumping to top level"){
           logs <- c(logs,run)
         }
       },
       cosa={
         run <- tryCatch({
           lambda <- as.numeric(args[3])
           niter <- as.numeric(args[4])
           noit <- as.numeric(args[5])
           
           name <- paste(method,counterSmallBig,lambda,niter,noit,sep = "_")
         
           cosaClust <- cosa2(data, lambda=lambda, niter = niter, noit = noit)
           hist <- hierclust(cosaClust$D, denplot = FALSE)
           clusters <- lapply(cut(hist$dendro, h=1.25)$lower,labels)
           clusters <- clusters[order(sapply(clusters,length),decreasing=T)]
           clusters <- clusters[sapply(clusters, function(x) length(x)>4)]
           
           nums <- sapply(clusters,length)[1:5]
           nums[is.na(nums)]<-0
           
           top5Diab <- rep(0.5,5)
           
           for (i in 1:5){
             if (i<=length(clusters)){
               currentCluster <- clusters[[i]]
               top5Diab[i] <- mean(pheno[currentCluster])
             }
           }
           top5Diab[is.na(top5Diab)]<-0.5
           
           result[counter, ] <- list(paste(method,counterSmallBig,lambda,niter,noit,sep = "_"),"NO INFO" ,top5Diab[[1]],top5Diab[[2]],top5Diab[[3]],top5Diab[[4]],top5Diab[[5]],paste0(nums,collapse = "_"),mean(abs(mean-top5Diab)))
           counter <- counter + 1
           return(NA)
           
           
         },
         error = function(cond){
           log <- c(paste(method,counterSmallBig,lambda,niter,noit,sep = "_"),cond$message,toString(cond$call))
           return(log)
         })
         if (run[[2]]!="keine Funktion abzubrechen, springe zum Top Level"&run[[2]]!="no function to return from, jumping to top level"){
           logs <- c(logs,run)
         }
       },
       doc={
         for (round in numberOfTimes){
             run <- tryCatch({
               alpha <- as.numeric(args[3])
               beta <- as.numeric(args[4])
               w <- calcW(data,as.numeric(args[5]))
               
               name <- paste(method,counterSmallBig,alpha,beta,w,sep = "_")
             
               docClust <- runDoc(data, alpha, beta, w)
               dimNums <- getDimsDoc(docClust)
               #dims <- lapply(dimNums,function(x) names(data)[x])
               dimNums <- dimNums[order(unlist(lapply(dimNums,"length")),decreasing = TRUE)]
               
               dims <- unlist(lapply(dimNums,"length"))[1:5]
               dims[is.na(dims)]<-0
               
               idsInClustersDoc <- getIdsDoc(docClust)
               idsInClustersDoc <- idsInClustersDoc[order(unlist(lapply(dimNums,"length")),decreasing = TRUE)]
               
               numInd <- unlist(lapply(idsInClustersDoc,length))[1:5]
               numInd[is.na(numInd)]<-0
               
               dimCondition <- sapply(dimNums, function(x) return(sum(x)>1))
               idsCondition <- sapply(idsInClustersDoc, function(x) return(length(x)>4))
               condition <- dimCondition&idsCondition
               
               dimNums <- dimNums[condition]
               idsInClustersDoc <- idsInClustersDoc[condition]
               
               dimsAndIndividuals <- paste0(c(dims,numInd),collapse="_")
               
               top5Metabs <- list()
               top5Diab <- rep(0.5,5)
               
               
               for (i in 1:5){
                 if (i<=length(idsInClustersDoc)&i<=length(dimNums)){
                   currentCluster <- idsInClustersDoc[[as.numeric(i)]]
                   top5Metabs[[i]] <- paste0(names(data)[dimNums[[as.numeric(i)]]], collapse = "_")
                   
                   top5Diab[i] <- mean(pheno[currentCluster])
                 }
               }
               top5Diab[is.na(top5Diab)]<-0.5
               
               #realDiabs <- abs(0.5-top5Diab)
               
               result[counter, ] <- list(paste(method,counterSmallBig,alpha,beta,w,sep = "_"),paste0(top5Metabs,collapse = ";") ,top5Diab[[1]],top5Diab[[2]],top5Diab[[3]],top5Diab[[4]],top5Diab[[5]],dimsAndIndividuals,mean(abs(mean-top5Diab)))
               counter <- counter + 1
               return(NA)
             },
             error = function(cond){
               log <- c(paste(method,counterSmallBig,alpha,beta,w,sep = "_"),cond$message,toString(cond$call))
               return(log)
             })
             if (run[[2]]!="keine Funktion abzubrechen, springe zum Top Level"&run[[2]]!="no function to return from, jumping to top level"){
               logs <- c(logs,run)
             }
          }
       },
       {
         # default is using Clique
         run <- tryCatch({
           xi <- as.numeric(args[3])
           tau <- as.numeric(args[4])
           
           name <- paste(method,counterSmallBig,xi,tau,sep = "_")
         
           cliqueClust <- CLIQUE(data, xi = xi,tau = tau)
           cliqueClust <- cliqueClust[order(unlist(lapply(lapply(cliqueClust, "[[", "subspace"),"sum")),decreasing = TRUE)]
           cliqueClust <- cliqueClust[sapply(cliqueClust, function(x) return(sum(x$subspace)>1&length(x$objects)>4))]
           
           
           dims <- unlist(lapply(lapply(cliqueClust, "[[", "subspace"),"sum"))[1:5]
           dims[is.na(dims)]<-0
           numInd <- unlist(lapply(lapply(cliqueClust, "[[", "objects"),"length"))[1:5]
           numInd[is.na(numInd)]<-0
           
           dimsAndIndividuals <- paste0(c(dims,numInd),collapse="_")

           top5Metabs <- list()
           top5Diab <- rep(0.5,5)
           
           
           for (i in 1:5){
             currentCluster <- cliqueClust[as.numeric(i)][[1]]
             top5Metabs[[i]] <- paste0(names(data)[currentCluster$subspace], collapse = "_")
             
             mean(pheno[currentCluster$objects])
           }
           top5Diab[is.na(top5Diab)]<-0.5
           
           #realDiabs <- abs(0.5-top5Diab)
           
           result[counter, ] <- list(paste(method,counterSmallBig,xi,tau,sep = "_"),paste0(top5Metabs,collapse = ";") ,top5Diab[[1]],top5Diab[[2]],top5Diab[[3]],top5Diab[[4]],top5Diab[[5]],dimsAndIndividuals,mean(abs(mean-top5Diab)))
           counter <- counter + 1
           return(NA)
         },
         error = function(cond){
           log <- c(paste(method,counterSmallBig,xi,tau,sep = "_"),cond$message,toString(cond$call))
           return(log)
         })
         if (run[[2]]!="keine Funktion abzubrechen, springe zum Top Level"&run[[2]]!="no function to return from, jumping to top level"){
           logs <- c(logs,run)
         }
    }
)

timeRes <- proc.time() - ptm



result <- result[complete.cases(result[[1]]),]
result$time <- timeRes[[3]]
write.table(result, file=paste0(name,'_output.tsv'), quote = FALSE, row.names = FALSE, sep='\t')
if (length(logs)!=0){
  logsToWrite <- data.frame(matrix(unlist(logs), nrow=length(unlist(logs))/3, byrow=T))
  write.table(logsToWrite, file=paste0(name,'_logs.tsv'), quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
}
