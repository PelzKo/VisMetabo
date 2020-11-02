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

getExpressionTableForCluster <- function(cluster,pheno){
  phenoClust <- length(cluster[pheno[cluster]==1])
  phenoNonClust <- length(pheno[pheno==1])-phenoClust
  nonPhenoClust <- length(cluster)-phenoClust
  nonPhenoNonClust <- length(pheno[pheno==0])-nonPhenoClust
  if (phenoClust+phenoNonClust+nonPhenoClust+nonPhenoNonClust!=length(pheno)){
    print("ERROR, table is not correct")
    return(NULL)
  }
  result <- data.frame(matrix(c(phenoClust,phenoNonClust,nonPhenoClust,nonPhenoNonClust),nrow = 2,ncol = 2))
  names(result) <- c("pheno","nonPheno")
  row.names(result) <- c("inCluster","outCluster")
  result
}

calcValues <- function(idList,dimList,time){
  significantClusters = numeric()
  insignificantClusters = numeric()
  pValues = numeric()
  enrichmentScore = numeric()
  
  phenoAvg <- numeric()
  phenoAvgSig <- numeric()
  phenoAvgInSig <- numeric()
  
  for (i in seq_len(length(idList))){
    expressionTableCluster <- getExpressionTableForCluster(idList[[i]],pheno)
    chi <- chisq.test(expressionTableCluster)
    pValues <- c(pValues,chi$p.value)
    enrichmentScore <- c(enrichmentScore,chi$statistic)
    phenoAvg <- c(phenoAvg,mean(pheno[idList[[i]]]))
    if (chi$p.value<0.05){
      significantClusters <- c(significantClusters,i)
      phenoAvgSig <- c(phenoAvgSig,mean(pheno[idList[[i]]]))
    } else {
      insignificantClusters <- c(insignificantClusters,i)
      phenoAvgInSig <- c(phenoAvgInSig,mean(pheno[idList[[i]]]))
    }
  }
  
  dims <- sapply(dimList,sum)
  numInd <- sapply(idList,length)
  
  dimsSig <- dims[significantClusters]
  numIndSig <- numInd[significantClusters]
  
  dimsInSig <- dims[-significantClusters]
  numIndInSig <- numInd[-significantClusters]
  
  phenoAvg <- abs(mean-phenoAvg)
  phenoAvgSig <- abs(mean-phenoAvgSig)
  phenoAvgInSig <- abs(mean-phenoAvgInSig)
  
  idMaxEnrichment <- seq_len(length(enrichmentScore))[enrichmentScore == max(enrichmentScore)][[1]]
  
  maxEnrichment <- paste(enrichmentScore[[idMaxEnrichment]],pValues[[idMaxEnrichment]],sep = "_")
  numberClusters <- paste(length(significantClusters),length(idList)-length(significantClusters),length(idList),sep ="_")
  meanEnrichment <- paste(mean(enrichmentScore[significantClusters]),mean(enrichmentScore[-significantClusters]),mean(enrichmentScore),sep = "_")
  meanPValue <- paste(mean(pValues[significantClusters]),mean(pValues[-significantClusters]),mean(pValues),sep = "_")
  sumEnrichment <- paste(sum(enrichmentScore[significantClusters]),sum(enrichmentScore[-significantClusters]),sum(enrichmentScore),sep = "_")
  sumPValue <- paste(sum(pValues[significantClusters]),sum(pValues[-significantClusters]),sum(pValues),sep = "_")
  meanInd <- paste(mean(numIndSig),mean(numIndInSig),mean(numInd),sep = "_")
  meanDims <- paste(mean(dimsSig),mean(dimsInSig),mean(dims),sep = "_")
  meanPheno <- paste(mean(phenoAvgSig),mean(phenoAvgInSig),mean(phenoAvg),sep = "_")
  
  return(list(name,maxEnrichment,numberClusters,meanEnrichment,meanPValue,sumEnrichment,sumPValue,meanInd,meanDims,meanPheno,time))
}

args <- commandArgs(trailingOnly = TRUE)
#args <- list("som","big","49","50","0.08,0.01","1Layer")
print(args)

method <- args[1]
counterSmallBig <- args[2]

N <- 1000  # total number of rows to preallocate--possibly an overestimate
numberOfTimes <- 1:3 #number of times to repeat algorithms which differ each iteration
counter <- 1 #current iteration

result <- data.frame(name=rep(NA, N), maxEnrichment=rep(NA, N),  # as many cols as you need
                     numberClusters=rep(NA, N),meanEnrichment=rep(NA, N),meanPValue=rep(NA, N),sumEnrichment=rep(NA, N),
                     sumPValue=rep(NA, N),meanInd=rep(NA, N),meanDims=rep(NA, N),meanPheno=rep(NA, N),time=rep(NA, N)) 
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
switch(method, 
       som={
         for (round in numberOfTimes){
           ptm <- proc.time()
           run <- tryCatch({
             gridSize <- as.numeric(args[3])
             rlen <- as.numeric(args[4])
             alphaSplit <- strsplit(as.character(args[5]),",")[[1]]
             alpha <- c(as.numeric(alphaSplit[[1]]),as.numeric(alphaSplit[[2]]))
             alphaAsString <- paste(alpha[[1]],alpha[[2]],sep = ";")
             oneLayer <- args[6]=="1Layer"
             
             name <- paste(method,counterSmallBig,gridSize,rlen,alphaAsString,sep = "_")
             
             if (oneLayer){
               somClust <- som(data.matrix(data), grid = somgrid(xdim = sqrt(gridSize),ydim = sqrt(gridSize), topo = "hexagonal"),rlen = rlen, alpha = alpha)
             } else {
               somClust <- xyf(data.matrix(data),data.matrix(data), grid = somgrid(xdim = sqrt(gridSize),ydim = sqrt(gridSize), topo = "hexagonal"),rlen = rlen, alpha = alpha)
             }
             timeRes <- proc.time() - ptm
             
             ids <- 1:nrow(data)
             clusters <- list()
             clusterCounts <- sort(table(somClust$unit.classif),decreasing = TRUE)
             if (counterSmallBig=="small"){
               clusterCounts <- clusterCounts[clusterCounts>4]
             } else{
               clusterCounts <- clusterCounts[clusterCounts>9]
             }
             
             minicount<-1
             for (i in names(clusterCounts)){
               clusters[[minicount]] <- ids[somClust$unit.classif==as.numeric(i)]
               minicount<-minicount+1
             }
             
             idList <- clusters
             dimList <- list()
             
             if (length(idList)==0){
               result[counter, ] <- list(name,NA,NA,NA,NA,NA,NA,NA,NA,NA)
               counter <- counter + 1
               return(NA)
             }
             
             result[counter, ] <- calcValues(idList,dimList,timeRes[[3]])
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
         }
       },
       cosa={
         ptm <- proc.time()
         run <- tryCatch({
           lambda <- as.numeric(args[3])
           niter <- as.numeric(args[4])
           noit <- as.numeric(args[5])
           
           name <- paste(method,counterSmallBig,lambda,niter,noit,sep = "_")
         
           cosaClust <- cosa2(data, lambda=lambda, niter = niter, noit = noit)
           timeRes <- proc.time() - ptm
           
           hist <- hierclust(cosaClust$D, denplot = FALSE)
           clusters <- lapply(cut(hist$dendro, h=1.25)$lower,labels)
           clusters <- clusters[order(sapply(clusters,length),decreasing=T)]
           if (counterSmallBig=="small"){
             clusters <- clusters[sapply(clusters, function(x) length(x)>4)]
           } else{
             clusters <- clusters[sapply(clusters, function(x) length(x)>9)]
           }
           
           idList <- clusters
           dimList <- list()
           
           if (length(idList)==0){
             result[counter, ] <- list(name,NA,NA,NA,NA,NA,NA,NA,NA,NA)
             counter <- counter + 1
             return(NA)
           }
           
           result[counter, ] <- calcValues(idList,dimList,timeRes[[3]])
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
             ptm <- proc.time()
             run <- tryCatch({
               alpha <- as.numeric(args[3])
               beta <- as.numeric(args[4])
               w <- calcW(data,as.numeric(args[5]))
               
               name <- paste(method,counterSmallBig,alpha,beta,w,sep = "_")
             
               docClust <- runDoc(data, alpha, beta, w)
               timeRes <- proc.time() - ptm
               
               dimNums <- getDimsDoc(docClust)
               dimNums <- dimNums[order(unlist(lapply(dimNums,"length")),decreasing = TRUE)]
               
               idsInClustersDoc <- getIdsDoc(docClust)
               idsInClustersDoc <- idsInClustersDoc[order(unlist(lapply(dimNums,"length")),decreasing = TRUE)]
               
               dimCondition <- sapply(dimNums, function(x) return(sum(x)>1))
               idsCondition <- sapply(idsInClustersDoc, function(x) return(length(x)>4))
               if (counterSmallBig=="small"){
                 idsCondition <- sapply(idsInClustersDoc, function(x) return(length(x)>4))
               } else{
                 idsCondition <- sapply(idsInClustersDoc, function(x) return(length(x)>9))
               }
               condition <- dimCondition&idsCondition
               
               idList <- idsInClustersDoc[condition]
               dimList <- dimNums[condition]
               
               if (length(idList)==0){
                 result[counter, ] <- list(name,NA,NA,NA,NA,NA,NA,NA,NA,NA)
                 counter <- counter + 1
                 return(NA)
               }
               
               result[counter, ] <- calcValues(idList,dimList,timeRes[[3]])
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
         ptm <- proc.time()
         run <- tryCatch({
           xi <- as.numeric(args[3])
           tau <- as.numeric(args[4])
           
           name <- paste(method,counterSmallBig,xi,tau,sep = "_")
         
           cliqueClust <- CLIQUE(data, xi = xi,tau = tau)
           timeRes <- proc.time() - ptm
           
           cliqueClust <- cliqueClust[order(unlist(lapply(lapply(cliqueClust, "[[", "subspace"),"sum")),decreasing = TRUE)]
           if (counterSmallBig=="small"){
             cliqueClust <- cliqueClust[sapply(cliqueClust, function(x) return(sum(x$subspace)>1&length(x$objects)>4))]
           } else{
             cliqueClust <- cliqueClust[sapply(cliqueClust, function(x) return(sum(x$subspace)>1&length(x$objects)>9))]
           }
           
           idList <- lapply(cliqueClust, "[[", "objects")
           dimList <- lapply(cliqueClust, "[[", "subspace")
           
           if (length(idList)==0){
             result[counter, ] <- list(name,NA,NA,NA,NA,NA,NA,NA,NA,NA)
             counter <- counter + 1
             return(NA)
           }
           
           result[counter, ] <- calcValues(idList,dimList,timeRes[[3]])
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

result <- result[complete.cases(result[[1]]),]
write.table(result, file=paste0("results/",name,'_output.tsv'), quote = FALSE, row.names = FALSE, sep='\t')
if (length(logs)!=0){
  logsToWrite <- data.frame(matrix(unlist(logs), nrow=length(unlist(logs))/3, byrow=T))
  write.table(logsToWrite, file=paste0("results/",name,'_logs.tsv'), quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
}
