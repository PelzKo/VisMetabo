library(subspace)
library(kohonen)
source("testdata.R")
#source("doc.R")

#setwd("C:\\Users\\Konstantin\\Desktop\\Uni\\6Semester\\BachelorArbeit\\paramTuning")

N <- 1000  # total number of rows to preallocate--possibly an overestimate

result <- data.frame(name=rep(NA, N), dimNames=rep("", N),  # as many cols as you need
                 diab1=rep(NA, N),diab2=rep(NA, N),diab3=rep(NA, N),diab4=rep(NA, N),
                 diab5=rep(NA, N),numberOfIndividuals=rep("", N),meanDist=rep(NA, N)) 
logs <- list()

small <- getData()
smallMetab <- small[[1]]
smallPheno <- small[[2]]
T2D <- smallPheno["T2D"]
big <- getData(FALSE)
bigMetab <- big[[1]]
bigPheno <- big[[2]]
genderTable <- bigPheno["GENDER"]
gender <- match(genderTable$GENDER,c("male","female"))-1 #0=male,1=female

counter <- 1
methods <- c("cosa")#c("clique","som","cosa","doc")
numberOfTimes <- 1:3

for (counterSmallBig in c("small","big")){
  if (counterSmallBig=="small"){
    data <- smallMetab
  } else {
    data <- bigMetab
  }
  
  for (method in methods){
    switch(method, 
           som={
             param1 <- 25#c(25,49,100)#25
             param2 <- 50#c(50,100,200)#50
             param3 <- "0.05,0.01"#c("0.05,0.01", "0.08,0.01", "0.05,0.001")#"0.05,0.01"
             
             combinations <- expand.grid(param1,param2,param3)
             for (row in 1:nrow(combinations)){
               for (round in numberOfTimes){
                 gridSize <- combinations[row,][[1]]
                 rlen <- combinations[row,][[2]]
                 alphaSplit <- strsplit(as.character(combinations[row, ][[3]]),",")[[1]]
                 alpha <- c(as.numeric(alphaSplit[[1]]),as.numeric(alphaSplit[[2]]))
      
                 print(paste("Starting",method,counterSmallBig,gridSize,rlen,alpha[[1]],alpha[[2]],sep = "_"))
      
                 run <- tryCatch({
                   somClust <- som(data.matrix(data), grid = somgrid(xdim = sqrt(gridSize),ydim = sqrt(gridSize), topo = "hexagonal"),rlen = rlen, alpha = alpha)
                   somClustTwoLayers <- xyf(data.matrix(data),data.matrix(data), grid = somgrid(xdim = sqrt(gridSize),ydim = sqrt(gridSize), topo = "hexagonal"),rlen = rlen, alpha = alpha)
                   
                   clusterCounts <- sort(table(somClust$unit.classif),decreasing = TRUE)
                   clusterCountsTwoLayers <- sort(table(somClustTwoLayers$unit.classif),decreasing = TRUE)
                   
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
                     if (counterSmallBig=="small"){
                        top5Diab[minicount] <- mean(T2D[[1]][currentCluster])
                     } else {
                        top5Diab[minicount] <- mean(gender[currentCluster])
                     }
                     minicount<-minicount+1
                   }
                   top5Diab[is.na(top5Diab)]<-0.5
                   minicount<-1
                   for (i in names(clusterCountsTwoLayers)[1:5]){
                     currentCluster <- ids[somClustTwoLayers$unit.classif==as.numeric(i)]
                     if (counterSmallBig=="small"){
                       top5DiabTwoLayers[minicount] <- mean(T2D[[1]][currentCluster])
                     } else {
                       top5DiabTwoLayers[minicount] <- mean(gender[currentCluster])
                     }
                     minicount<-minicount+1
                   }
                   top5DiabTwoLayers[is.na(top5DiabTwoLayers)]<-0.5
                   
                   result[counter, ] <- list(paste(method,counterSmallBig,"1Layer",gridSize,rlen,alpha[[1]],alpha[[2]],sep = "_"),"NO INFO" ,top5Diab[[1]],top5Diab[[2]],top5Diab[[3]],top5Diab[[4]],top5Diab[[5]],paste0(nums,collapse = "_"),mean(abs(0.5-top5Diab)))
                   counter <- counter + 1
                   result[counter, ] <- list(paste(method,counterSmallBig,"2Layer",gridSize,rlen,alpha[[1]],alpha[[2]],sep = "_"),"NO INFO" ,top5DiabTwoLayers[[1]],top5DiabTwoLayers[[2]],top5DiabTwoLayers[[3]],top5DiabTwoLayers[[4]],top5DiabTwoLayers[[5]],paste0(numsTwoLayers,collapse = "_"),mean(abs(0.5-top5DiabTwoLayers)))
                   counter <- counter + 1
                   return(NA)
                 },
                 error = function(cond){
                   log <- c(paste(method,counterSmallBig,gridSize,rlen,alpha[[1]],alpha[[2]],sep = "_"),cond$message,toString(cond$call))
                   return(log)
                 })
                 if (run[[2]]!="no function to return from, jumping to top level"){
                   logs <- c(logs,run)
                 }
                 print(paste(counter-2,counter-1,method,counterSmallBig,gridSize,rlen,alpha[[1]],alpha[[2]],sep = "_"))
               }
             }
           },
           cosa={
             param1 <- 0.2#c(0.1,0.2,0.4)
             param2 <- 1#c(2,4,7)
             param3 <- 1#c(30,80,150)
             
             combinations <- expand.grid(param1,param2,param3)
             for (row in 1:nrow(combinations)){
               for (round in 1){
                 lambda <- combinations[row,][[1]]
                 niter <- combinations[row,][[2]]
                 noit <- combinations[row,][[3]]
                 
                 run <- tryCatch({
                   cosaClust <- cosa2(data, lambda=lambda, niter = niter, noit = noit)
                   hist <- hierclust(cosaClust$D, denplot = FALSE)
                   clusters <- lapply(cut(hist$dendro, h=1.25)$lower,labels)
                   clusters <- clusters[order(sapply(clusters,length),decreasing=T)]
                   
                   nums <- sapply(clusters,length)[1:5]
                   nums[is.na(nums)]<-0
                   
                   top5Diab <- rep(0.5,5)
                   
                   for (i in 1:5){
                     if (i<=length(clusters)){
                       currentCluster <- clusters[[i]]
                       if (counterSmallBig=="small"){
                         top5Diab[i] <- mean(T2D[[1]][currentCluster])
                       } else {
                         top5Diab[i] <- mean(gender[currentCluster])
                       }
                     }
                   }
                   top5Diab[is.na(top5Diab)]<-0.5
                   
                   result[counter, ] <- list(paste(method,counterSmallBig,lambda,niter,noit,sep = "_"),"NO INFO" ,top5Diab[[1]],top5Diab[[2]],top5Diab[[3]],top5Diab[[4]],top5Diab[[5]],paste0(nums,collapse = "_"),mean(abs(0.5-top5Diab)))
                   counter <- counter + 1
                   return(NA)
                   
                   
                 },
                 error = function(cond){
                   log <- c(paste(method,counterSmallBig,lambda,niter,noit,sep = "_"),cond$message,toString(cond$call))
                   return(log)
                 })
                 if (run[[2]]!="keine Funktion abzubrechen, springe zum Top Level"){
                   logs <- c(logs,run)
                 }
               }
             }
            
             
             
           },
           doc={
             #param1 <- 0.6#c(0.1,0.2,0.4)
             #param2 <- 0.25#c(0.25,0.5,0.75)
             #
             #combinations <- expand.grid(param1,param2)
             #for (row in 1:nrow(combinations)){
             #  alpha <- 0.1#combinations[row,][[1]]
             #  beta <- 0.25#combinations[row,][[2]]
             #  w <- 15
             #  
             #  run <- tryCatch({
             #    test <- testDataFrame()
             #    docClust <- runDoc(test, alpha, beta, w)
             #    idsInClustersDoc <- getIdsDoc(docClust)
             #    dims <- getDimsDoc(docClust)
             #    
             #    #-------#
             #    cliqueClust <- cliqueClust[order(unlist(lapply(lapply(cliqueClust, "[[", "subspace"),"sum")),decreasing = TRUE)]
             #    nums <- unlist(lapply(lapply(cliqueClust, "[[", "subspace"),"sum"))
             #    top5 <- nums[1:5]
             #    top5Metabs <- list()
             #    top5Diab <- rep(0.5,5)
             #    
             #    
             #    for (i in 1:5){
             #      currentCluster <- idsInClustersDoc[[as.numeric(i)]]
             #      usedDimensions <- dims[as.numeric(i),]
             #      top5Metabs[[i]] <- paste0(names(data)[usedDimensions], collapse = "_")
             #      top5Diab[i] <- mean(T2D[[1]][currentCluster])
             #    }
             #    
             #    #realDiabs <- abs(0.5-top5Diab)
             #    
             #    result[counter, ] <- list(paste(method,xi,tau,sep = "_"),paste0(top5Metabs,collapse = ";") ,top5Diab[[1]],top5Diab[[2]],top5Diab[[3]],top5Diab[[4]],top5Diab[[5]])
             #    counter <- counter + 1
             #    return(NA)
             #  },
             #  error = function(cond){
             #    log <- c(paste(method,xi,tau,sep = "_"),cond$message,toString(cond$call))
             #    return(log)
             #  })
             #  if (run[[2]]!="keine Funktion abzubrechen, springe zum Top Level"){
             #    logs <- c(logs,run)
             #  }
             #}
           },
           {
             # default is using Clique
             param1 <- 13#c(2,5,8,11,14)#13
             param2 <- 0.2#c(0.04,0.08,0.1,0.3)#0.2
             
             combinations <- expand.grid(param1,param2)
             for (row in 1:nrow(combinations)){
               for (round in 1){
                 xi <- combinations[row,][[1]]
                 tau <- combinations[row,][[2]]
                 
                 print(paste("Starting",method,counterSmallBig,xi,tau,sep = "_"))
                 
                 run <- tryCatch({
                   cliqueClust <- CLIQUE(data, xi = xi,tau = tau)
                   cliqueClust <- cliqueClust[order(unlist(lapply(lapply(cliqueClust, "[[", "subspace"),"sum")),decreasing = TRUE)]
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
                     
                     
                     if (counterSmallBig=="small"){
                       top5Diab[i] <- mean(T2D[[1]][currentCluster$objects])
                     } else {
                       top5Diab[i] <- mean(gender[currentCluster$objects])
                     }
                   }
                   top5Diab[is.na(top5Diab)]<-0.5
                   
                   #realDiabs <- abs(0.5-top5Diab)
                   
                   result[counter, ] <- list(paste(method,counterSmallBig,xi,tau,sep = "_"),paste0(top5Metabs,collapse = ";") ,top5Diab[[1]],top5Diab[[2]],top5Diab[[3]],top5Diab[[4]],top5Diab[[5]],dimsAndIndividuals,mean(abs(0.5-top5Diab)))
                   counter <- counter + 1
                   return(NA)
                 },
                 error = function(cond){
                   log <- c(paste(method,counterSmallBig,xi,tau,sep = "_"),cond$message,toString(cond$call))
                   return(log)
                 })
                 if (run[[2]]!="no function to return from, jumping to top level"){
                   logs <- c(logs,run)
                 }
                 print(paste(counter-1,method,counterSmallBig,xi,tau,sep = "_"))
               }
             }
        }
    )
  }
}

if (length(logs)==0){
  logs<-c("No","Errors","found")
}

result <- result[complete.cases(result[[1]]),]
logsToWrite <- data.frame(matrix(unlist(logs), nrow=length(unlist(logs))/3, byrow=T))
write.table(result, file='output.tsv', quote = FALSE, row.names = FALSE, sep='\t')
write.table(logsToWrite, file='logs.tsv', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
