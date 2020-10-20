library(subspace)
library(kohonen)
source("C:\\Users\\Konstantin\\Desktop\\Uni\\6Semester\\BachelorArbeit\\code\\testdata.R")
source("C:\\Users\\Konstantin\\Desktop\\Uni\\6Semester\\BachelorArbeit\\code\\doc.R")

N <- 1000  # total number of rows to preallocate--possibly an overestimate

result <- data.frame(name=rep(NA, N), dimNames=rep("", N),  # as many cols as you need
                 diab1=rep(NA, N),diab2=rep(NA, N),diab3=rep(NA, N),diab4=rep(NA, N),
                 diab5=rep(NA, N),numberOfIndividuals=rep("", N),meanDist=rep(NA, N)) 
logs <- list()

small <- getData()
smallMetab <- small[[1]]
smallPheno <- small[[2]]
T2D <- smallPheno["T2D"]
#big <- getData(FALSE)
#bigMetab <- big[[1]]
#bigPheno <- big[[2]]

counter<-1
methods <- c("clique","som")#c("clique","som","cosa","doc")
numberOfTimes <- 1:3

for (method in methods){
  switch(method, 
         som={
           param1 <- 25#c(25,49,100)
           param2 <- 50#c(50,100,200)
           param3 <- "0.05,0.01"#c("0.05,0.01", "0.08,0.01", "0.05,0.001")
           
           combinations <- expand.grid(param1,param2,param3)
           for (row in 1:nrow(combinations)){
             for (round in numberOfTimes){
               gridSize <- combinations[row,][[1]]
               rlen <- combinations[row,][[2]]
               alphaSplit <- strsplit(as.character(combinations[row, ][[3]]),",")[[1]]
               alpha <- c(as.numeric(alphaSplit[[1]]),as.numeric(alphaSplit[[2]]))
               
               run <- tryCatch({
                 somClust <- som(data.matrix(smallMetab), grid = somgrid(xdim = sqrt(gridSize),ydim = sqrt(gridSize), topo = "hexagonal"),rlen = rlen, alpha = alpha)
                 somClustTwoLayers <- xyf(data.matrix(smallMetab),data.matrix(smallMetab), grid = somgrid(xdim = sqrt(gridSize),ydim = sqrt(gridSize), topo = "hexagonal"),rlen = rlen, alpha = alpha)
                 
                 clusterCounts <- sort(table(somClust$unit.classif),decreasing = TRUE)
                 clusterCountsTwoLayers <- sort(table(somClustTwoLayers$unit.classif),decreasing = TRUE)
                 
                 nums <- clusterCounts[1:5]
                 numsTwoLayers <- clusterCountsTwoLayers[1:5]
                 
                 top5 <- nums[1:5]
                 top5TwoLayers <- numsTwoLayers[1:5]
                 
                 top5Metabs <- list()
                 top5MetabsTwoLayers <- list()
                 
                 top5Diab <- rep(0.5,5)
                 top5DiabTwoLayers <- rep(0.5,5)
                 
                 ids <- 1:nrow(smallMetab)
                 minicount<-1
                 for (i in names(clusterCounts)[1:5]){
                   currentCluster <- ids[somClust$unit.classif==as.numeric(i)]
                   top5Diab[minicount] <- mean(T2D[[1]][currentCluster])
                   minicount<-minicount+1
                 }
                 minicount<-1
                 for (i in names(clusterCountsTwoLayers)[1:5]){
                   currentCluster <- ids[somClustTwoLayers$unit.classif==as.numeric(i)]
                   top5DiabTwoLayers[minicount] <- mean(T2D[[1]][currentCluster])
                   minicount<-minicount+1
                 }
                 
                 result[counter, ] <- list(paste(method,"1Layer",gridSize,rlen,alpha[[1]],alpha[[2]],sep = "_"),"NO INFO" ,top5Diab[[1]],top5Diab[[2]],top5Diab[[3]],top5Diab[[4]],top5Diab[[5]],paste0(nums,collapse = "_"),mean(abs(0.5-top5Diab)))
                 counter <- counter + 1
                 result[counter, ] <- list(paste(method,"2Layer",gridSize,rlen,alpha[[1]],alpha[[2]],sep = "_"),"NO INFO" ,top5DiabTwoLayers[[1]],top5DiabTwoLayers[[2]],top5DiabTwoLayers[[3]],top5DiabTwoLayers[[4]],top5DiabTwoLayers[[5]],paste0(numsTwoLayers,collapse = "_"),mean(abs(0.5-top5DiabTwoLayers)))
                 counter <- counter + 1
                 return(NA)
               },
               error = function(cond){
                 log <- c(paste(method,gridSize,rlen,alpha,sep = "_"),cond$message,toString(cond$call))
                 return(log)
               })
               if (run[[2]]!="keine Funktion abzubrechen, springe zum Top Level"){
                 logs <- c(logs,run)
               }
             }
           }
         },
         cosa={
           #param1 <- 13#c(2,5,8,11,14)
           #param2 <- c(0.2,0.9,0.8)#c(0.04,0.08,0.1,0.3)
           #
           #combinations <- expand.grid(param1,param2)
           #for (row in 1:nrow(combinations)){
           #  xi <- combinations[row,][[1]]
           #  tau <- combinations[row,][[2]]
           #  
           #  run <- tryCatch({
           #    cliqueClust <- CLIQUE(smallMetab, xi = xi,tau = tau)
           #    cliqueClust <- cliqueClust[order(unlist(lapply(lapply(cliqueClust, "[[", "subspace"),"sum")),decreasing = TRUE)]
           #    nums <- unlist(lapply(lapply(cliqueClust, "[[", "subspace"),"sum"))
           #    top5 <- nums[1:5]
           #    top5Metabs <- list()
           #    top5Diab <- rep(0.5,5)
           #    
           #    
           #    for (i in 1:5){
           #      currentCluster <- cliqueClust[as.numeric(i)][[1]]
           #      top5Metabs[[i]] <- paste0(names(smallMetab)[currentCluster$subspace], collapse = "_")
           #      top5Diab[i] <- mean(T2D[[1]][currentCluster$objects])
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
         doc={
           param1 <- 0.6#c(0.1,0.2,0.4)
           param2 <- 0.25#c(0.25,0.5,0.75)
           
           combinations <- expand.grid(param1,param2)
           for (row in 1:nrow(combinations)){
             alpha <- 0.1#combinations[row,][[1]]
             beta <- 0.25#combinations[row,][[2]]
             w <- 15
             
             run <- tryCatch({
               test <- testDataFrame()
               docClust <- runDoc(test, alpha, beta, w)
               idsInClustersDoc <- getIdsDoc(docClust)
               dims <- getDimsDoc(docClust)
               
               #-------#
               cliqueClust <- cliqueClust[order(unlist(lapply(lapply(cliqueClust, "[[", "subspace"),"sum")),decreasing = TRUE)]
               nums <- unlist(lapply(lapply(cliqueClust, "[[", "subspace"),"sum"))
               top5 <- nums[1:5]
               top5Metabs <- list()
               top5Diab <- rep(0.5,5)
               
               
               for (i in 1:5){
                 currentCluster <- idsInClustersDoc[[as.numeric(i)]]
                 usedDimensions <- dims[as.numeric(i),]
                 top5Metabs[[i]] <- paste0(names(smallMetab)[usedDimensions], collapse = "_")
                 top5Diab[i] <- mean(T2D[[1]][currentCluster])
               }
               
               #realDiabs <- abs(0.5-top5Diab)
               
               result[counter, ] <- list(paste(method,xi,tau,sep = "_"),paste0(top5Metabs,collapse = ";") ,top5Diab[[1]],top5Diab[[2]],top5Diab[[3]],top5Diab[[4]],top5Diab[[5]])
               counter <- counter + 1
               return(NA)
             },
             error = function(cond){
               log <- c(paste(method,xi,tau,sep = "_"),cond$message,toString(cond$call))
               return(log)
             })
             if (run[[2]]!="keine Funktion abzubrechen, springe zum Top Level"){
               logs <- c(logs,run)
             }
           }
         },
         {
           # default is using Clique
           param1 <- 13#c(2,5,8,11,14)
           param2 <- 0.2#c(0.04,0.08,0.1,0.3)
           
           combinations <- expand.grid(param1,param2)
           for (row in 1:nrow(combinations)){
             for (round in numberOfTimes){
               xi <- combinations[row,][[1]]
               tau <- combinations[row,][[2]]
               
               run <- tryCatch({
                 cliqueClust <- CLIQUE(smallMetab, xi = xi,tau = tau)
                 cliqueClust <- cliqueClust[order(unlist(lapply(lapply(cliqueClust, "[[", "subspace"),"sum")),decreasing = TRUE)]
                 nums <- unlist(lapply(lapply(cliqueClust, "[[", "subspace"),"sum"))[1:5]
                 top5Metabs <- list()
                 top5Diab <- rep(0.5,5)
                 
                 
                 for (i in 1:5){
                   currentCluster <- cliqueClust[as.numeric(i)][[1]]
                   top5Metabs[[i]] <- paste0(names(smallMetab)[currentCluster$subspace], collapse = "_")
                   top5Diab[i] <- mean(T2D[[1]][currentCluster$objects])
                 }
                 
                 #realDiabs <- abs(0.5-top5Diab)
                 
                 result[counter, ] <- list(paste(method,xi,tau,sep = "_"),paste0(top5Metabs,collapse = ";") ,top5Diab[[1]],top5Diab[[2]],top5Diab[[3]],top5Diab[[4]],top5Diab[[5]],paste0(nums,collapse = "_"),mean(abs(0.5-top5Diab)))
                 counter <- counter + 1
                 return(NA)
               },
               error = function(cond){
                 log <- c(paste(method,xi,tau,sep = "_"),cond$message,toString(cond$call))
                 return(log)
               })
               if (run[[2]]!="keine Funktion abzubrechen, springe zum Top Level"){
                 logs <- c(logs,run)
               }
             }
           }
      }
  )
}


result <- result[complete.cases(result[[1]]),]
logsToWrite <- data.frame(matrix(unlist(logs), nrow=length(unlist(logs))/3, byrow=T))
write.table(result, file='output.tsv', quote = FALSE, row.names = FALSE, sep='\t')
write.table(logsToWrite, file='logs.tsv', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
