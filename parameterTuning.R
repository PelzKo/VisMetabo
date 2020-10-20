library(subspace)
library(kohonen)
source("C:\\Users\\Konstantin\\Desktop\\Uni\\6Semester\\BachelorArbeit\\code\\testdata.R")
source("C:\\Users\\Konstantin\\Desktop\\Uni\\6Semester\\BachelorArbeit\\code\\doc.R")

N <- 1000  # total number of rows to preallocate--possibly an overestimate

result <- data.frame(name=rep(NA, N), dimNames=rep("", N),  # as many cols as you need
                 diab1=rep(NA, N),diab2=rep(NA, N),diab3=rep(NA, N),diab4=rep(NA, N),
                 diab5=rep(NA, N)) 
logs <- list()

small <- getData()
smallMetab <- small[[1]]
smallPheno <- small[[2]]
T2D <- smallPheno["T2D"]
#big <- getData(FALSE)
#bigMetab <- big[[1]]
#bigPheno <- big[[2]]

counter<-1
currentCombinationCounter <- 1

for (method in c("clique","som","cosa","doc")){
  switch(method, 
         som={
           param1 <- 25#c(25,49,100)
           param2 <- 50#c(50,100,200)
           param3 <- c(0.05,0.01)#c(c(0.05,0.01), c(0.08,0.01), c(0.05,0.001)) <---GEHT NOCH NICHT!!
           
           combinations <- expand.grid(param1,param2,param3)
           for (row in 1:nrow(combinations)){
             gridSize <- combinations[row,][[1]]
             rlen <- combinations[row,][[2]]
             alpha <- combinations[row,][[3]]
             
             run <- tryCatch({
               somClust <- som(data.matrix(smallMetab), grid = somgrid(xdim = sqrt(gridSize),ydim = sqrt(gridSize), topo = "hexagonal"),rlen = rlen, alpha = alpha)
               somClustTwoLayers <- xyf(data.matrix(smallMetab),data.matrix(smallMetab), grid = somgrid(xdim = sqrt(gridSize),ydim = sqrt(gridSize), topo = "hexagonal"),rlen = rlen, alpha = alpha)
               
               cliqueClust <- cliqueClust[order(unlist(lapply(lapply(cliqueClust, "[[", "subspace"),"sum")),decreasing = TRUE)]
               nums <- unlist(lapply(lapply(cliqueClust, "[[", "subspace"),"sum"))
               top5 <- nums[1:5]
               top5Metabs <- list()
               top5Diab <- rep(0.5,5)
               
               
               clusterNumbers[clusteringData$SOM$unit.classif==input$clusterId] <- 17
               output$metabUsed <- renderUI(HTML(sprintf("Cluster %s selected",input$clusterId)))
               finalValues$idsFromCluster <- finalValues$numFromId[clusteringData$SOM$unit.classif==input$clusterId]
               
               
               for (i in 1:5){
                 currentCluster <- cliqueClust[as.numeric(i)][[1]]
                 top5Metabs[[i]] <- paste0(names(smallMetab)[currentCluster$subspace], collapse = "_")
                 top5Diab[i] <- mean(T2D[[1]][currentCluster$objects])
               }
               
               #realDiabs <- abs(0.5-top5Diab)
               
               result[counter, ] <- list(paste(method,xi,tau,sep = "_"),paste0(top5Metabs,collapse = ";") ,top5Diab[[1]],top5Diab[[2]],top5Diab[[3]],top5Diab[[4]],top5Diab[[5]])
               counter <- counter + 1
               return(NA)
             },
             error = function(cond){
               log <- c(paste(method,xi,tau,sep = "_"),cond$message,toString(cond$call))
               print(log)
               return(log)
             })
             if (run[[2]]!="keine Funktion abzubrechen, springe zum Top Level"){
               logs <- c(logs,run)
             }
             #currentCombinationCounter <- currentCombinationCounter + 1
           }
         },
         cosa={
           param1 <- 13#c(2,5,8,11,14)
           param2 <- c(0.2,0.9,0.8)#c(0.04,0.08,0.1,0.3)
           
           combinations <- expand.grid(param1,param2)
           for (row in 1:nrow(combinations)){
             xi <- combinations[row,][[1]]
             tau <- combinations[row,][[2]]
             
             print(xi)
             print(tau)
             run <- tryCatch({
               cliqueClust <- CLIQUE(smallMetab, xi = xi,tau = tau)
               cliqueClust <- cliqueClust[order(unlist(lapply(lapply(cliqueClust, "[[", "subspace"),"sum")),decreasing = TRUE)]
               nums <- unlist(lapply(lapply(cliqueClust, "[[", "subspace"),"sum"))
               top5 <- nums[1:5]
               top5Metabs <- list()
               top5Diab <- rep(0.5,5)
               
               
               for (i in 1:5){
                 currentCluster <- cliqueClust[as.numeric(i)][[1]]
                 top5Metabs[[i]] <- paste0(names(smallMetab)[currentCluster$subspace], collapse = "_")
                 top5Diab[i] <- mean(T2D[[1]][currentCluster$objects])
               }
               
               #realDiabs <- abs(0.5-top5Diab)
               
               result[counter, ] <- list(paste(method,xi,tau,sep = "_"),paste0(top5Metabs,collapse = ";") ,top5Diab[[1]],top5Diab[[2]],top5Diab[[3]],top5Diab[[4]],top5Diab[[5]])
               counter <- counter + 1
               return(NA)
             },
             error = function(cond){
               log <- c(paste(method,xi,tau,sep = "_"),cond$message,toString(cond$call))
               print(log)
               return(log)
             })
             if (run[[2]]!="keine Funktion abzubrechen, springe zum Top Level"){
               logs <- c(logs,run)
             }
             #currentCombinationCounter <- currentCombinationCounter + 1
           }
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
               docClust <- runDoc(smallMetab, alpha, beta, w)
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
               print(log)
               return(log)
             })
             if (run[[2]]!="keine Funktion abzubrechen, springe zum Top Level"){
               logs <- c(logs,run)
             }
             #currentCombinationCounter <- currentCombinationCounter + 1
           }
         },
         {
           # default is using Clique
           param1 <- 13#c(2,5,8,11,14)
           param2 <- 0.2#c(0.04,0.08,0.1,0.3)
           
           combinations <- expand.grid(param1,param2)
           for (row in 1:nrow(combinations)){
             xi <- combinations[row,][[1]]
             tau <- combinations[row,][[2]]
             
             run <- tryCatch({
               cliqueClust <- CLIQUE(smallMetab, xi = xi,tau = tau)
               cliqueClust <- cliqueClust[order(unlist(lapply(lapply(cliqueClust, "[[", "subspace"),"sum")),decreasing = TRUE)]
               nums <- unlist(lapply(lapply(cliqueClust, "[[", "subspace"),"sum"))
               top5 <- nums[1:5]
               top5Metabs <- list()
               top5Diab <- rep(0.5,5)
               
               
               for (i in 1:5){
                 currentCluster <- cliqueClust[as.numeric(i)][[1]]
                 top5Metabs[[i]] <- paste0(names(smallMetab)[currentCluster$subspace], collapse = "_")
                 top5Diab[i] <- mean(T2D[[1]][currentCluster$objects])
               }
               
               #realDiabs <- abs(0.5-top5Diab)
               
               result[counter, ] <- list(paste(method,xi,tau,sep = "_"),paste0(top5Metabs,collapse = ";") ,top5Diab[[1]],top5Diab[[2]],top5Diab[[3]],top5Diab[[4]],top5Diab[[5]])
               counter <- counter + 1
               return(NA)
             },
             error = function(cond){
               log <- c(paste(method,xi,tau,sep = "_"),cond$message,toString(cond$call))
               print(log)
               return(log)
             })
             if (run[[2]]!="keine Funktion abzubrechen, springe zum Top Level"){
               logs <- c(logs,run)
             }
             #currentCombinationCounter <- currentCombinationCounter + 1
           }
      }
  )
}


result <- result[complete.cases(result[[1]]),]
logsToWrite <- data.frame(matrix(unlist(logs), nrow=length(unlist(logs))/3, byrow=T))
write.table(result, file='output.tsv', quote = FALSE, row.names = FALSE, sep='\t')
write.table(logsToWrite, file='logs.tsv', quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
