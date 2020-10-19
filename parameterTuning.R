source("C:\\Users\\Konstantin\\Desktop\\Uni\\6Semester\\BachelorArbeit\\code\\testdata.R")

small <- getData()
smallMetab <- small[[1]]
smallPheno <- small[[2]]
T2D <- smallPheno["T2D"]
#big <- getData(FALSE)
#bigMetab <- big[[1]]
#bigPheno <- big[[2]]

cliqueClust <- CLIQUE(smallMetab, xi = 5,tau = 0.1)
cliqueClust <- cliqueClust[order(unlist(lapply(lapply(cliqueClust, "[[", "subspace"),"sum")),decreasing = TRUE)]
nums <- unlist(lapply(lapply(cliqueClust, "[[", "subspace"),"sum"))
top5 <- nums[1:5]
top5Metabs <- list()
top5Diab <- rep(0.5,5)


for (i in 1:5){
  currentCluster <- cliqueClust[as.numeric(i)][[1]]
  top5Metabs[[i]] <- names(smallMetab)[currentCluster$subspace]
  top5Diab[i] <- mean(T2D[[1]][currentCluster$objects])
}

realDiabs <- abs(0.5-top5Diab)

if(!is.null(input$clusterId)){
  if (input$clusterId!=0){
    currentCluster <- clusteringData$Clique[as.numeric(input$clusterId)][[1]]
    metabs <- names(metabComplete)[currentCluster$subspace]
    usedMetabs <- sprintf("This cluster was calculed using the following metabolites: <br/>%s", paste(metabs, collapse = '<br/>'))
    output$metabUsed <- renderUI(HTML(usedMetabs))
    
    clusterNumbers[currentCluster$objects] <- 17
    finalValues$idsFromCluster <- currentCluster$objects
  } else {
    output$metabUsed <- renderUI(HTML("No cluster selected"))
  }
}



cliqueMDSValues <- plotFromClusters(lapply(clusteringData$Clique, `[[`, "objects"), returnMDS = TRUE)
clusteringData$CliqueMDS <- data.frame(cbind(seq_len(nrow(cliqueMDSValues)),cliqueMDSValues))
names(clusteringData$CliqueMDS) <- c("id","x","y")