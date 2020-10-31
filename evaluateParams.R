setwd("C:\\Users\\Konstantin\\Desktop\\Uni\\6Semester\\BachelorArbeit\\paramTuning")

cliqueSom <- read.csv("outputCliqueSom.tsv",sep = "\t")
doc <- read.csv("outputDoc.tsv",sep = "\t")


doc <- doc[order(doc$meanDist,decreasing = TRUE),]

docDataset <- getAvgId(doc,2)
getAvgId(doc,2,TRUE)
docAlpha <- getAvgId(doc,3)
docBeta <- getAvgId(doc,4)
docW <- getAvgId(doc,5)




clique <- cliqueSom[grepl("clique",cliqueSom$name,fixed = TRUE),]
clique <- clique[order(clique$meanDist,decreasing = TRUE),]

cliqueDataset <- getAvgId(clique,2)
cliqueXi <- getAvgId(clique,3)
cliqueTau <- getAvgId(clique,4)


som <- cliqueSom[grepl("som",cliqueSom$name,fixed = TRUE),]
som <- som[order(som$meanDist,decreasing = TRUE),]
som$name <- gsub("(.*)_([:digit:]*.*[:digit:]*)","\\1;\\2",som$name)


somDataset <- getAvgId(som,2)
somLayer <- getAvgId(som,3)
somGrid <- getAvgId(som,4)
somRLen <- getAvgId(som,5)
somRAlpha <- getAvgId(som,6)
#lapply(c(2,3,4,5,6),function(x) getAvgId(som,x))


rowsWith <- function(data,propertyNum,justBools=FALSE){
  if (length(propertyNum)==1){
    valueTable <- data.frame(matrix(unlist(strsplit(data$name,"_")), nrow=length(strsplit(data$name,"_")), byrow=T))
    uniques <- unique(valueTable[,propertyNum])
    contains <- lapply(uniques,function(x) valueTable[,propertyNum]==x)
    names(contains) <- uniques
    if (justBools){
      return(contains)
    }
    return(lapply(contains,function(x) data[x,]))
  } 
}

getAvgId <- function(data,propNum,numInd=FALSE){
  ids <- lapply(rowsWith(data,propNum,TRUE),function(x) seq_len(length(data$name))[x])
  if (numInd){
    return(length(ids))
  }
  return(lapply(ids,function(x) sum(x)/length(x)))
}

makeDf <- function(x,name){
  result <- data.frame(matrix(unlist(docBeta),byrow=T),row.names = names(x))
  colnames(result) <- name
  result
}
