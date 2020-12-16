#library(gt)
setwd("C:\\Users\\Konstantin\\Desktop\\Uni\\6Semester\\BachelorArbeit\\paramTuning")

#data <- readFiles("resultsServer\\results")
data <- read.csv("fullParamTuning.tsv",sep = "\t")

doc <- data[grepl("doc",data$name,fixed = TRUE),]
doc <- doc[order(as.numeric(sapply(strsplit(doc$meanEnrichment,"_"),"[[",1)),decreasing = TRUE),]

#docDataset <- getAvgId(doc,2)
#getAvgId(doc,2,TRUE)
#docAlpha <- getAvgId(doc,3)
#docBeta <- getAvgId(doc,4)
#docW <- getAvgId(doc,5)




clique <- data[grepl("clique",data$name,fixed = TRUE),]
clique <- clique[order(as.numeric(sapply(strsplit(clique$meanPheno,"_"),"[[",1)),decreasing = TRUE),]

#cliqueDataset <- getAvgId(clique,2)
#cliqueXi <- getAvgId(clique,3)
#cliqueTau <- getAvgId(clique,4)



som <- data[grepl("som",data$name,fixed = TRUE),]
som <- som[order(as.numeric(sapply(strsplit(som$meanEnrichment,"_"),"[[",1)),decreasing = TRUE),]


#somDataset <- getAvgId(som,2)
#somLayer <- getAvgId(som,3)
#somGrid <- getAvgId(som,4)
#somRLen <- getAvgId(som,5)
#somRAlpha <- getAvgId(som,6)
#lapply(c(2,3,4,5,6),function(x) getAvgId(som,x))


cosa <- data[grepl("cosa",data$name,fixed = TRUE),]
cosa <- cosa[order(as.numeric(sapply(strsplit(cosa$meanEnrichment,"_"),"[[",1)),decreasing = TRUE),]



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
readFiles <- function(pwd){
  files <- list.files(pwd,full.names = TRUE)
  result <- data.frame()
  for (fileName in files){
    tsv <- read.csv(fileName,sep = "\t")
    result <- rbind(result,tsv)
  }
  result
}

tableFromDataFrame <- function(data,v1=3,v2=4,property="maxEnrichment",NumInSplit=1,seperate=NULL){
  attributes <- data.frame(matrix(unlist(strsplit(data$name,"_")), nrow=length(strsplit(data$name,"_")), byrow=T))
  valuesOne <- unique(attributes[,v1])
  valuesTwo <- unique(attributes[,v2])
  result <- data.frame(matrix(rep(NA,length(valuesOne)*length(valuesTwo)),nrow = length(valuesTwo),ncol = length(valuesOne)),row.names = valuesTwo)
  colnames(result) <- valuesOne
  if (!is.null(seperate)){
    valuesSeperate <- unique(attributes[,seperate])
    result <- rep(result,length(valuesSeperate))
    names(result) <- valuesSeperate
  }
  
  for (row in seq_len(nrow(data))){
    line <- data[row,]
    split <- strsplit(line$name,"_")[[1]]
    valOne <- split[v1]
    valTwo <- split[v2]
    value <- as.numeric(strsplit(line[[property]],"_")[[1]][NumInSplit])
    if (is.null(seperate)){
      if (is.na(result[valTwo,valOne])){
        result[valTwo,valOne] <- value
      } else if (!is.na(value)){
        result[valTwo,valOne] <- mean(result[valTwo,valOne],value)
      }
    } else {
      valSeperate <- split[seperate]
      if (is.na(result[valSeperate][valTwo,valOne])){
        result[valSeperate][valTwo,valOne] <- value
      } else if (!is.na(value)){
        result[valSeperate][valTwo,valOne] <- mean(result[valSeperate][valTwo,valOne],value)
      }
    }
  }
  result
  
}
