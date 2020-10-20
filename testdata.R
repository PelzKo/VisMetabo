source("C:\\Users\\Konstantin\\Desktop\\Uni\\6Semester\\BachelorArbeit\\code\\ReadingData.R")

testdata <- function(){
  data <- readFile("C:\\Users\\Konstantin\\Desktop\\Uni\\6Semester\\BachelorArbeit\\BeispielDaten\\metabExampleMaleFemaleSMALLSET50.xlsx")
  #data <- readFile("C:\\Users\\Konstantin\\Desktop\\Uni\\6.Semester\\BachelorArbeit\\data\\QMDiab_aminoGlucose.xlsx")
  inputData <- data$values[complete.cases(data$values), ]
  
  tempId <- 1
  tempStart <- 4
  tempEnd <- 18
  
  
  id <- inputData[[tempId]]
  metab <- inputData[c(tempStart:tempEnd)]
  pheno <- inputData[-c(tempId,c(tempStart:tempEnd))]
  pheno$None <- numeric(length(id[[1]]))
  
  metabComplete <- data.frame(scale(metab))
  metabComplete
}

getData <- function(small=TRUE){
  if (small){
    data <- readFile("C:\\Users\\Konstantin\\Desktop\\Uni\\6Semester\\BachelorArbeit\\BeispielDaten\\paramTuningSmall.xlsx")
    #data <- readFile("C:\\Users\\Konstantin\\Desktop\\Uni\\6.Semester\\BachelorArbeit\\data\\paramTuningSmall.xlsx")
    
    tempId <- 1
    tempStart <- 2
    tempEnd <- 37
  } else {
    data <- readFile("C:\\Users\\Konstantin\\Desktop\\Uni\\6Semester\\BachelorArbeit\\BeispielDaten\\Metabolomics.data.f4.published2011.xlsx")
    #data <- readFile("C:\\Users\\Konstantin\\Desktop\\Uni\\6.Semester\\BachelorArbeit\\data\\Metabolomics.data.f4.published2011.xlsx")
    
    tempId <- 1
    tempStart <- 5
    tempEnd <- 155
  }
  
  inputData <- data$values[complete.cases(data$values), ]
  
  
  
  id <- inputData[[tempId]]
  metab <- inputData[c(tempStart:tempEnd)]
  pheno <- inputData[-c(tempId,c(tempStart:tempEnd))]
  pheno$None <- numeric(length(id[[1]]))
  
  metabComplete <- data.frame(scale(metab))
  res <- list(metabComplete,pheno)
  res
}

testDataFrame <- function(){
  result <- data.frame(one=rep(NA,10),two=rep(NA,10),three=rep(NA,10),four=rep(NA,10),five=rep(NA,10))
  result[1,] <- list(1,1,1,1,1)
  result[2,] <- list(0,0,0,0,0)
  result[3,] <- list(0,0,0,0,0)
  result[4,] <- list(1,1,1,1,1)
  result[5,] <- list(1,1,1,1,1)
  result[6,] <- list(1,1,1,1,1)
  result[7,] <- list(0,0,0,0,0)
  result[8,] <- list(0,0,0,0,0)
  result[9,] <- list(1,1,1,1,1)
  result[10,] <- list(1,1,1,1,1)
  result
}
