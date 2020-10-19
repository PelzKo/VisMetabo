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
