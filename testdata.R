testdata <- function(){
  #data <- readFile("C:\\Users\\Konstantin\\Desktop\\Uni\\6Semester\\BachelorArbeit\\BeispielDaten\\metabExampleMaleFemale.xlsx")
  data <- readFile("C:\\Users\\Konstantin\\Desktop\\Uni\\6.Semester\\BachelorArbeit\\data\\QMDiab_aminoGlucose.xlsx")
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
