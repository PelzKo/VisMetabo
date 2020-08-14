library("openxlsx")
#dftest <- "C:/Users/Konstantin/Desktop/Uni/6Semester/BachelorArbeit/BeispielDaten/metabExampleMaleFemale.xlsx"
#df <- read.xlsx("C:/Users/Konstantin/Desktop/Uni/6Semester/BachelorArbeit/BeispielDaten/QMDiab_metabolomics_Preprocessed.xlsx", sheet = 1)

readFile <- function(filePath, fileSheet = 1){
  df <- read.xlsx(filePath, sheet = fileSheet)
  df <- df[complete.cases(df),]
  
  id_col <- 0
  need_id <- TRUE
  metab_start <- 1
  metab_end <- 1
  metab_start_final <- 1
  metab_end_final <- 1
  
  #iterating through all columns
  for (col in seq_len(ncol(df))){
    #if column is unique and we have not found and id column yet
    if (need_id&&length(df[[col]])==length(unique(df[[col]]))){
      id_col<-col
      need_id<-FALSE
    }
    if (typeof(df[[col]])=="double"){
      #append the current metabolite interval if the last column was a double value as well
      if (col==metab_end+1){
        metab_end<-col
      }
      else {
        #if longer than current longest metabolites
        if (metab_end-metab_start>metab_end_final-metab_start_final){
          metab_start_final<-metab_start
          metab_end_final<-metab_end
        }
        metab_start<-col
        metab_end<-col
      }
    }
    #columns start with a non-double column, otherwise it would be included in first column interval
    else if (col==1){
      metab_start <- 0
      metab_end <- -1
      metab_start_final <- 0
      metab_end_final <- -1
    }
  }
  
  #check once at the end
  if (metab_end-metab_start>metab_end_final-metab_start_final){
    metab_start_final<-metab_start
    metab_end_final<-metab_end
  }
  
  if (id_col==0){
    stop("NO ID COLUMN FOUND")
  }
  if (id_col>=metab_start_final&&id_col<=metab_end_final){
    stop("ID COLUMN FOUND IN METABOLITES")
  }
  
  
  print("Done")
  
  #metab <- df[c(metab_start_final:metab_end_final)]
  #pheno <- df[-c(metab_start_final:metab_end_final)]
  #return(list(metab,pheno))
  return (list(values = df,id = id_col,metab_start = metab_start_final,metab_end = metab_end_final))
}

readPhenoFile <- function(filePath, fileSheet = 1){
  df <- read.xlsx(filePath, sheet = fileSheet)
  df <- df[complete.cases(df),]
 
  return (df)
}
