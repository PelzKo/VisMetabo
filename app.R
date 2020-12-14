library(shiny)
library(shinyBS)
library(subspace)
library(ggplot2)
library(ggfortify)
library(rCOSA)
library(kohonen)
library(rJava)
library(RANN)
library(dummies)
library(ggdendro)

source("ReadingData.R")
source("utility.R")
source("smacofFixed.R")
source("doc.R")

# Allow Uploads until 200MB
options(shiny.maxRequestSize=200*1024^2)
options(java.parameters = "-Xmx16000m")


# Define UI for app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Analysis of Metabolite data"),
  
  tabsetPanel(id = "tabs", type = "tabs",
    # Input Panel
    tabPanel(title = "Data Input", value = "0",
      titlePanel("Please upload your data"),
      "Important things about the data:",
      "- Rows are individuals, Columns are features",
      fileInput(inputId = "inputData",
                label = "Data:",
                accept = c(".xlsx",".csv",".tsv")),
      checkboxInput(inputId = "externalPheno",
                    label = "Upload phenotypes from another file"),
      uiOutput(outputId = "phenoInputField"),
      actionButton("toTypes", "Upload")
             
    ),
    # Input Validation Panel
    tabPanel(title = "Specifying Column Types", value = "1",
      selectInput(inputId = "idField",
                 label = "Unique ID:",
                 choices = list()),
      textOutput(outputId = "idError"),
      selectInput(inputId = "metabStart",
                 label = "First Metabolite:",
                 choices = list()),
      selectInput(inputId = "metabEnd",
                 label = "Last Metabolite:",
                 choices = list()),
      actionButton("toOut", "Cluster Now")
    ),
    # Output Panel
    tabPanel(title = "Output", value = "2",
      # Sidebar layout with input and output definitions ----
      sidebarLayout(
        
        # Sidebar panel for inputs ----
        sidebarPanel(
          
          # Input: Slider for the number of bins ----
          selectInput(inputId = "selectedPhenotype",
                      label = "Phenotype:",
                      choices = list()),
          selectInput(inputId = "clusteringType",
                      label = "Clustering Algorithm used:",
                      choices = list("Nothing","Clique","SOM", "COSA","DOC")),
          actionButton("resetClustering", "Reset the Clustering"),
          checkboxInput("pcaSwitch", "Display the Clusters in the PCA Plot"),
          uiOutput(outputId = "min"),
          uiOutput(outputId = "max"),
          checkboxInput(inputId = "reverse",
                        label = "Reverse Color Strength"),
          checkboxInput(inputId = "removeOutliers",
                        label = "Remove the highest and lowest 5% (Only of number phenotypes)"),
          textOutput(outputId = "phenoInfo"),
          plotOutput(outputId = "phenoHist"),
          
          
          width = 3
        ),
        
        # Main panel for displaying outputs ----
        mainPanel(
          
          bsModal("cosaHist", "Histogram of the rCosa", "go", size = "large",plotOutput("hist", click = "clickHist")),
          uiOutput("clusteringPlot"),
          uiOutput("pcaAll"),
          uiOutput("downloadButton"),
          htmlOutput("info"),
          
          width = 9
          
        )
      )
    )
  )
)


# Define server logic required ----
server <- function(input, output, session) {
  plotsSize <- "1000px"
  # Data uploaded by the user
  data <- reactive({
    if (!is.null(input$inputData)){
      readFile(input$inputData$datapath)
    }
    else{
      list("Error" = "Please upload a file first!")
    }
  })
  # Transformed/normalized data
  finalValues <- reactiveValues(id = 1, metab = list(), pheno = list())
  # Clusters
  clusteringData <- reactiveValues(Clique = list(), CliqueMDS = list(), SOM = list(), COSA = list(), DOC = list(), DocMDS = list(), PhenoPValues = list())
  # PCAData
  pcaValues <- reactiveValues()
  #Temp Values to distinguish which click/brush changed
  temps <- reactiveValues(clusteringClick = NULL,clusteringBrush = NULL, pcaBrush = NULL, idsFromClusterOld = NULL, index = list(),
                          oldk = NULL, oldx = NULL, hc = NULL, grps = c(0), id = 1, clickHist = NULL, rects = data.frame(x1=NULL,y1=NULL,x2=NULL,y2=NULL,border=NULL))
  
  
  output$pca <- renderPlot({
    if (length(data())==4&input$clusteringType!="Nothing"){
      firstRunForClusteringMethod <- length(clusteringData[[input$clusteringType]])==0
      firstPCA <- length(pcaValues$visualisation)==0
      
      colorClusterNotPheno <- FALSE
      if (!is.null(input$colorClusterNoPheno)){
        colorClusterNotPheno <- input$colorClusterNoPheno[[1]]
      }
      
      
      metabComplete <- data.frame(scale(finalValues$metab))
      
      if (firstRunForClusteringMethod){
        switch(input$clusteringType, 
               SOM={
                 clusteringData$SOM <- som(data.matrix(metabComplete), grid = somgrid(xdim = 5,ydim = 5, topo = "hexagonal"),rlen = 100, alpha = c(0.08,0.01))
               },
               COSA={
                 clusteringData$COSA <- cosa2(metabComplete, lambda = 0.4, niter = 2, noit = 30)
                 clusteringData$COSA$smacof <- smacof(clusteringData$COSA$D, niter = 30, interc = 1, VERBOSE = FALSE, PLOT = FALSE)
                 #clusteringData$COSA <- cosa2(metabComplete, niter = 1, noit = 1)
                 #clusteringData$COSA$smacof <- smacof(clusteringData$COSA$D, niter = 1, interc = 1, VERBOSE = FALSE, PLOT = FALSE)
                 clusteringData$COSA$idsAndSmacof <- data.frame(cbind(names(finalValues$idFromNum),clusteringData$COSA$smacof$X))
                 names(clusteringData$COSA$idsAndSmacof) <- c("id","x","y")
                 
                 toggleModal(session, "cosaHist", toggle = "toggle")
                   
                 
               },
               DOC={
                 clusteringData$DOC <- tryCatch(runDoc(metabComplete, 0.1, 0.75, calcW(metabComplete, 1.3)),error = function(x){return(list())})
                 
                 docMDSValues <- plotFromClusters(getIdsDoc(clusteringData$DOC), returnMDS = TRUE)
                 clusteringData$DocMDS <- data.frame(cbind(seq_len(nrow(docMDSValues)),docMDSValues))
                 names(clusteringData$DocMDS) <- c("id","x","y")
               },
               {
                 # default is using Clique
                 inds <- nrow(metabComplete)
                 xi <- 8
                 if (inds>800){
                   xi <- xi + floor((inds-800)/200)
                   if (xi>14){
                     xi <- 14
                   }
                 }
                 clusteringData$Clique <- CLIQUE(metabComplete, xi = xi,tau = 0.1) 
                 clusteringData$Clique <- clusteringData$Clique[order(unlist(lapply(lapply(clusteringData$Clique, "[[", "subspace"),"sum")),decreasing = TRUE)]
                 
                 cliqueMDSValues <- plotFromClusters(lapply(clusteringData$Clique, `[[`, "objects"), returnMDS = TRUE)
                 clusteringData$CliqueMDS <- data.frame(cbind(seq_len(nrow(cliqueMDSValues)),cliqueMDSValues))
                 names(clusteringData$CliqueMDS) <- c("id","x","y")
                 
               }
        )
      }
      
      origPheno <- finalValues$pheno[[as.numeric(input$selectedPhenotype)]]
      phenotype <- origPheno
      phenotypeNoNa <- phenotype[!is.na(phenotype)]
      
      
      if (class(phenotype)=="character"){
        origPheno <- factor(origPheno)
        phenotype <- as.numeric(factor(phenotype))
        phenotypeNoNa <- as.numeric(factor(phenotypeNoNa))
      } else {
        if (input$removeOutliers[[1]]){
          n <- 5
          lowestTooHigh <- min(phenotypeNoNa[phenotypeNoNa >= quantile(phenotypeNoNa,prob=1-n/100)])
          highestTooLow <- max(phenotypeNoNa[phenotypeNoNa <= quantile(phenotypeNoNa,prob=n/100)])
          phenotype[phenotype>lowestTooHigh] <- lowestTooHigh
          phenotype[phenotype<highestTooLow] <- highestTooLow
        }
      }
      
      if (length(unique(finalValues$pheno[[as.numeric(input$selectedPhenotype)]]))==1&&finalValues$pheno[[as.numeric(input$selectedPhenotype)]][[1]]==0){
        coloredPoints <- rgb(0, 0, 0, maxColorValue=255, alpha=255)
        
        
        output$phenoInfo <- NULL
        output$phenoHist <- NULL
        
      } else {
        colorPalette <- getColors(length(unique(phenotypeNoNa)))
        if (input$reverse[[1]]){
          output$min <- renderUI(HTML(sprintf(paste("<div style='background-color: ",colorPalette[1],";height: 18px;width: 18px;float: left;margin-right: 3px;'></div> Min: %s",sep = ""), fourDeci(min(phenotypeNoNa)))))
          output$max <- renderUI(HTML(sprintf(paste("<div style='background-color: ",colorPalette[length(colorPalette)],";height: 18px;width: 18px;float: left;margin-right: 3px;'></div> Max: %s",sep = ""), fourDeci(max(phenotypeNoNa)))))
          coloredPoints <- vecToCol(phenotype,colorPalette)
        } else {
          colorPalette <- rev(colorPalette)
          output$min <- renderUI(HTML(sprintf(paste("<div style='background-color: ",colorPalette[1],";height: 18px;width: 18px;float: left;margin-right: 3px;'></div> Min: %s",sep = ""), fourDeci(min(phenotypeNoNa)))))
          output$max <- renderUI(HTML(sprintf(paste("<div style='background-color: ",colorPalette[length(colorPalette)],";height: 18px;width: 18px;float: left;margin-right: 3px;'></div> Max: %s",sep = ""), fourDeci(max(phenotypeNoNa)))))
          coloredPoints <- vecToCol(phenotype,colorPalette)
        }
        output$phenoInfo <- renderText(sprintf("Mean of phenotypes after normalization between 0 and 1: %s", fourDeci(mean(range01(phenotypeNoNa)))))
        if (class(origPheno)=="factor"){
          #print(colorPalette)
          output$phenoHist <- renderPlot({
                                hist(phenotypeNoNa,xlab="Phenotype as factors",main="Frequency of levels", col=colorPalette,breaks=seq(0,length(levels(origPheno))))
                                legend("topright",levels(origPheno),col=colorPalette,pch=15,pt.cex = 2)
                              })
        } else {
          output$phenoHist <- renderPlot(hist(phenotype))
        }
      }
      
      clusterNumbers <- rep(16,nrow(metabComplete))
      clusterColors <- rep("blue",nrow(metabComplete))
      switch(input$clusteringType,
             SOM={
               output$clusteringPlot <- renderUI(plotOutput(outputId = "clustering", click = "clustering_click"))
               #output$clustering <- renderPlot(plot(clusteringData$SOM, type="mapping", classif=predict(clusteringData$SOM)
               #                                     , pchs = c(1,2,3,4,5)))
               if (as.numeric(input$selectedPhenotype)==length(finalValues$pheno)){
                 output$clustering <- renderPlot(plotHeatMap(clusteringData$SOM,finalValues$pheno,0))
               } else {
                 output$clustering <- renderPlot(plotHeatMap(clusteringData$SOM,finalValues$pheno,as.numeric(input$selectedPhenotype)))
               }
               
               #clusters <- as.character(sort(unique(clusteringData$SOM$unit.classif)))
               #colorClusters <- rainbow(length(clusters))
               #names(colorClusters) <- clusters
               #clusterColors <- clusteringData$SOM$unit.classif
               #for (num in clusters){
               #  clusterColors <- replace(clusterColors, clusterColors==as.numeric(num), colorClusters[[num]])
               #}
               if(!is.null(input$clusterId)){
                 if (input$clusterId!=0){
                   clusterNumbers[clusteringData$SOM$unit.classif==input$clusterId] <- 17
                   clusterColors[clusteringData$SOM$unit.classif==input$clusterId] <- "red"
                   currentCluster <- getCurrentClusters()[[as.numeric(input$clusterId)]]
                   enrichedPhenos <- ""
                   for (i in seq_len(length(finalValues$pheno)-1)){
                     currentPheno <- finalValues$pheno[[i]]
                     pVal <- calcPValueForCluster(currentCluster,currentPheno)
                     if (pVal<0.05){
                       enrichedPhenos <- sprintf("%s%s with a pValue of %s<br/>",enrichedPhenos,names(finalValues$pheno)[[i]],pVal)
                     }
                   }
                   numberOfClusters <- length(unique(clusteringData$SOM$unit.classif))
                   enrichedPhenos <- sprintf("%s<br/>Bonferroni cut-off with %s clusters would be %s<br/>",enrichedPhenos,numberOfClusters,0.05/numberOfClusters)
                   if (!is.na(clusteringData$PhenoPValues[[input$clusteringType]][[as.numeric(input$selectedPhenotype)]])){
                     sumNormal <- sum(clusteringData$PhenoPValues[[input$clusteringType]][[as.numeric(input$selectedPhenotype)]]<0.05,na.rm = TRUE)
                     sumBonferroni <- sum(clusteringData$PhenoPValues[[input$clusteringType]][[as.numeric(input$selectedPhenotype)]]<(0.05/numberOfClusters),na.rm = TRUE)
                     enrichedPhenos <- sprintf("%s%s Clusters are significant (for this phenotype) with the cutoff <0.05<br/>%s Clusters are significant (for this phenotype) with the Bonferroni cutoff<br/>",enrichedPhenos,sumNormal,sumBonferroni)
                   }
                   output$metabUsed <- renderUI(HTML(sprintf("Cluster %s selected<br/>Significant Phenotypes:<br/>%s",input$clusterId,enrichedPhenos)))
                   finalValues$idsFromCluster <- finalValues$numFromId[clusteringData$SOM$unit.classif==input$clusterId]
                 } else {
                   output$metabUsed <- renderUI(HTML("No cluster selected"))
                 }
               }
               
             },
             COSA={
               output$clusteringPlot <- renderUI(plotOutput(outputId = "clustering", brush = "clustering_brush", height = plotsSize, width = plotsSize))
               
               if(!is.null(input$clusterId)){
                 if (input$clusterId!=0){
                   currentCluster <- temps$index[[as.numeric(input$clusterId)]]
                   clusterNumbers[currentCluster] <- 17
                   clusterColors[currentCluster] <- "red"
                   enrichedPhenos <- ""
                   for (i in seq_len(length(finalValues$pheno)-1)){
                     currentPheno <- finalValues$pheno[[i]]
                     pVal <- calcPValueForCluster(currentCluster,currentPheno)
                     if (pVal<0.05){
                       enrichedPhenos <- sprintf("%s%s with a pValue of %s<br/>",enrichedPhenos,names(finalValues$pheno)[[i]],pVal)
                     }
                   }
                   
                   numberOfClusters <- length(temps$index)
                   enrichedPhenos <- sprintf("%s<br/>Bonferroni cut-off with %s clusters would be %s<br/>",enrichedPhenos,numberOfClusters,0.05/numberOfClusters)
                   if (!is.na(clusteringData$PhenoPValues[[input$clusteringType]][[as.numeric(input$selectedPhenotype)]])){
                     sumNormal <- sum(clusteringData$PhenoPValues[[input$clusteringType]][[as.numeric(input$selectedPhenotype)]]<0.05,na.rm = TRUE)
                     sumBonferroni <- sum(clusteringData$PhenoPValues[[input$clusteringType]][[as.numeric(input$selectedPhenotype)]]<(0.05/numberOfClusters),na.rm = TRUE)
                     enrichedPhenos <- sprintf("%s%s Clusters are significant (for this phenotype) with the cutoff <0.05<br/>%s Clusters are significant (for this phenotype) with the Bonferroni cutoff<br/>",enrichedPhenos,sumNormal,sumBonferroni)
                   }
                   output$metabUsed <- renderUI(HTML(sprintf("Cluster %s selected<br/>Significant Phenotypes:<br/>%s",input$clusterId,enrichedPhenos)))
                   finalValues$idsFromCluster <- currentCluster
                 } else {
                   output$metabUsed <- renderUI(HTML("No cluster selected"))
                 }
               }
               
               
               if(colorClusterNotPheno){
                 coloredPoints <- clusterColors
               }
               
               output$clustering <- renderPlot(plotSmacof(clusteringData$COSA$smacof[["X"]], cols = coloredPoints, pch = clusterNumbers))
             },
             DOC={
               output$clusteringPlot <- renderUI(plotOutput(outputId = "clustering", brush = "clustering_brush", height = plotsSize, width = plotsSize))
               
               idsInClustersDoc <- getIdsDoc(clusteringData$DOC)
               
               if(!is.null(input$clusterId)){#&input$pcaSwitch[[1]]){
                 if (input$clusterId!=0){
                   currentCluster <- idsInClustersDoc[[as.numeric(input$clusterId)]]
                   enrichedPhenos <- ""
                   for (i in seq_len(length(finalValues$pheno)-1)){
                     currentPheno <- finalValues$pheno[[i]]
                     pVal <- calcPValueForCluster(currentCluster,currentPheno)
                     if (pVal<0.05){
                       enrichedPhenos <- sprintf("%s%s with a pValue of %s<br/>",enrichedPhenos,names(finalValues$pheno)[[i]],pVal)
                     }
                   }
                   
                   numberOfClusters <- length(idsInClustersDoc)
                   enrichedPhenos <- sprintf("%s<br/>Bonferroni cut-off with %s clusters would be %s<br/>",enrichedPhenos,numberOfClusters,0.05/numberOfClusters)
                   if (!is.na(clusteringData$PhenoPValues[[input$clusteringType]][[as.numeric(input$selectedPhenotype)]])){
                     sumNormal <- sum(clusteringData$PhenoPValues[[input$clusteringType]][[as.numeric(input$selectedPhenotype)]]<0.05,na.rm = TRUE)
                     sumBonferroni <- sum(clusteringData$PhenoPValues[[input$clusteringType]][[as.numeric(input$selectedPhenotype)]]<(0.05/numberOfClusters),na.rm = TRUE)
                     enrichedPhenos <- sprintf("%s%s Clusters are significant (for this phenotype) with the cutoff <0.05<br/>%s Clusters are significant (for this phenotype) with the Bonferroni cutoff<br/>",enrichedPhenos,sumNormal,sumBonferroni)
                   }
                   
                   usedDimensions <- getDimsDoc(clusteringData$DOC)[[as.numeric(input$clusterId)]]
                   avgs <- getAvgsDoc(clusteringData$DOC)[as.numeric(input$clusterId),]
                   metabs <- names(metabComplete)[usedDimensions]
                   usedMetabs <- sprintf("This cluster was calculed using the following metabolites: <br/>%s<br/><br/>Significant Phenotypes:<br/>%s", paste(metabs, collapse = '<br/>'),enrichedPhenos)
                   averagesFormatted <- mapply(function(x,y) paste(x, round(as.numeric(y), digits=4), sep=": "), metabs, avgs, SIMPLIFY=FALSE)
                   metabAvgs <- sprintf("<br/><br/>The metabolites have the following averages: <br/>%s", paste(averagesFormatted, collapse = '<br/>'))
                   output$metabUsed <- renderUI(HTML(paste(usedMetabs,metabAvgs)))
                   
                   clusterNumbers[currentCluster] <- 17
                   clusterColors[currentCluster] <- "red"
                   finalValues$idsFromCluster <- currentCluster
                 } else {
                   output$metabUsed <- renderUI(HTML("No cluster selected"))
                 }  
               }
               
               if(colorClusterNotPheno){
                 coloredPoints <- clusterColors
               }
               output$clustering <- renderPlot(plotFromClusters(idsInClustersDoc, label = FALSE, colors = coloredPoints, pch = clusterNumbers))
             },
             {
               # default is using Clique
               output$clusteringPlot <- renderUI(plotOutput(outputId = "clustering", brush = "clustering_brush", height = plotsSize, width = plotsSize))
               
               if(!is.null(input$clusterId)){
                 if (input$clusterId!=0&input$clusterId!=""){
                   currentCluster <- clusteringData$Clique[as.numeric(input$clusterId)][[1]]
                   enrichedPhenos <- ""
                   for (i in seq_len(length(finalValues$pheno)-1)){
                     currentPheno <- finalValues$pheno[[i]]
                     #print(names(finalValues$pheno)[[i]])
                     pVal <- calcPValueForCluster(currentCluster$objects,currentPheno)
                     if (pVal<0.05){
                       enrichedPhenos <- sprintf("%s%s with a pValue of %s<br/>",enrichedPhenos,names(finalValues$pheno)[[i]],pVal)
                     }
                   }
                   numberOfClusters <- length(clusteringData$Clique)
                   enrichedPhenos <- sprintf("%s<br/>Bonferroni cut-off with %s clusters would be %s<br/>",enrichedPhenos,numberOfClusters,0.05/numberOfClusters)
                   if (!is.na(clusteringData$PhenoPValues[[input$clusteringType]][[as.numeric(input$selectedPhenotype)]])){
                     sumNormal <- sum(clusteringData$PhenoPValues[[input$clusteringType]][[as.numeric(input$selectedPhenotype)]]<0.05,na.rm = TRUE)
                     sumBonferroni <- sum(clusteringData$PhenoPValues[[input$clusteringType]][[as.numeric(input$selectedPhenotype)]]<(0.05/numberOfClusters),na.rm = TRUE)
                     enrichedPhenos <- sprintf("%s%s Clusters are significant (for this phenotype) with the cutoff <0.05<br/>%s Clusters are significant (for this phenotype) with the Bonferroni cutoff<br/>",enrichedPhenos,sumNormal,sumBonferroni)
                   }
                   metabs <- names(metabComplete)[currentCluster$subspace]
                   usedMetabs <- sprintf("This cluster was calculed using the following metabolites: <br/>%s<br/><br/>Significant Phenotypes:<br/>%s", paste(metabs, collapse = '<br/>'),enrichedPhenos)
                   output$metabUsed <- renderUI(HTML(usedMetabs))
                   
                   clusterNumbers[currentCluster$objects] <- 17
                   clusterColors[currentCluster$objects] <- "red"
                   finalValues$idsFromCluster <- currentCluster$objects
                 } else {
                   output$metabUsed <- renderUI(HTML("No cluster selected"))
                 }
               }
               
               if(colorClusterNotPheno){
                 coloredPoints <- clusterColors
               }
               
               output$clustering <- renderPlot(plotFromClusters(lapply(clusteringData$Clique, `[[`, "objects"), label = FALSE, colors = coloredPoints, pch = clusterNumbers))
               
             }
      )
      
      
      if (firstPCA){
        pca_data <- prcomp(metabComplete, scale. = TRUE)
        ## Let us calculat the variances covered by components.
        pca_data_perc=round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)
        pcaValues$percentage <- pca_data_perc
        
        ## create a data frame with principal component 1 (PC1), PC2, Conditions and sample names
        df_pca_data = data.frame(PC1 = pca_data$x[,1], PC2 = pca_data$x[,2])
        pcaValues$visualisation <- df_pca_data
        pcaValues$visWithId <- data.frame(cbind(names(finalValues$idFromNum),pcaValues$visualisation))
        names(pcaValues$visWithId) <- c("id","x","y")
        
      }
      
      coloredPointsPCA <- coloredPoints
      
      if(colorClusterNotPheno){
        coloredPointsPCA <- clusterColors
      }
      
      plotNoLims(pcaValues$visualisation, "PCA - PC1 vs PC2", sprintf("PC1 (%s%%)",pcaValues$percentage[[1]]),sprintf("PC2 (%s%%)",pcaValues$percentage[[2]]), cols = coloredPointsPCA, pch = clusterNumbers)
    }
  })
  
  #height = function() {session$clientData$output_pca_width}
  
  # Goto Input Validation Panel
  observeEvent(input$toTypes, {
    if (length(data())==4){
      columnNames <- colnames(data()$values)
      columns <- seq_len(length(columnNames))
      names(columns)<-columnNames
      updateSelectInput(session, "idField", choices = columns, selected = columns[data()$id])
      updateSelectInput(session, "metabStart", choices = columns, selected = columns[data()$metab_start])
      updateSelectInput(session, "metabEnd", choices = columns, selected = columns[data()$metab_end])
    }
    
    updateTabsetPanel(session, "tabs",
                      selected = "1")
    output$pcaAll <- renderUI({
      plotOutput(outputId = "pca",
                 brush = "pca_brush", height = plotsSize, width = plotsSize)
    })
  })
  
  # Check if ID Field is unique
  observeEvent(input$idField, {
    if (length(data())==4){
      dataNoNa <- data()$values[complete.cases(data()$values), ]
      tempId <- as.numeric(input$idField)
      if (length(dataNoNa[[tempId]])!=length(unique(dataNoNa[[tempId]]))){
        output$idError <- renderText("Identifier is not unique")
      } else {
        output$idError <- NULL
      }
    }
    
  })
  
  # Reset the saved Clusterings
  observeEvent(input$resetClustering, {
      switch(input$clusteringType, 
             SOM={
               clusteringData$SOM = list()
             },
             COSA={
               clusteringData$COSA = list()
               temps$index = list()
               temps$oldk = NULL
               temps$oldx = NULL
               temps$hc = NULL
               temps$grps = c(0)
               temps$clickHist = NULL
               temps$rects = data.frame(x1=NULL,y1=NULL,x2=NULL,y2=NULL,border=NULL)
               },
             DOC={  
               clusteringData$DOC = list()
             },
             Clique={
               clusteringData$Clique = list()
             }
      )
  })
      
  # Goto Output Panel
  observeEvent(input$toOut, {
    if (length(data())==4){
      clusteringData$Clique = list()
      clusteringData$SOM = list()
      clusteringData$COSA = list()
      clusteringData$DOC = list()
      
      tempId <- as.numeric(input$idField)
      tempStart <- as.numeric(input$metabStart)
      tempEnd <- as.numeric(input$metabEnd)
      
      
      inputData <- data()$values[complete.cases(data()$values[c(tempStart:tempEnd)]), ]
      
      if (length(inputData[[tempId]])!=length(unique(inputData[[tempId]]))){
        output$idError <- renderText("Identifier is not unique")
        return()
      }
      
      if (input$externalPheno[[1]]&&!is.null(input$inputPheno)){
        dataTemp <- merge(inputData, readPhenoFile(input$inputPheno$datapath), by = names(inputData)[[tempId]], all.x = TRUE)
        if (nrow(dataTemp)==0){
          showNotification("Could not match any metab ids with phenotype ids. Please check you have the correct 
                           ID column, they are named the same and the values are written in the same format.",type = "warning")
          showNotification("Continuing with the phenotypes from the metabolite data sheet.",type = "message")
        } else {
          inputData <- dataTemp
        }
      }
      
      finalValues$tempInfo <- "Nothing clicked/brushed yet"
      finalValues$completeData <- inputData
      finalValues$id <- inputData[[tempId]]
      
      finalValues$numFromId <- seq_len(length(inputData[[tempId]]))
      names(finalValues$numFromId) <- inputData[[tempId]]
      finalValues$idFromNum <- inputData[[tempId]]
      names(finalValues$idFromNum) <- seq_len(length(inputData[[tempId]]))
      
      finalValues$metab <- inputData[c(tempStart:tempEnd)]
      finalValues$pheno <- inputData[-c(tempId,c(tempStart:tempEnd))]
      finalValues$pheno$None <- numeric(length(finalValues$id[[1]]))
      
      columnNames <- names(finalValues$pheno)
      #clusteringData$PhenoPValues <- as.list(rep(NA,length(columnNames)))
      columns <- seq_len(length(columnNames))
      names(columns)<-columnNames
      if (length(columns)>0){
        updateSelectInput(session, "selectedPhenotype", choices = columns, selected = columns[["None"]])
      }
      
      
    }
    updateTabsetPanel(session, "tabs",
                      selected = "2")
  })
  
  #Take phenotypes from different file
  output$phenoInputField <- renderUI({
    if (input$externalPheno[[1]]){
      return(fileInput(inputId = "inputPheno",
                       label = "Phenotypes (they need to have the same unique id as the metab data):",
                       accept = c(".xlsx",".csv",".tsv")))
    } else {
      return(NULL)
    }
  })
  
  
  output$info <- renderUI({
    clusteringClick <- input$clustering_click
    clusteringBrush <- input$clustering_brush
    pcaBrush <- input$pca_brush
    idsFromCluster <- finalValues$idsFromCluster
    
    if (input$clusteringType=="Nothing"){
      return(HTML(finalValues$tempInfo))
    }
    
    if (!is.null(clusteringBrush)&!equalsBrush(temps$clusteringBrush,clusteringBrush)){
      switch(input$clusteringType, 
             SOM={
             },
             COSA={
               points <- brushedPoints(clusteringData$COSA$idsAndSmacof, clusteringBrush, xvar = "x", yvar = "y")
               finalValues$currentIds <- finalValues$idFromNum[as.character(points$id)]
               phenotypes <- finalValues$pheno[[as.numeric(input$selectedPhenotype)]]
               phenoInPoints <- phenotypes[points$id]
               phenotypeAverage <- mean(phenoInPoints[!is.na(phenoInPoints)])
               averageSelected <- colMeans(finalValues$metab[points$id,])
               
               ids <- sprintf("The area you selected (phenotype average of %s), contains the following ids: <br/>%s", round(phenotypeAverage, digits = 2), paste(finalValues$currentIds, collapse = ', '))
               averagesFormatted <- mapply(function(x,y) paste(x, round(as.numeric(y), digits=4), sep=": "), names(finalValues$metab), averageSelected, SIMPLIFY=FALSE)
               average <- sprintf("The average values in this area are: <br/>%s", paste(averagesFormatted, collapse = '<br/>'))
               
               finalValues$tempInfo <- paste(ids, average, sep = '<br/>')
             },
             DOC={
               points <- brushedPoints(clusteringData$DocMDS, clusteringBrush, xvar = "x", yvar = "y")
               finalValues$currentIds <- finalValues$idFromNum[as.character(points$id)]
               phenotypes <- finalValues$pheno[[as.numeric(input$selectedPhenotype)]]
               phenoInPoints <- phenotypes[points$id]
               phenotypeAverage <- mean(phenoInPoints[!is.na(phenoInPoints)])
               averageSelected <- colMeans(finalValues$metab[points$id,])
               
               ids <- sprintf("The area you selected (phenotype average of %s), contains the following ids: <br/>%s", round(phenotypeAverage, digits = 2), paste(finalValues$currentIds, collapse = ', '))
               averagesFormatted <- mapply(function(x,y) paste(x, round(as.numeric(y), digits=4), sep=": "), names(finalValues$metab), averageSelected, SIMPLIFY=FALSE)
               average <- sprintf("The average values in this area are: <br/>%s", paste(averagesFormatted, collapse = '<br/>'))
               
               finalValues$tempInfo <- paste(ids, average, sep = '<br/>')
             },
             {
               # default is using Clique
               points <- brushedPoints(clusteringData$CliqueMDS, clusteringBrush, xvar = "x", yvar = "y")
               finalValues$currentIds <- finalValues$idFromNum[as.character(points$id)]
               phenotypes <- finalValues$pheno[[as.numeric(input$selectedPhenotype)]]
               phenoInPoints <- phenotypes[points$id]
               phenotypeAverage <- mean(phenoInPoints[!is.na(phenoInPoints)])
               averageSelected <- colMeans(finalValues$metab[points$id,])
               
               ids <- sprintf("The area you selected (phenotype average of %s), contains the following ids: <br/>%s", round(phenotypeAverage, digits = 2), paste(finalValues$currentIds, collapse = ', '))
               averagesFormatted <- mapply(function(x,y) paste(x, round(as.numeric(y), digits=4), sep=": "), names(finalValues$metab), averageSelected, SIMPLIFY=FALSE)
               average <- sprintf("The average values in this area are: <br/>%s", paste(averagesFormatted, collapse = '<br/>'))
               
               finalValues$tempInfo <- paste(ids, average, sep = '<br/>')
             }
      )
      
    } else if (!is.null(pcaBrush)&!equalsBrush(temps$pcaBrush,pcaBrush)){
      pointsPCA <- brushedPoints(pcaValues$visWithId, pcaBrush, xvar = "x", yvar = "y")
      finalValues$currentIds <- finalValues$idFromNum[as.character(pointsPCA$id)]
      phenotypes <- finalValues$pheno[[as.numeric(input$selectedPhenotype)]]
      phenoInPoints <- phenotypes[pointsPCA$id]
      phenotypeAverage <- mean(phenoInPoints[!is.na(phenoInPoints)])
      averageSelected <- colMeans(finalValues$metab[pointsPCA$id,])
      
      ids <- sprintf("The area you selected (phenotype average of %s), contains the following ids: <br/>%s", round(phenotypeAverage, digits = 2), paste(finalValues$currentIds, collapse = ', '))
      averagesFormatted <- mapply(function(x,y) paste(x, round(as.numeric(y), digits=4), sep=": "), names(finalValues$metab), averageSelected, SIMPLIFY=FALSE)
      average <- sprintf("The average values in this area are: <br/>%s", paste(averagesFormatted, collapse = '<br/>'))
      
      finalValues$tempInfo <- paste(ids, average, sep = '<br/>')
      
    } else if (!is.null(clusteringClick)&!equalsClick(temps$clusteringClick,clusteringClick)){
      switch(input$clusteringType, 
             SOM={
               if (length(clusteringData$SOM)==0)
                 return(" ")
               xValues <- clusteringData$SOM$grid$pts[,1]
               distX <- abs(xValues-clusteringClick$x)
               yValues <- clusteringData$SOM$grid$pts[,2]
               distY <- abs(yValues-clusteringClick$y)
               allBubbleIds <- 1:nrow(clusteringData$SOM$grid$pts)
               currentBubbleId <- allBubbleIds[distX==min(distX)&distY==min(distY)]
               if (length(currentBubbleId)==0)
                 return(" ")
               
               averageInBubble <- clusteringData$SOM$codes[[1]][currentBubbleId,]
               
               idsInBubble <- seq_len(nrow(clusteringData$SOM$data[[1]]))[clusteringData$SOM$unit.classif==currentBubbleId]
               finalValues$currentIds <- finalValues$idFromNum[idsInBubble]
               
               phenotypes <- finalValues$pheno[[as.numeric(input$selectedPhenotype)]]
               phenoInPoints <- phenotypes[idsInBubble]
               phenotypeAverage <- mean(phenoInPoints[!is.na(phenoInPoints)])
               ids <- sprintf("This node (%s, phenotype average of %s) contains the following ids: <br/>%s", currentBubbleId, round(phenotypeAverage, digits = 2), paste(finalValues$currentIds, collapse = ', '))
               averagesFormatted <- mapply(function(x,y) paste(x, round(as.numeric(y), digits=4), sep=": "), names(finalValues$metab), averageInBubble, SIMPLIFY=FALSE)
               average <- sprintf("The average values in this node are: <br/>%s", paste(averagesFormatted, collapse = '<br/>'))
               
               finalValues$tempInfo <- paste(ids, average, sep = '<br/>')
             },
             COSA={
             },
             DOC={  
             },
             {
               # default is using Clique
             }
      )
    } else if (!is.null(idsFromCluster)&!identical(idsFromCluster,temps$idsFromClusterOld)){
      finalValues$currentIds <- finalValues$idFromNum[as.character(idsFromCluster)]
      phenotypes <- finalValues$pheno[[as.numeric(input$selectedPhenotype)]]
      phenoInPoints <- phenotypes[idsFromCluster]
      phenotypeAverage <- mean(phenoInPoints[!is.na(phenoInPoints)])
      averageSelected <- colMeans(finalValues$metab[idsFromCluster,])
      #if(!is.null(input$clusterId)){
      
      ids <- sprintf("The cluster you selected (%s, phenotype average of %s), contains the following ids: <br/>%s", input$clusterId, round(phenotypeAverage, digits = 2), paste(finalValues$currentIds, collapse = ', '))
      averagesFormatted <- mapply(function(x,y) paste(x, round(as.numeric(y), digits=4), sep=": "), names(finalValues$metab), averageSelected, SIMPLIFY=FALSE)
      average <- sprintf("The average values in this area are: <br/>%s", paste(averagesFormatted, collapse = '<br/>'))
      
      finalValues$tempInfo <- paste(ids, average, sep = '<br/>')
    }
    temps$clusteringClick <- clusteringClick
    temps$clusteringBrush <- clusteringBrush
    temps$pcaBrush <- pcaBrush
    temps$idsFromClusterOld <- idsFromCluster
    
    return(HTML(finalValues$tempInfo))
    
  })
  
  output$downloadButton <- renderUI({
    if (length(finalValues$currentIds)==0){
      return(NULL)
    } else {
      downloadButton("downloadValues", "Download selected metabolites")
    }
  })
  
  # Downloads the current subset displayed in info
  output$downloadValues <- downloadHandler(
    filename = function() {
      paste("metabolites-", Sys.Date(), ".tsv", sep="")
    },
    content = function(file) {
      dataToDownload <- finalValues$completeData[finalValues$numFromId[as.character(finalValues$currentIds)],]
      write.table(dataToDownload, file, sep = "\t", dec = ",", row.names=FALSE)
    },
    contentType = "text/tab-separated-values"
  )
  
  PCASwitchObserver <- observe({
    if (input$clusteringType=="Nothing"){
      return()
    }
    if(input$pcaSwitch[[1]]){
      output$pcaAll <- renderUI({
        return({
            fluidRow(
              column(9,
                     plotOutput(outputId = "pca",
                                brush = "pca_brush", height = plotsSize, width = plotsSize)
              ),
              column(3,
                     selectInput(inputId = "clusterId",
                                 label = "Select the cluster:",
                                 choices = list()),
                     checkboxInput(inputId = "colorClusterNoPheno",
                                   label = "Color Cluster instead of phenotype"),
                     htmlOutput("metabUsed")
              )
            )
          })
        })
      current <- getCurrentClusters()
      len <- length(current)
      
      columns <- seq(0,len)
      if (input$clusteringType=="SOM"){
        columns <- sort(c(0,unique(clusteringData$SOM$unit.classif)))
      } 
      
      if (len>0){
        if (is.null(clusteringData$PhenoPValues[[input$clusteringType]][[as.numeric(input$selectedPhenotype)]])){
          clusteringData$PhenoPValues[[input$clusteringType]][[as.numeric(input$selectedPhenotype)]] <- unlist(sapply(current, function(x) calcPValueForCluster(x,finalValues$pheno[[as.numeric(input$selectedPhenotype)]])))
        }
        
        pValues <- clusteringData$PhenoPValues[[input$clusteringType]][[as.numeric(input$selectedPhenotype)]]
        pValues <- c(max(pValues)+1,pValues)
        
        if (length(columns)!=length(pValues)){
          print("Something is wrong")
          print(paste("Columns:",length(columns),"pValues:",length(pValues)))
        }
        columns <- columns[order(pValues)]
      }
      
      updateSelectInput(session, "clusterId", choices = columns)
    } else {
      output$pcaAll <- renderUI({
        plotOutput(outputId = "pca",
                   brush = "pca_brush", height = plotsSize, width = plotsSize)
      })
    }
  })
  
  getCurrentClusters <- function(){
    if (input$clusteringType=="Nothing"){
      return()
    }
    if (input$clusteringType=="Clique"){
      return(lapply(clusteringData$Clique, "[[", "objects"))
    } else if (input$clusteringType=="DOC"){
      return(getIdsDoc(clusteringData$DOC))
    } else if (input$clusteringType=="SOM"){
      ids <- 1:nrow(finalValues$metab)
      clusters <- list()
      clusterCounts <- sort(table(clusteringData$SOM$unit.classif),decreasing = TRUE)
      
      minicount<-1
      for (i in names(clusterCounts)){
        clusters[[minicount]] <- ids[clusteringData$SOM$unit.classif==as.numeric(i)]
        minicount<-minicount+1
      }
      return(clusters)
    }
    return(temps$index)
  }
  
  output$hist <- renderPlot({
    if (length(clusteringData$COSA)>0){
      if (is.null(clusteringData$COSA$hist)){
        clusteringData$COSA$hist <- hierclust(clusteringData$COSA$D, denplot = FALSE)
      }
      
      if (!is.null(input$clickHist)){
        if (is.null(clusteringData$COSA$hist)|equalsClick(input$clickHist,temps$clickHist)){
          return(clusteringData$COSA$hist)
        }
        rec.col = "blue"
        old.col = "blue"
        
        retval <- temps$index #list()
        oldk <- temps$oldk #NULL
        oldx <- temps$oldx #NULL
        id <- temps$id #1
        cat("Showing dynamic visualisation. Press Escape/Ctrl + C to stop.")
        
        x <- input$clickHist
        if (is.null(x)) {
          break
        }
        k <- min(which(rev(clusteringData$COSA$hist$height) < x$y))
        k <- max(k, 2)
        if (!is.null(oldx)) {
          modifiedRect(clusteringData$COSA$hist, k = oldk, x = oldx, border = old.col)
        }
        retval[[temps$id]] <- unlist(modifiedRect(clusteringData$COSA$hist, k = k, x = x$x, 
                                                 border = rec.col))
        temps$oldx <- x$x
        temps$oldk <- k
        
        grps <- temps$grps
        grps[retval[[temps$id]]] <- temps$id
        
        names(retval) <- paste("grp", 1:temps$id, sep = "")
        temps$id <- temps$id+1
        temps$clickHist <- x
        temps$grps <- grps
        temps$index <- unique(retval)
        temps$index <- temps$index[!sapply(temps$index, is.null)]
        #invisible(list(grps = grps, index = retval))
      }
      outPlot <- ggdendrogram(clusteringData$COSA$hist$dendro)
      if (nrow(temps$rects)>0){
        outPlot <- outPlot + geom_rect(data = temps$rects, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="blue")
      }
      outPlot
    }
  })
  
  modifiedRect <- function (tree, k = NULL, which = NULL, x = NULL, h = NULL, 
                            border = 2, cluster = NULL) {
    if (length(h) > 1L | length(k) > 1L) 
      stop("'k' and 'h' must be a scalar")
    if (!is.null(h)) {
      if (!is.null(k)) 
        stop("specify exactly one of 'k' and 'h'")
      k <- min(which(rev(tree$height) < h))
      k <- max(k, 2)
    }
    else if (is.null(k)) 
      stop("specify exactly one of 'k' and 'h'")
    if (k < 2 | k > length(tree$height)) 
      stop(gettextf("k must be between 2 and %d", length(tree$height)), 
           domain = NA)
    if (is.null(cluster)) 
      cluster <- cutree(tree, k = k)
    clustab <- table(cluster)[unique(cluster[tree$order])]
    m <- c(0, cumsum(clustab))
    if (!is.null(x)) {
      if (!is.null(which)) 
        stop("specify exactly one of 'which' and 'x'")
      which <- x
      for (n in seq_along(x)) which[n] <- max(which(m < x[n]))
    }
    else if (is.null(which)) 
      which <- 1L:k
    if (any(which > k)) 
      stop(gettextf("all elements of 'which' must be between 1 and %d", 
                    k), domain = NA)
    border <- rep_len(border, length(which))
    retval <- list()
    for (n in seq_along(which)) {
      newRow <- data.frame(x1=m[which[n]] + 0.66, y1=par("usr")[3L], x2=m[which[n] + 1] + 0.33, y2=mean(rev(tree$height)[(k - 1):k]), border=border[n])
      temps$rects <- rbind(temps$rects, newRow)
      retval[[n]] <- which(cluster == as.integer(names(clustab)[which[n]]))
    }
    temps$rects <- unique.data.frame(temps$rects)
    invisible(retval)
  }
  
  
}

#Scales vector to the range 0 to 1
range01 <- function(x){(x-min(x[!is.na(x)]))/(max(x[!is.na(x)])-min(x[!is.na(x)]))}
# round to k decimals
specify_decimal <- function(k) { function(x){trimws(format(round(x, k), nsmall=k))} }
# round to 4 decimals
fourDeci <- specify_decimal(4)



# Transforms a numbervector into colors
vecToCol <- function(data, colors){
  dataNoNa <- data[!is.na(data)]
  result <- data
  counter = 1
  for (i in sort(unique(dataNoNa))){
    result[result==i] <- colors[counter]
    counter <- counter + 1
  }
  result[is.na(result)]<-"black"
  result
}
#Equa
equalsBrush <- function(one,two){
  if (is.null(one)){
    if (is.null(two)){
      return(TRUE)
    }
    return(FALSE)
  } else if (is.null(two)) {
    return(FALSE)
  }
  return(one$xmin==two$xmin&one$xmax==two$xmax&one$ymin==two$ymin&one$ymax==two$ymax)
}
equalsClick <- function(one,two){
  if (is.null(one)){
    if (is.null(two)){
      return(TRUE)
    }
    return(FALSE)
  } else if (is.null(two)) {
    return(FALSE)
  }
  return(one$x==two$x&one$y==two$y)
}

#Calculates the w parameter for doc
calcW <- function(data,c=1){
  neighbors <- nn2(data,k=2)
  pairs <- neighbors$nn.idx[,2]
  counter <- 1
  wLocals <- lapply(pairs, function(x){
    one <- data[counter,]
    counter <- counter + 1
    two <- data[x,]
    wLocal <- sum(abs(one-two))/length(one)
    wLocal
  })
  w <- (sum(unlist(wLocals))/length(wLocals))*c
  w
}

calcPValueForCluster <- function(cluster,phenotype){
  phenoNoNa <- phenotype[!is.na(phenotype)]
  uniNoNa <- unique(phenoNoNa)
  if (length(uniNoNa)==2){
    return(getChiForCluster(cluster,phenotype,TRUE))
  } else {
    cluster <- as.numeric(cluster)
    if (class(phenotype)=="character"){
      phenotype <- factor(phenotype)
    }
    phenotype <- as.numeric(phenotype)
    clusterPhenos <- phenotype[cluster]
    clusterPhenos <- clusterPhenos[!is.na(clusterPhenos)]
    nonClusterPhenos <- phenotype[-cluster]
    nonClusterPhenos <- nonClusterPhenos[!is.na(nonClusterPhenos)]
    if(length(unique(c(clusterPhenos,nonClusterPhenos)))==1||length(clusterPhenos)==0||length(nonClusterPhenos)==0){
      return(1)
    } else if(length(unique(clusterPhenos))==1 & length(unique(nonClusterPhenos))==1){
      return(0)
    } else {
      return(t.test(clusterPhenos,nonClusterPhenos)$p.value)
    }
    
  }
}

#Calculate chi-squared value for cluster
getChiForCluster <- function(cluster,phenotype,pvalue=FALSE){
  phenoNoNa <- phenotype[!is.na(phenotype)]
  one <- unique(phenoNoNa)[[1]]
  two <- unique(phenoNoNa)[[2]]
  cluster <- as.numeric(unlist(cluster))
  pheno <- as.numeric(unlist(phenotype))
  phenoClust <- sum(pheno[cluster]==one,na.rm = T)
  phenoNonClust <- sum(pheno[-cluster]==one,na.rm = T)
  nonPhenoClust <- sum(pheno[cluster]==two,na.rm = T)
  nonPhenoNonClust <- sum(pheno[-cluster]==two,na.rm = T)
  #phenoNonClust <- length(pheno[pheno==one])-phenoClust
  #nonPhenoClust <- length(cluster)-phenoClust
  #nonPhenoNonClust <- length(pheno[pheno==two])-nonPhenoClust
  if (phenoClust+phenoNonClust+nonPhenoClust+nonPhenoNonClust!=length(phenoNoNa)){
    print("ERROR, table is not correct")
    return(NULL)
  }
  result <- data.frame(matrix(c(phenoClust,phenoNonClust,nonPhenoClust,nonPhenoNonClust),nrow = 2,ncol = 2))
  names(result) <- c("pheno","nonPheno")
  row.names(result) <- c("inCluster","outCluster")
  chi <- chisq.test(result)
  if (pvalue){
    return(chi$p.value)
  }
  return(chi$statistic)
}



# Create Shiny app ----
shinyApp(ui = ui, server = server)