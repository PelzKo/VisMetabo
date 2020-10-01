library(shiny)
library(shinyBS)
library(subspace)
library(ggplot2)
library(ggfortify)
library(rCOSA)
library(kohonen)
library(rJava)
#library(ggmap)
#library(ggdendro)

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
                accept = ".xlsx"),
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
                      choices = list("Clique","SOM", "COSA","DOC")),
          actionButton("resetClustering", "Reset the Clustering"),
          checkboxInput("pcaSwitch", "Color the Clusters in the PCA Plot"),
          uiOutput(outputId = "min"),
          uiOutput(outputId = "max"),
          checkboxInput(inputId = "reverse",
                        label = "Reverse Color Strength"),
          checkboxInput(inputId = "removeOutliers",
                        label = "Remove the highest and lowest 5% (Only of number phenotypes)"),
          #checkboxInput(inputId = "clusterRemove",
          #              label = "Remove them before the clustering"),
          textOutput(outputId = "phenoInfo"),
          plotOutput(outputId = "phenoHist")
          
          #selectInput(inputId = "coloredPurple",
          #            label = "Color in Purple:",
          #            choices = list("None","Cluster 0","Cluster 1"))
          
        ),
        
        # Main panel for displaying outputs ----
        mainPanel(
          
          # Output: Clustering & PCA ----
          #plotOutput(outputId = "clustering",
          #           click = "clustering_click"),
          bsModal("cosaHist", "Histogram of the rCosa", "go", size = "large",plotOutput("hist", click = "histClick")),
          uiOutput("clusteringPlot"),
          uiOutput("pcaAll"),
          uiOutput("downloadButton"),
          htmlOutput("info")
          
        )
      )
    )
  )
)


# Define server logic required ----
server <- function(input, output, session) {
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
  clusteringData <- reactiveValues(Clique = list(), CliqueMDS = list(), SOM = list(), COSA = list(), DOC = list(), DocMDS = list())
  # PCAData
  pcaValues <- reactiveValues()
  #Temp Values to distinguish which click/brush changed
  temps <- reactiveValues(clusteringClick = NULL,clusteringBrush = NULL, pcaBrush = NULL)
  
  
  output$pca <- renderPlot({
    if (length(data())==4){
      firstRunForClusteringMethod <- length(clusteringData[[input$clusteringType]])==0
      firstPCA <- length(pcaValues$visualisation)==0
      colorPalette <- redWhiteBlue(20)
      
      dev.interactive()
      inter <- dev.interactive()
      
      metabComplete <- data.frame(scale(finalValues$metab))
      
      if (firstRunForClusteringMethod){
        switch(input$clusteringType, 
               SOM={
                 clusteringData$SOM <- som(data.matrix(metabComplete))
               },
               COSA={
                 #clusteringData$COSA <- cosa2(metabComplete, niter = 7, noit = 15)
                 #clusteringData$COSA$smacof <- smacof(clusteringData$COSA$D, niter = 30, interc = 1, VERBOSE = FALSE, PLOT = FALSE)
                 clusteringData$COSA <- cosa2(metabComplete, niter = 1, noit = 1)
                 clusteringData$COSA$smacof <- smacof(clusteringData$COSA$D, niter = 1, interc = 1, VERBOSE = FALSE, PLOT = FALSE)
                 clusteringData$COSA$idsAndSmacof <- data.frame(cbind(names(finalValues$idFromNum),clusteringData$COSA$smacof$X))
                 names(clusteringData$COSA$idsAndSmacof) <- c("id","x","y")
                 
                 toggleModal(session, "cosaHist", toggle = "toggle")
                 output$hist <- renderPlot({ #AHHHHHH
                   hclst.cosa <- hierclust(clusteringData$COSA$D)
                   grps.cosa <- getclust(hclst.cosa)
                   #ggdendrogram(hclst.cosa)
                 })
                 #grps <- gggetClust(hclst.cosa)
                   
                 
               },
               DOC={
                 clusteringData$DOC <- runDoc(metabComplete)
                 
                 docMDSValues <- plotFromClusters(getIdsDoc(clusteringData$DOC), returnMDS = TRUE)
                 clusteringData$DocMDS <- data.frame(cbind(seq_len(nrow(docMDSValues)),docMDSValues))
                 names(clusteringData$DocMDS) <- c("id","x","y")
               },
               {
                 # default is using Clique
                 clusteringData$Clique <- CLIQUE(metabComplete) 
                 
                 cliqueMDSValues <- plotFromClusters(lapply(clusteringData$Clique, `[[`, "objects"), returnMDS = TRUE)
                 clusteringData$CliqueMDS <- data.frame(cbind(seq_len(nrow(cliqueMDSValues)),cliqueMDSValues))
                 names(clusteringData$CliqueMDS) <- c("id","x","y")
                 
               }
        )
      }
      
      phenotype <- finalValues$pheno[[as.numeric(input$selectedPhenotype)]]
      
      
      if (class(phenotype)=="character"){
        phenotype <- as.numeric(factor(phenotype))
      } else {
        if (input$removeOutliers[[1]]){
          n <- 5
          lowestTooHigh <- min(phenotype[phenotype >= quantile(phenotype,prob=1-n/100)])
          highestTooLow <- max(phenotype[phenotype <= quantile(phenotype,prob=n/100)])
          phenotype[phenotype>lowestTooHigh] <- lowestTooHigh
          phenotype[phenotype<highestTooLow] <- highestTooLow
        }
      }
      
      if (length(unique(finalValues$pheno[[as.numeric(input$selectedPhenotype)]]))==1&&finalValues$pheno[[as.numeric(input$selectedPhenotype)]][[1]]==0){
        coloredPoints <- rgb(0, 0, 0, maxColorValue=255, alpha=255)
        
        
        output$phenoInfo <- NULL
        output$phenoHist <- NULL
        
      } else {
        if (input$reverse[[1]]){
          output$min <- renderUI(HTML(sprintf(paste("<div style='background-color: ",colorPalette[1],";height: 18px;width: 18px;float: left;margin-right: 3px;'></div> Min: %s",sep = ""), fourDeci(min(phenotype)))))
          output$max <- renderUI(HTML(sprintf(paste("<div style='background-color: ",colorPalette[length(colorPalette)],";height: 18px;width: 18px;float: left;margin-right: 3px;'></div> Max: %s",sep = ""), fourDeci(max(phenotype)))))
          coloredPoints <- vecToCol(phenotype,colorPalette)
        } else {
          output$min <- renderUI(HTML(sprintf(paste("<div style='background-color: ",colorPalette[length(colorPalette)],";height: 18px;width: 18px;float: left;margin-right: 3px;'></div> Min: %s",sep = ""), fourDeci(min(phenotype)))))
          output$max <- renderUI(HTML(sprintf(paste("<div style='background-color: ",colorPalette[1],";height: 18px;width: 18px;float: left;margin-right: 3px;'></div> Max: %s",sep = ""), fourDeci(max(phenotype)))))
          coloredPoints <- vecToCol(phenotype,rev(colorPalette))
        }
        output$phenoInfo <- renderText(sprintf("Mean of phenotypes after normalization between 0 and 1: %s", fourDeci(mean(range01(phenotype)))))
        output$phenoHist <- renderPlot(hist(phenotype))
      }
      
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
               
               clusters <- as.character(sort(unique(clusteringData$SOM$unit.classif)))
               colorClusters <- rainbow(length(clusters))
               names(colorClusters) <- clusters
               clusterColors <- clusteringData$SOM$unit.classif
               for (num in clusters){
                 clusterColors <- replace(clusterColors, clusterColors==as.numeric(num), colorClusters[[num]])
               }
               
             },
             COSA={
               output$clusteringPlot <- renderUI(plotOutput(outputId = "clustering", brush = "clustering_brush"))
               output$clustering <- renderPlot(
                 plotSmacof(clusteringData$COSA$smacof[["X"]], cols = coloredPoints))
             },
             DOC={
               output$clusteringPlot <- renderUI(plotOutput(outputId = "clustering", brush = "clustering_brush"))
               
               idsInClustersDoc <- getIdsDoc(clusteringData$DOC)
               output$clustering <- renderPlot(plotFromClusters(idsInClustersDoc, label = FALSE, colors = coloredPoints))
               
               
               clusterColors <- rep("black",max(unlist(idsInClustersDoc)))
               if(!is.null(input$clusterId)&input$pcaSwitch[[1]]){
                 currentCluster <- idsInClustersDoc[[as.numeric(input$clusterId)]]
                 usedDimensions <- getDimsDoc(clusteringData$DOC)[as.numeric(input$clusterId),]
                 avgs <- getAvgsDoc(clusteringData$DOC)[as.numeric(input$clusterId),]
                 metabs <- names(metabComplete)[usedDimensions]
                 usedMetabs <- sprintf("This cluster was calculed using the following metabolites: <br/>%s", paste(metabs, collapse = '<br/>'))
                 averagesFormatted <- mapply(function(x,y) paste(x, round(as.numeric(y), digits=4), sep=": "), metabs, avgs, SIMPLIFY=FALSE)
                 metabAvgs <- sprintf("<br/><br/>The metabolites have the following averages: <br/>%s", paste(averagesFormatted, collapse = '<br/>'))
                 output$metabUsed <- renderUI(HTML(paste(usedMetabs,metabAvgs)))
                 
                 clusterColors[currentCluster] <- "red"
               } 
             },
             {
               output$clusteringPlot <- renderUI(plotOutput(outputId = "clustering", brush = "clustering_brush"))
               # default is using Clique
               
               output$clustering <- renderPlot(plotFromClusters(lapply(clusteringData$Clique, `[[`, "objects"), label = FALSE, colors = coloredPoints))
               
               
               clusterColors <- rep("black",max(unlist(clusteringData$Clique)))
               if(!is.null(input$clusterId)){
                 currentCluster <- clusteringData$Clique[as.numeric(input$clusterId)][[1]]
                 metabs <- names(metabComplete)[currentCluster$subspace]
                 usedMetabs <- sprintf("This cluster was calculed using the following metabolites: <br/>%s", paste(metabs, collapse = '<br/>'))
                 output$metabUsed <- renderUI(HTML(usedMetabs))
                 
                 clusterColors[currentCluster$objects] <- "red"
               }
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
      if (input$pcaSwitch[[1]]){
        coloredPointsPCA <- clusterColors
      }
      
      plotNoLims(pcaValues$visualisation, cols = coloredPointsPCA)
    }
  })
  
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
                 brush = "pca_brush")
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
      clusteringData$Clique = list()
      clusteringData$SOM = list()
      clusteringData$COSA = list()
      clusteringData$DOC = list()
  })
      
  # Goto Output Panel
  observeEvent(input$toOut, {
    if (length(data())==4){
      clusteringData$Clique = list()
      clusteringData$SOM = list()
      clusteringData$COSA = list()
      clusteringData$DOC = list()
      inputData <- data()$values[complete.cases(data()$values), ]
      
      tempId <- as.numeric(input$idField)
      tempStart <- as.numeric(input$metabStart)
      tempEnd <- as.numeric(input$metabEnd)
      
      if (length(inputData[[tempId]])!=length(unique(inputData[[tempId]]))){
        output$idError <- renderText("Identifier is not unique")
        return()
      }
      
      if (input$externalPheno[[1]]&&!is.null(input$inputPheno)){
        dataTemp <- merge(inputData, readPhenoFile(input$inputPheno$datapath), by = names(inputData)[[tempId]])
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
                       accept = ".xlsx"))
    } else {
      return(NULL)
    }
  })
  
  
  output$info <- renderUI({
    clusteringClick <- input$clustering_click
    clusteringBrush <- input$clustering_brush
    pcaBrush <- input$pca_brush
    
    if (!is.null(clusteringBrush)&!equalsBrush(temps$clusteringBrush,clusteringBrush)){
      switch(input$clusteringType, 
             SOM={
             },
             COSA={
               points <- brushedPoints(clusteringData$COSA$idsAndSmacof, clusteringBrush, xvar = "x", yvar = "y")
               finalValues$currentIds <- finalValues$idFromNum[as.character(points$id)]
               phenotypes <- finalValues$pheno[[as.numeric(input$selectedPhenotype)]]
               phenotypeAverage <- mean(phenotypes[points$id])
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
               phenotypeAverage <- mean(phenotypes[points$id])
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
               phenotypeAverage <- mean(phenotypes[points$id])
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
      phenotypeAverage <- mean(phenotypes[pointsPCA$id])
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
               yValues <- clusteringData$SOM$grid$pts[,2]
               allBubbleIds <- 1:nrow(clusteringData$SOM$grid$pts)
               currentBubbleId <- allBubbleIds[xValues==round(as.numeric(clusteringClick$x))&yValues==round(as.numeric(clusteringClick$y))]
               if (length(currentBubbleId)==0)
                 return(" ")
               
               averageInBubble <- clusteringData$SOM$codes[[1]][currentBubbleId,]
               
               idsInBubble <- seq_len(nrow(clusteringData$SOM$data[[1]]))[clusteringData$SOM$unit.classif==currentBubbleId]
               finalValues$currentIds <- finalValues$idFromNum[idsInBubble]
               
               phenotypes <- finalValues$pheno[[as.numeric(input$selectedPhenotype)]]
               phenotypeAverage <- mean(phenotypes[idsInBubble])
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
    }
    temps$clusteringClick <- clusteringClick
    temps$clusteringBrush <- clusteringBrush
    temps$pcaBrush <- pcaBrush
    
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
    if(input$pcaSwitch[[1]]&(input$clusteringType=="Clique"|input$clusteringType=="DOC")){
      output$pcaAll <- renderUI({
        return({
            fluidRow(
              column(8,
                     plotOutput(outputId = "pca",
                                brush = "pca_brush")
              ),
              column(4,
                     selectInput(inputId = "clusterId",
                                 label = "Select the cluster to be colored:",
                                 choices = list()),
                     htmlOutput("metabUsed")
              )
            )
          })
        })
      if (input$clusteringType=="Clique"){
        len <- length(clusteringData$Clique)
      } else {
        len <- length(getIdsDoc(clusteringData$DOC))
      }
      columns <- seq_len(len)
      updateSelectInput(session, "clusterId", choices = columns)
    } else {
      output$pcaAll <- renderUI({
        plotOutput(outputId = "pca",
                   brush = "pca_brush")
      })
    }
  })
  
}

#Scales vector to the range 0 to 1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# round to k decimals
specify_decimal <- function(k) { function(x){trimws(format(round(x, k), nsmall=k))} }
# round to 4 decimals
fourDeci <- specify_decimal(4)



# Transforms a numbervector into colors
vecToCol <- function(data, colors){
  result <- data
  result[] <- colors[length(colors)]
  range = max(data)-min(data)
  intervals = range/(length(colors)-1)
  lastValue = -Inf
  counter = 1
  for (i in seq(min(data), max(data), by=intervals)){
    result[data>lastValue&data<=i] <- colors[counter]
    lastValue <- i
    counter <- counter + 1
  }
  result
}

gggetClust <- function (hc, ngr = 20, rec.col = "blue", old.col = "blue") 
{
  retval <- list()
  oldk <- NULL
  oldx <- NULL
  cat("Showing dynamic visualisation. Press Escape/Ctrl + C to stop.")
  for (n in 1:ngr) {
    x <- gglocator(1)
    if (is.null(x)) {
      break
    }
    k <- min(which(rev(hc$height) < x$y))
    k <- max(k, 2)
    if (!is.null(oldx)) {
      rect.hclust(hc, k = oldk, x = oldx, border = old.col)
    }
    retval[[n]] <- unlist(rect.hclust(hc, k = k, x = x$x, 
                                      border = rec.col))
    oldx <- x$x
    oldk <- k
    ngr <- n
  }
  grps <- rep(0, length(hc$order))
  for (n in 1:ngr) {
    grps[retval[[n]]] <- n
  }
  names(retval) <- paste("grp", 1:n, sep = "")
  invisible(list(grps = grps, index = retval))
}

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


# Create Shiny app ----
shinyApp(ui = ui, server = server)