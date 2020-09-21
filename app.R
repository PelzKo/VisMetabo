library(shiny)
library(subspace)
library(ggplot2)
library(ggfortify)
library(rCOSA)
library(kohonen)

source("ReadingData.R")
source("plotHeatMap.R")
source("smacofFixed.R")

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
          uiOutput("clusteringPlot"),
          plotOutput(outputId = "pca",
                     brush = "pca_brush"),
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
  clusteringData <- reactiveValues(Clique = list(), SOM = list(), COSA = list(), DOC = list())
  # PCAData
  pcaValues <- reactiveValues()
  
  
  output$pca <- renderPlot({
    if (length(data())==4){
      firstRunForClusteringMethod <- length(clusteringData[[input$clusteringType]])==0
      firstPCA <- length(pcaValues$visualisation)==0
      colorPalette <-colorRampPalette(c("red","white","blue"), space="Lab")(20)
      
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
                 clusteringData$COSA$idsAndSmacof <- data.frame(cbind(finalValues$id,clusteringData$COSA$smacof$X))
                 names(clusteringData$COSA$idsAndSmacof) <- c("id","x","y")
               },
               DOC={
                 #clusteringData$DOC <- doc(metabComplete)   
               },
               {
                 # default is using Clique
                 clusteringData$Clique <- CLIQUE(metabComplete) 
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
             },
             COSA={
               output$clusteringPlot <- renderUI(plotOutput(outputId = "clustering", brush = "clustering_brush"))
               output$clustering <- renderPlot(
                 plotSmacof(clusteringData$COSA$smacof[["X"]], cols = coloredPoints))
             },
             DOC={
               output$clusteringPlot <- renderUI(plotOutput(outputId = "clustering", brush = "clustering_brush"))
               #clusteringData$DOC <- doc(metabComplete)   
             },
             {
               output$clusteringPlot <- renderUI(plotOutput(outputId = "clustering", brush = "clustering_brush"))
               # default is using Clique
               #output$clustering <- renderPlot(plot(clusteringData$Clique,metabComplete))
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
        pcaValues$visWithId <- data.frame(cbind(finalValues$id,pcaValues$visualisation))
        names(pcaValues$visWithId) <- c("id","x","y")
        
      }
      
      #output$clustering <- renderPlot({
      #  ggplot(clusteringData[[input$clusteringType]])
      #})
      
      
      #lim <- c(min(clusterInfo$visualisation),max(clusterInfo$visualisation))
      #coordinates <- data.frame(x=clusterInfo$visualisation[,1],y=clusterInfo$visualisation[,2]) #Does not work
      #plot(clusterInfo$visualisation[, 1:2],col=clusters,bg = clusters,pch = 21,xlim = lim,ylim = lim)
      #autoplot(clusterInfo$visualisation,col=clusters, shape = 16)
      #finalValues$metabForBrush <- cbind(metabComplete,pcaValues$visualisation)
      plotNoLims(pcaValues$visualisation, cols = coloredPoints)
      #ggplot(pcaValues$visualisation, aes(PC1,PC2, color = coloredPoints))+
      #  geom_point()+
      #  theme(legend.position="none")+
      #  labs(x=paste0("PC1 (",pcaValues$percentage[1],")"), y=paste0("PC2 (",pcaValues$percentage[2],")")) 
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
      finalValues$realIds <- seq_len(length(inputData[[tempId]]))
      names(finalValues$realIds) <- inputData[[tempId]]
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
    
    if (!is.null(clusteringBrush)){
      switch(input$clusteringType, 
             SOM={
             },
             COSA={
               points <- brushedPoints(clusteringData$COSA$idsAndSmacof, clusteringBrush, xvar = "x", yvar = "y")
               finalValues$currentIds <- finalValues$realIds[as.character(points$id)]
               phenotypes <- finalValues$pheno[[as.numeric(input$selectedPhenotype)]]
               phenotypeAverage <- mean(phenotypes[finalValues$currentIds])
               averageSelected <- colMeans(finalValues$metab[finalValues$currentIds,])
               
               ids <- sprintf("The area you selected (phenotype average of %s), contains the following ids: <br/>%s", round(phenotypeAverage, digits = 2), paste(points$id, collapse = ', '))
               averagesFormatted <- mapply(function(x,y) paste(x, round(as.numeric(y), digits=4), sep=": "), names(finalValues$metab), averageSelected, SIMPLIFY=FALSE)
               average <- sprintf("The average values in this area are: <br/>%s", paste(averagesFormatted, collapse = '<br/>'))
               
               finalValues$tempInfo <- paste(ids, average, sep = '<br/>')
             },
             DOC={
                 
             },
             {
               # default is using Clique
               #output$clustering <- renderPlot(plot(clusteringData$Clique,metabComplete))
             }
      )
      
    } else if (!is.null(pcaBrush)){
      pointsPCA <- brushedPoints(pcaValues$visWithId, pcaBrush, xvar = "x", yvar = "y")
      finalValues$currentIds <- finalValues$realIds[as.character(pointsPCA$id)]
      phenotypes <- finalValues$pheno[[as.numeric(input$selectedPhenotype)]]
      phenotypeAverage <- mean(phenotypes[finalValues$currentIds])
      averageSelected <- colMeans(finalValues$metab[finalValues$currentIds,])
      
      ids <- sprintf("The area you selected (phenotype average of %s), contains the following ids: <br/>%s", round(phenotypeAverage, digits = 2), paste(pointsPCA$id, collapse = ', '))
      averagesFormatted <- mapply(function(x,y) paste(x, round(as.numeric(y), digits=4), sep=": "), names(finalValues$metab), averageSelected, SIMPLIFY=FALSE)
      average <- sprintf("The average values in this area are: <br/>%s", paste(averagesFormatted, collapse = '<br/>'))
      
      finalValues$tempInfo <- paste(ids, average, sep = '<br/>')
      
    } else if (!is.null(clusteringClick)){
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
               finalValues$currentIds <- finalValues$realIds[as.character(idsInBubble)]
               
               phenotypes <- finalValues$pheno[[as.numeric(input$selectedPhenotype)]]
               phenotypeAverage <- mean(phenotypes[finalValues$currentIds])
               ids <- sprintf("This node (%s, phenotype average of %s) contains the following ids: <br/>%s", currentBubbleId, round(phenotypeAverage, digits = 2), paste(idsInBubble, collapse = ', '))
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
      dataToDownload <- finalValues$completeData[finalValues$currentIds,]
      write.table(dataToDownload, file, sep = "\t", dec = ",", row.names=FALSE)
    },
    contentType = "text/tab-separated-values"
  )
  
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

doc <- function(data, alpha, beta, w){
  arr <- java_object_from_data(data)
  #Now that the data is in the correct format, we can call into our Java Code that will then call into the
  #actual implementation of the Algorithm
  res <- rJava::.jcall("ClusteringApplier",returnSig="[Li9/subspace/base/Cluster;",method="fires",arr,
                       base_dbscan_epsilon,
                       as.integer(base_dbscan_minpts),
                       minimumpercent,
                       as.integer(k),
                       as.integer(mu),
                       as.integer(minclu),
                       split,
                       post_dbscan_epsilon,
                       as.integer(post_dbscan_minpts),
                       evalArray=F)
  #We can then turn the Java Clustering Object that was returned into an R-Friendly S3-Object
  res <- r_clusters_from_java_clusters(res)
  return(res)
  
}


# Create Shiny app ----
shinyApp(ui = ui, server = server)