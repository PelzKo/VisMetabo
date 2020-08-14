library(shiny)
library(subspace)
library(ggplot2)
library(ggfortify)
library(rCOSA)
library(kohonen)

source("ReadingData.R")

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
          
          # Output: Histogram ----
          plotOutput(outputId = "distPlot",
                     brush = "plot_brush"),
          verbatimTextOutput("info")
          
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
  
  output$distPlot <- renderPlot({
    if (length(data())==4){
      firstRunForClusteringMethod <- length(clusteringData$as.character(input$selectedPhenotype))==0
      colorPalette <-colorRampPalette(c("red","white","blue"), space="Lab")(20)
      
      #metabComplete <- finalValues$metab[complete.cases(finalValues$metab), ]
      metabComplete <- scale(finalValues$metab)
      
      if (firstRunForClusteringMethod){
        clusterInfo$cluster <- CLIQUE(metabComplete)
        for (i in seq_len(length(clusterInfo$cluster))){
          metabComplete[clusterInfo$cluster[[i]][[2]],clusterInfo$cluster[[i]][[1]]]<-10#metabComplete[clusterInfo$cluster[[i]][[2]],clusterInfo$cluster[[i]][[1]]]*100
          #clusters[clusterInfo$cluster[[i]][[2]]]=colors[[i]]
        }
      }
      
      #clusters <- rep("grey",length(metabComplete[,1]))
      #colors <- c("red","blue","green","purple")
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
        clusters <- rgb(0, 0, 0, maxColorValue=255, alpha=255)
        
        
        output$phenoInfo <- NULL
        output$phenoHist <- NULL
        
      } else {
        if (input$reverse[[1]]){
          output$min <- renderUI(HTML(sprintf(paste("<div style='background-color: ",colorPalette[1],";height: 18px;width: 18px;float: left;margin-right: 3px;'></div> Min: %s",sep = ""), fourDeci(min(phenotype)))))
          output$max <- renderUI(HTML(sprintf(paste("<div style='background-color: ",colorPalette[length(colorPalette)],";height: 18px;width: 18px;float: left;margin-right: 3px;'></div> Max: %s",sep = ""), fourDeci(max(phenotype)))))
          clusters <- vecToCol(phenotype,colorPalette)
        } else {
          output$min <- renderUI(HTML(sprintf(paste("<div style='background-color: ",colorPalette[length(colorPalette)],";height: 18px;width: 18px;float: left;margin-right: 3px;'></div> Min: %s",sep = ""), fourDeci(min(phenotype)))))
          output$max <- renderUI(HTML(sprintf(paste("<div style='background-color: ",colorPalette[1],";height: 18px;width: 18px;float: left;margin-right: 3px;'></div> Max: %s",sep = ""), fourDeci(max(phenotype)))))
          clusters <- vecToCol(phenotype,rev(colorPalette))
        }
        
        
        
        output$phenoInfo <- renderText(sprintf("Mean of phenotypes after normalization between 0 and 1: %s", fourDeci(mean(range01(phenotype)))))
        output$phenoHist <- renderPlot(hist(phenotype))
        
        
      }
      
      
      if (firstRunForClusteringMethod){
        #test <- dist(metabComplete)
        #clusterInfo$visualisation <- pcoa(test)
        #clusterInfo$visualisation <- umap(metabComplete)
        pca_data <- prcomp(metabComplete, scale. = TRUE)
        ## Let us calculat the variances covered by components.
        pca_data_perc=round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)
        clusterInfo$percentage <- pca_data_perc
        
        ## create a data frame with principal component 1 (PC1), PC2, Conditions and sample names
        df_pca_data = data.frame(PC1 = pca_data$x[,1], PC2 = pca_data$x[,2])
        clusterInfo$visualisation <- df_pca_data
        
      }
      
      
      #lim <- c(min(clusterInfo$visualisation),max(clusterInfo$visualisation))
      #coordinates <- data.frame(x=clusterInfo$visualisation[,1],y=clusterInfo$visualisation[,2]) #Does not work
      #plot(clusterInfo$visualisation[, 1:2],col=clusters,bg = clusters,pch = 21,xlim = lim,ylim = lim)
      #autoplot(clusterInfo$visualisation,col=clusters, shape = 16)
      finalValues$metabForBrush <- cbind(metabComplete,clusterInfo$visualisation)
      ggplot(clusterInfo$visualisation, aes(PC1,PC2, color = clusters))+
        geom_point()+
        labs(x=paste0("PC1 (",clusterInfo$percentage[1],")"), y=paste0("PC2 (",clusterInfo$percentage[2],")")) 
    }
  })
  
  # Info about the selected points
  output$info <- renderPrint({
    # With ggplot2, no need to tell it what the x and y variables are.
    # threshold: set max distance, in pixels
    # maxpoints: maximum number of rows to return
    # addDist: add column with distance, in pixels
    #nearPoints(clusterInfo$visualisation$vectors, input$plot_hover, threshold = 10, maxpoints = 1,
               #addDist = FALSE)
    brushedPoints(finalValues$metabForBrush, input$plot_brush)
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
      dataNoNa <- data()$values[complete.cases(data()$values), ]
      
      tempId <- as.numeric(input$idField)
      tempStart <- as.numeric(input$metabStart)
      tempEnd <- as.numeric(input$metabEnd)
      
      if (length(dataNoNa[[tempId]])!=length(unique(dataNoNa[[tempId]]))){
        output$idError <- renderText("Identifier is not unique")
        return()
      }
      
      if (input$externalPheno[[1]]&&!is.null(input$inputPheno)){
        dataNoNaTemp <- merge(dataNoNa, readPhenoFile(input$inputPheno$datapath), by = names(dataNoNa)[[tempId]])
        if (nrow(dataNoNaTemp)==0){
          showNotification("Could not match any metab ids with phenotype ids. Please check you have the correct 
                           ID column, they are named the same and the values are written in the same format.",type = "warning")
          showNotification("Continuing with the phenotypes from the metabolite data sheet.",type = "message")
        } else {
          dataNoNa <- dataNoNaTemp
        }
      }
      
      
      finalValues$id <- dataNoNa[[tempId]]
      finalValues$metab <- dataNoNa[c(tempStart:tempEnd)]
      #none <- numeric(length(finalValues$id[[1]]))
      #names(none) <- "None"
      finalValues$pheno <- dataNoNa[-c(tempId,c(tempStart:tempEnd))]
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


# Create Shiny app ----
shinyApp(ui = ui, server = server)