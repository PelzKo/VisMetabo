library(shiny)
library(subspace)
library(ape)

source("C:/Users/Konstantin/Desktop/Uni/6Semester/BachelorArbeit/Code/ReadingData.R")

options(shiny.maxRequestSize=200*1024^2)
options(java.parameters = "-Xmx16000m")


# Define UI for app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Analysis of Metabolite data"),
  
  tabsetPanel(id = "tabs", type = "tabs",
    tabPanel(title = "Data Input", value = "0",
      titlePanel("Please upload your data"),
      "Important things about the data:",
      "- Rows are individuals, Columns are features",
      fileInput(inputId = "inputData",
                label = "Data:",
                accept = ".xlsx"),
      actionButton("toTypes", "Upload")
             
    ),
    
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
      tableOutput(outputId = "dataAnalysis"),
      actionButton("toOut", "Cluster Now")
    ),
    
    tabPanel(title = "Output", value = "2",
      # Sidebar layout with input and output definitions ----
      sidebarLayout(
        
        # Sidebar panel for inputs ----
        sidebarPanel(
          
          # Input: Slider for the number of bins ----
          selectInput(inputId = "coloredRed",
                      label = "Color in Red:",
                      choices = list()),
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
          plotOutput(outputId = "distPlot")
          
        )
      )
    )
  )
)


# Define server logic required ----
server <- function(input, output, session) {
  data <- reactive({
    #validate(
    #  need(!is.null(input$inputData),"Please upload a file first")
    #)
    #readFile(input$inputData$datapath)
    if (!is.null(input$inputData)){
      readFile(input$inputData$datapath)
    }
    else{
      list("Error" = "Please upload a file first!")
    }
  })
  
  finalValues <- reactiveValues(id = 1, metab = list(), pheno = list())
  clusterInfo <- reactiveValues(cluster = list(), res = list())
  
  output$distPlot <- renderPlot({
    if (length(data())==4){
      firstRun <- length(clusterInfo$cluster)==0
      
      #metabComplete <- finalValues$metab[complete.cases(finalValues$metab), ]
      metabComplete <- scale(finalValues$metab)
      
      if (firstRun){
        clusterInfo$cluster <- SubClu(metabComplete)
        for (i in seq_len(length(clusterInfo$cluster))){
          metabComplete[clusterInfo$cluster[[i]][[2]],clusterInfo$cluster[[i]][[1]]]<-10#metabComplete[clusterInfo$cluster[[i]][[2]],clusterInfo$cluster[[i]][[1]]]*100
          #clusters[clusterInfo$cluster[[i]][[2]]]=colors[[i]]
        }
      }
      
      #clusters <- rep("grey",length(metabComplete[,1]))
      #colors <- c("red","blue","green","purple")
      phenotype <- finalValues$pheno[[as.numeric(input$coloredRed)]]
      
      
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
      
      if (length(unique(finalValues$pheno[[as.numeric(input$coloredRed)]]))==1&&finalValues$pheno[[as.numeric(input$coloredRed)]][[1]]==0){
        clusters <- rgb(0, 0, 0, maxColorValue=255, alpha=255)
        
        
        output$phenoInfo <- NULL
        output$phenoHist <- NULL
        
      } else {
        if (input$reverse[[1]]){
          output$min <- renderUI(HTML(sprintf("<div style='background-color: #FF0000FF;height: 18px;width: 18px;float: left;margin-right: 3px;'></div> Min: %s", fourDeci(min(phenotype)))))
          output$max <- renderUI(HTML(sprintf("<div style='background-color: #000000FF;height: 18px;width: 18px;float: left;margin-right: 3px;'></div> Max: %s", fourDeci(max(phenotype)))))
          clusters <- rgb((1-range01(phenotype))*255, 0, 0, maxColorValue=255, alpha=255)
        } else {
          output$min <- renderUI(HTML(sprintf("<div style='background-color: #000000FF;height: 18px;width: 18px;float: left;margin-right: 3px;'></div> Min: %s", fourDeci(min(phenotype)))))
          output$max <- renderUI(HTML(sprintf("<div style='background-color: #FF0000FF;height: 18px;width: 18px;float: left;margin-right: 3px;'></div> Max: %s", fourDeci(max(phenotype)))))
          clusters <- rgb(range01(phenotype)*255, 0, 0, maxColorValue=255, alpha=255)
        }
        
        
        
        output$phenoInfo <- renderText(sprintf("Mean of phenotypes after normalization between 0 and 1: %s", fourDeci(mean(range01(phenotype)))))
        output$phenoHist <- renderPlot(hist(phenotype))
        
        
      }
      
      
      if (firstRun){
        test <- dist(metabComplete)
        clusterInfo$res <- pcoa(test)
      }
      
      
      lim <- c(min(clusterInfo$res$vectors),max(clusterInfo$res$vectors))
      plot(clusterInfo$res$vectors,col=clusters,bg = clusters,pch = 21,xlim = lim,ylim = lim)
    }
  })
  
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
  
  observeEvent(input$toOut, {
    if (length(data())==4){
      clusterInfo$cluster = list()
      dataNoNa <- data()$values[complete.cases(data()$values), ]
      
      
      tempId <- as.numeric(input$idField)
      tempStart <- as.numeric(input$metabStart)
      tempEnd <- as.numeric(input$metabEnd)
      #validate(
      #  need(length(dataNoNa[[tempId]])==length(unique(dataNoNa[[tempId]])),message = "Identifier is not unique")
      #)
      if (length(dataNoNa[[tempId]])!=length(unique(dataNoNa[[tempId]]))){
        output$idError <- renderText("Identifier is not unique")
        return()
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
        updateSelectInput(session, "coloredRed", choices = columns, selected = columns[["None"]])
      }
      
      
    }
    updateTabsetPanel(session, "tabs",
                      selected = "2")
  })
  
  output$dataAnalysis <- renderTable({
    if (length(data())==4){
      data()[c(2,3,4)]
    }
    else{
      data()
    }
    
  })
  
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
specify_decimal <- function(k) { function(x){trimws(format(round(x, k), nsmall=k))} }
fourDeci <- specify_decimal(4)


# Create Shiny app ----
shinyApp(ui = ui, server = server)