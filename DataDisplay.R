library(shiny)
library(ape)

#testMatrix <- matrix(c(3,3,1,2.1,2,2,1,1,3),3,3,dimnames = list(c("S1","S2","S3"),c("Gen1","Gen2","Gen3")))
#test <- dist(testMatrix)
test <- dist(metabTest)

res <- pcoa(test)
#res2 <- pcoa(mite.D)
#clusters <- c(0,0,1)
#clusters[clusters==0]="red"
#clusters[clusters==1]="blue"
#plot(res$vectors,col=clusters,bg = clusters,pch = 21)

clusters <- rep("black",length(metabTest[,1]))
colors <- c("red","blue","green","purple")
for (i in seq_len(length(clusterInfo))){
  metabTest[clusterInfo[[i]][[2]],clusterInfo[[i]][[1]]]<-metabTest[clusterInfo[[i]][[2]],clusterInfo[[i]][[1]]]*5
  clusters[clusterInfo[[i]][[2]]]=colors[[i]]
}

plot(res$vectors,col=clusters,bg = clusters,pch = 21)

'
# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Analysis of Metabolite data"),
  
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
      
      # Sidebar panel for inputs ----
      sidebarPanel(
        
        # Input: Slider for the number of bins ----
        selectInput(inputId = "coloredRed",
                    label = "Color in red:",
                    list("None","Cluster 0","Cluster 1")),
        
        selectInput(inputId = "coloredBlue",
                    label = "Color in Blue:",
                    list("None","Cluster 0","Cluster 1"))
        
      ),
      
      # Main panel for displaying outputs ----
      mainPanel(
        
        # Output: Histogram ----
        plotOutput(outputId = "distPlot")
        
      )
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  output$distPlot <- renderPlot({
    
    #x    <- faithful$waiting
    #bins <- seq(min(x), max(x), length.out = input$bins + 1)
    #
    #hist(x, breaks = bins, col = "#75AADB", border = "white",
    #     xlab = "Waiting time to next eruption (in mins)",
    #     main = "Histogram of waiting times")
    testMatrix <- matrix(c(3,3,1,2.1,2,2,1,1,3),3,3,dimnames = list(c("S1","S2","S3"),c("Gen1","Gen2","Gen3")))
    test <- dist(testMatrix)
    res <- pcoa(test)
    
    clusters <- c(0,0,1)
    colors <- c("black","black","black")
    
    clusterId <- strsplit(input$coloredRed," ")[[1]]
    
    if (length(clusterId)>1){
      colors[clusters==clusterId[[2]]]="red"
    }
    
    
    clusterId <- strsplit(input$coloredBlue," ")[[1]]
    
    if (length(clusterId)>1){
      colors[clusters==clusterId[[2]]]="blue"
    }
    
    
    lim <- c(min(res$vectors),max(res$vectors))
    plot(res$vectors,col=colors,bg = colors,pch = 21,xlim = lim,ylim = lim)
    
  })
  
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)
'
