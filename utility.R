getColors <- function(n) {
  if (n>9){
    return (colorRampPalette(c("red","white","blue"), space="Lab")(n))
  } else if (n==2){
    return(c("red","blue"))
  }
  return (rainbow(n))
}

# Plot SOM Heatmap for variable 
#
# Shane Lynn 13/1/2014

# this function is to plot the heatmap of a variable

plotHeatMap <- function(som_model, phenotype_data, variable=0){    
  # Plot a heatmap for any variable from the data set "phenotype_data".
  # If variable is 0, an interactive window will be provided to choose the variable.
  # If not, the variable in "variable" will be plotted.
  
  require(dummies)
  require(kohonen)
  
  interactive <- TRUE
  
  while (interactive == TRUE){
    
    if (variable == 0){
      #show interactive window.
      #color_by_var <- select.list(names(phenotype_data), multiple=FALSE,
      #                            graphics=TRUE, 
      #                            title="Choose variable to color map by.")
      # check for user finished.
      #if (color_by_var == ""){ # if user presses Cancel - we quit function        
      #  return(TRUE)
      #}
      #interactive <- TRUE
      #color_variable <- data.frame(phenotype_data[, color_by_var])
      
      return(plot(som_model, type = "quality", main="Quality score (mean distance of objects to codebook vector of that unit, smaller is better)", palette.name=getColors))
      
    } else {
      color_variable <- data.frame(phenotype_data[, variable])
      color_by_var <- names(phenotype_data)[variable]
      interactive <- FALSE
    }
    
    #if the variable chosen is a string or factor - 
    #Get the levels and ask the user to choose which one they'd like.
    
    if (class(color_variable[,1]) %in% c("character", "factor", "logical")){
      #want to spread this out into dummy factors - but colour by one of those.
      temp_data <- dummy.data.frame(color_variable, sep="_")
      chosen_factor <- select.list(names(temp_data), 
                                   multiple=FALSE,
                                   graphics=TRUE, 
                                   title="Choose level of variable for colouring")
      color_variable <- temp_data[, chosen_factor]
      rm(temp_data, chosen_factor)
      color_by <- color_variable
    } else {      
      #impute the missing values with the mean.
      color_variable[is.na(color_variable[,1]),1] <- mean(color_variable[,1], na.rm=TRUE)
      #color_by <- capVector(color_variable[,1])
      #color_by <- scale(color_by)  
      color_by <- color_variable[,1]
    }
    unit_colors <- aggregate(color_by, by=list(som_model$unit.classif), FUN=mean, simplify=TRUE)
    plot(som_model, type = "property", property=unit_colors[,2], main=color_by_var, palette.name=getColors)    
  }
}

plotFromClusters <- function(clus, label = TRUE, colorCluster = 0, colors = NULL, returnMDS = FALSE, main = "Classical (Metric) Multidimensional Scaling",xlab = "Dimension 1",ylab = "Dimension 2",pch = 16){
  factor <- 1
  numberNodes <- max(unlist(clus))
  n <- (factor*length(clus))+1
  m <- matrix(n,nrow=numberNodes,ncol = numberNodes)
  diag(m)<-0
  
  cols <- rep("black",numberNodes)
  if (colorCluster!=0&colorCluster<n){
    cols[clus[[colorCluster]]] <- "red"
  }
  if (!is.null(colors)){
    cols<-colors
  }
  
  for (cluster in clus){
    for (numOne in seq_len(length(cluster)-1)) {
      for (numTwo in seq(numOne+1,length(cluster))) {
        #print(paste(cluster[[numOne]],cluster[[numTwo]]))
        m[cluster[[numOne]],cluster[[numTwo]]]<-m[cluster[[numOne]],cluster[[numTwo]]]-factor
        m[cluster[[numTwo]],cluster[[numOne]]]<-m[cluster[[numTwo]],cluster[[numOne]]]-factor
      }
    }
  }
  mds <- cmdscale(m, k = 2)
  if (returnMDS){
    return(mds)
  }
  
  lims <- c(
    ll = floor(100 * min(c(mds[, 1], mds[, 2]))) / 100,
    ul = ceiling(100 * max(c(mds[, 1], mds[, 2]))) / 100
  )
  
  size <- 2
  par(mar = c(5,5,4,2)+0.1)
  pl <- plot(mds, 
             xlab = xlab, ylab = ylab,
             xlim = lims, ylim = lims,
             xaxs = "r", yaxs = "r",
             cex.lab = size,
             cex.axis = size,
             cex.main = size,
             main = main, pch = pch, col = cols)
  if (label){
    text(mds, row.names(mds), cex=0.6, pos=4, col="red")
  }
  
  return(pl)
}