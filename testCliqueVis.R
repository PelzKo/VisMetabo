library(maps)

data(us.cities)

#create distance matrix
d <- dist(us.cities[, c("lat", "long")])

#multidimensional scaling so we can plot and retain distance relationships
foo <- cmdscale(d, k = 2)

#everything is upside down and backwards
#plot(foo)

#plot(-foo)
#plot(-foo, cex = scale(us.cities$pop))

test <- list(c(1,2,3,4),c(1,3,4),c(1,2))

plotFromClusters <- function(clus, label = TRUE, colorCluster = 0, colors = NULL){
  numberNodes <- max(unlist(clus))
  n <- length(clus)+1
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
        m[cluster[[numOne]],cluster[[numTwo]]]<-m[cluster[[numOne]],cluster[[numTwo]]]-1
        m[cluster[[numTwo]],cluster[[numOne]]]<-m[cluster[[numTwo]],cluster[[numOne]]]-1
      }
    }
  }
  mds <- cmdscale(m, k = 2)
  pl <- plot(mds, pch = 16, col = cols)
  if (label){
    text(mds, row.names(mds), cex=0.6, pos=4, col="red")
  }
  return(pl)
}
