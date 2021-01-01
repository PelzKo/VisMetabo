library(rJava)

#setwd("C:\\Users\\Konstantin\\Desktop\\Uni\\6Semester\\BachelorArbeit\\code")

runDoc <- function(data, alpha=0.2, beta=0.8, w=0.05){
  data <- apply(data, 2, function(x) {
            if (suppressWarnings(is.na(as.numeric(x[[1]])))){
              return(as.numeric(factor(x)))
            }
            return(as.numeric(x))
          })
  initDoc()
  clusterApplier <- .jnew("ClusterApplier")
  arrayDouble <- .jcall(clusterApplier, "[[D", "matrix_from_array",as.vector(as.matrix(data)),ncol(data),evalArray=F)
  result <- .jcall(clusterApplier, "[Lde/lmu/ifi/dbs/elki/data/Cluster;", "doc",.jarray(as.matrix(data), dispatch = TRUE),alpha,beta,w,evalArray = F)
  result
}

getProtoDoc <- function(cluster){
  if (length(cluster)==0){
    return()
  }
  clusterApplier <- .jnew("ClusterApplier")
  protos <- .jcall(clusterApplier, "[S", "getPrototypeTypes",cluster, simplify = TRUE)
  
}

getIdsDoc <- function(cluster){
  if (length(cluster)==0){
    return()
  }
  clusterApplier <- .jnew("ClusterApplier")
  #ids <- .jcall(clusterApplier, "[[D", "getIds",cluster, simplify = TRUE)
  ids <- .jcall(clusterApplier, "[D", "getIdsOneDim",cluster, simplify = TRUE)
  #if (class(ids)[[1]]!='list'){
  #  ids <- split(ids, rep(1:nrow(ids), each = ncol(ids)))
  #}
  .jcheck()
  listIds <- lapply(ids,function(x){max(x-min(unlist(ids[ids!=0]))+1,0)})
  result <- list()
  counter <- 1
  temp <- numeric()
  for (elem in unlist(listIds)){
    if (elem == 0){
      result[[counter]] <- temp
      counter <- counter + 1
      temp <- numeric()
    } else {
      temp <- c(temp,elem)
    }
  }
  result[[counter]] <- temp
  result
}

getDimsDoc <- function(cluster){
  if (length(cluster)==0){
    return()
  }
  clusterApplier <- .jnew("ClusterApplier")
  dims <- .jcall(clusterApplier, "[[D", "getDims",cluster, simplify = T)
  .jcheck()
  dims
}

getAvgsDoc <- function(cluster){
  if (length(cluster)==0){
    return()
  }
  clusterApplier <- .jnew("ClusterApplier")
  avgs <- .jcall(clusterApplier, "[[D", "getAverages",cluster, simplify = T)
  .jcheck()
  avgs
}

initDoc <- function(){
  .jinit(".")
  .jaddClassPath("dependencies/elki-clustering-0.7.5.jar")
  .jaddClassPath("dependencies/elki-core-0.7.5.jar")
  .jaddClassPath("dependencies/elki-core-api-0.7.5.jar")
  .jaddClassPath("dependencies/elki-core-data-0.7.5.jar")
  .jaddClassPath("dependencies/elki-core-dbids-0.7.5.jar")
  .jaddClassPath("dependencies/elki-core-dbids-int-0.7.5.jar")
  .jaddClassPath("dependencies/elki-core-distance-0.7.5.jar")
  .jaddClassPath("dependencies/elki-core-math-0.7.5.jar")
  .jaddClassPath("dependencies/elki-core-util-0.7.5.jar")
  .jaddClassPath("dependencies/elki-database-0.7.5.jar")
  .jaddClassPath("dependencies/elki-input-0.7.5.jar")
  .jaddClassPath("dependencies/elki-logging-0.7.5.jar")
  .jaddClassPath("dependencies/fastutil-8.2.2.jar")
  .jaddClassPath("dependencies/jafama-2.3.1.jar")
  .jaddClassPath("dependencies/Clust.jar")
}
