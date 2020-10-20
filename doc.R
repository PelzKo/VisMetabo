library(rJava)

runDoc <- function(data, alpha=0.2, beta=0.8, w=0.05){
  initDoc()
  clusterApplier <- .jnew("ClusterApplier")
  arrayDouble <- .jcall(clusterApplier, "[[D", "matrix_from_array",as.vector(as.matrix(data)),ncol(data),evalArray=F)
  result <- .jcall(clusterApplier, "[Lde/lmu/ifi/dbs/elki/data/Cluster;", "doc",arrayDouble,alpha,beta,w,evalArray = F)
  result
}

getIdsDoc <- function(cluster){
  clusterApplier <- .jnew("ClusterApplier")
  ids <- .jcall(clusterApplier, "[[D", "getIds",cluster, simplify = T)
  if (class(ids)[[1]]!='list'){
    ids <- split(ids, rep(1:nrow(ids), each = ncol(ids)))
  }
  .jcheck()
  lapply(ids,function(x){x-min(unlist(ids))+1})
  
}

getDimsDoc <- function(cluster){
  clusterApplier <- .jnew("ClusterApplier")
  dims <- .jcall(clusterApplier, "[[D", "getDims",cluster, simplify = T)
  .jcheck()
  dims
}

getAvgsDoc <- function(cluster){
  clusterApplier <- .jnew("ClusterApplier")
  avgs <- .jcall(clusterApplier, "[[D", "getAverages",cluster, simplify = T)
  .jcheck()
  avgs
}

initDoc <- function(){
  .jinit()
  .jaddClassPath("dependencies/Clust.jar")
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
  
  #.jaddClassPath("C:\\Users\\Konstantin\\Programmierung\\javaBA\\out\\production\\javaBA\\Clust.jar")
  #.jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\elki\\elki-clustering-0.7.5.jar")
  #.jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\elki\\elki-core-0.7.5.jar")
  #.jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\elki\\elki-core-api-0.7.5.jar")
  #.jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\elki\\elki-core-data-0.7.5.jar")
  #.jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\elki\\elki-core-dbids-0.7.5.jar")
  #.jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\elki\\elki-core-dbids-int-0.7.5.jar")
  #.jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\elki\\elki-core-distance-0.7.5.jar")
  #.jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\elki\\elki-core-math-0.7.5.jar")
  #.jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\elki\\elki-core-util-0.7.5.jar")
  #.jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\elki\\elki-database-0.7.5.jar")
  #.jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\elki\\elki-input-0.7.5.jar")
  #.jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\elki\\elki-logging-0.7.5.jar")
  #.jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\dependency\\fastutil-8.2.2.jar")
  #.jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\dependency\\jafama-2.3.1.jar")
}
