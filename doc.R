runDoc <- function(data){
  initDoc()
  clusterApplier <- .jnew("ClusterApplier")
  arrayDouble <- .jcall(clusterApplier, "[[D", "matrix_from_array",as.vector(as.matrix(data)),ncol(data),evalArray=F)
  result <- .jcall(clusterApplier, "[Lde/lmu/ifi/dbs/elki/data/Cluster;", "doc",arrayDouble,0.2,0.8,0.05,evalArray = F)
  result
}

getIdsDoc <- function(cluster){
  clusterApplier <- .jnew("ClusterApplier")
  ids <- .jcall(clusterApplier, "[[D", "getIds",cluster, simplify = T)
  .jcheck()
  unlist(lapply(ids,function(x){x-min(unlist(ids))+1}))
  
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
  .jaddClassPath("C:\\Users\\Konstantin\\Programmierung\\javaBA\\out\\production\\javaBA\\Clust.jar")
  .jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\elki\\elki-clustering-0.7.5.jar")
  .jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\elki\\elki-core-0.7.5.jar")
  .jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\elki\\elki-core-api-0.7.5.jar")
  .jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\elki\\elki-core-data-0.7.5.jar")
  .jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\elki\\elki-core-dbids-0.7.5.jar")
  .jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\elki\\elki-core-dbids-int-0.7.5.jar")
  .jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\elki\\elki-core-distance-0.7.5.jar")
  .jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\elki\\elki-core-math-0.7.5.jar")
  .jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\elki\\elki-core-util-0.7.5.jar")
  .jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\elki\\elki-database-0.7.5.jar")
  .jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\elki\\elki-input-0.7.5.jar")
  .jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\elki\\elki-logging-0.7.5.jar")
  .jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\dependency\\fastutil-8.2.2.jar")
  .jaddClassPath("C:\\Users\\Konstantin\\Downloads\\elki-0.7.5\\dependency\\jafama-2.3.1.jar")
}
