## A function to use identify to select points, and overplot the
## points with another symbol as they are selected
identifyPch <- function(x, y = NULL, n = length(x), plot = FALSE, pch = 19, ...)
{
  xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
  sel <- rep(FALSE, length(x))
  while(sum(sel) < n) {
    ans <- identify(x[!sel], y[!sel], labels = which(!sel), n = 1, plot = plot, ...)
    if(!length(ans)) break
    ans <- which(!sel)[ans]
    points(x[ans], y[ans], pch = pch)
    sel[ans] <- TRUE
  }
  ## return indices of selected points
  which(sel)
}

if(dev.interactive()) { ## use it
  x <- rnorm(50); y <- rnorm(50)
  plot(x,y); identifyPch(x,y) # how fast to get all?
}


somData$codes[[1]][27,]
vecs <- identify(somData,plot=FALSE, n=1)
somData$unit.classif==27
seq_len(1756)[somData$unit.classif==27]


#install.packages("https://cran.r-project.org/src/contrib/Archive/rJava/rJava_0.9-12.tar.gz", repos=NULL, method="libcurl")
#arr <- c(c(1,1,1,1),c(1,1,1,1),c(9,9,9,9),c(1,1,1,1),c(1,1,1,1),c(9,9,9,9))
arr <- c(1,1,1,1,9,9,1,1,1,1,9,9,1,1,1,1,9,9,1,2,3,4,5,6)
dataM <- matrix(arr,6,4)
data <- testdata()#data.frame(dataM)
test <- as.matrix(arr)
library(rJava)
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
.jclassPath()
clusterAppl <- .jnew("ClusterApplier")
.jmethods(clusterAppl)
#stringArrayTest <- .jcall(clusterAppl, "[S", "getStringArray")
arrayDouble <- .jcall(clusterAppl, "[[D", "matrix_from_array",as.vector(as.matrix(data)),ncol(data),evalArray=F)
#res <- .jcall(clusterAppl, returnSig="[[Ljava/lang/String", method="doc",arrayDouble,0.2,0.8,0.05,evalArray=T)
res <- .jcall(clusterAppl, "[Lde/lmu/ifi/dbs/elki/data/Cluster;", "doc",arrayDouble,0.2,0.8,0.05,evalArray = F)
ids <- .jcall(clusterAppl, "[[D", "getIds",res, simplify = T)
dims <- .jcall(clusterAppl, "[[D", "getDims",res, simplify = T)
avg <- .jcall(clusterAppl, "[[D", "getAverages",res, simplify = T)
lapply(ids,function(x){x-min(unlist(ids))+1})

res <- rJava::.jcall("ClusteringApplier",returnSig="[Lde.lmu.ifi.dbs.elki.data.Cluster;",method="doc",arr,0.2,0.8,0.05,evalArray=F)
.rs.restartR()
