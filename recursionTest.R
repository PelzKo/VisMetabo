library(rlist)
getClusters <- function(dendro,limit=10){
  if (nobs(dendro)<limit){
    return(NULL)
  }
  return(list(unlist(dendro),getClusters(dendro[[1]],limit),getClusters(dendro[[2]],limit)))
}
