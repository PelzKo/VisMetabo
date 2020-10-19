library(subspace)
source("C:\\Users\\Konstantin\\Desktop\\Uni\\6Semester\\BachelorArbeit\\code\\testdata.R")

N <- 1000  # total number of rows to preallocate--possibly an overestimate

result <- data.frame(name=rep("", N), dimNames=rep("", N),  # as many cols as you need
                 diab1=rep(NA, N),diab2=rep(NA, N),diab3=rep(NA, N),diab4=rep(NA, N),
                 diab5=rep(NA, N)) 


small <- getData()
smallMetab <- small[[1]]
smallPheno <- small[[2]]
T2D <- smallPheno["T2D"]
#big <- getData(FALSE)
#bigMetab <- big[[1]]
#bigPheno <- big[[2]]
counter<-1
currentCombinationCounter <- 1


method <- "clique"
xi <- 40
tau <- 0.1

combinations <- expand.grid(c(40,60,100), c(0.1,0.2,0.4))

cliqueClust <- CLIQUE(smallMetab, xi = xi,tau = tau)
cliqueClust <- cliqueClust[order(unlist(lapply(lapply(cliqueClust, "[[", "subspace"),"sum")),decreasing = TRUE)]
nums <- unlist(lapply(lapply(cliqueClust, "[[", "subspace"),"sum"))
top5 <- nums[1:5]
top5Metabs <- list()
top5Diab <- rep(0.5,5)


for (i in 1:5){
  currentCluster <- cliqueClust[as.numeric(i)][[1]]
  top5Metabs[[i]] <- names(smallMetab)[currentCluster$subspace]
  top5Diab[i] <- mean(T2D[[1]][currentCluster$objects])
}

realDiabs <- abs(0.5-top5Diab)

result[counter, ] <- list(paste(method,xi,tau,sep = "_"),do.call("paste", c(top5Metabs, sep = "_")) ,top5Diab[[1]],top5Diab[[2]],top5Diab[[3]],top5Diab[[4]],top5Diab[[5]])


result <- result[complete.cases(result),]

write.table(result, file='output.tsv', quote = FALSE, row.names = FALSE, sep='\t')
