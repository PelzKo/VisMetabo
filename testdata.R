testdata <- function(){
  data <- readFile("C:\\Users\\Konstantin\\Desktop\\Uni\\6Semester\\BachelorArbeit\\BeispielDaten\\metabExampleMaleFemale.xlsx")
  inputData <- data$values[complete.cases(data$values), ]
  
  tempId <- 1
  tempStart <- 4
  tempEnd <- 18
  
  
  id <- inputData[[tempId]]
  metab <- inputData[c(tempStart:tempEnd)]
  pheno <- inputData[-c(tempId,c(tempStart:tempEnd))]
  pheno$None <- numeric(length(id[[1]]))
  
  metabComplete <- data.frame(scale(metab))
  metabComplete
}

data(wines)
set.seed(7)

kohmap <- xyf(scale(wines), vintages,grid = somgrid(5, 5, "hexagonal"), rlen=100)

plot(kohmap, type="changes")
counts <- plot(kohmap, type="counts", shape = "straight")

## show both sets of codebook vectors in the map
par(mfrow = c(1,2))
plot(kohmap, type="codes", main = c("Codes X", "Codes Y"))

par(mfrow = c(1,1))
similarities <- plot(kohmap, type="quality", palette.name = terrain.colors)
plot(kohmap, type="mapping",labels = as.integer(vintages), col = as.integer(vintages),main = "mapping plot")
#plot(kohmap, type="mapping",labels = seq_len(177), col = seq_len(177),main = "mapping plot")

## add background colors to units according to their predicted class labels
xyfpredictions <- classmat2classvec(getCodes(kohmap, 2))
bgcols <- c("gray", "pink", "lightgreen")
plot(kohmap, type="mapping", col = as.integer(vintages),pchs = as.integer(vintages), bgcol = bgcols[as.integer(xyfpredictions)],main = "another mapping plot", shape = "straight", border = NA)

## Show'component planes'
set.seed(7)
sommap <- som(scale(wines), grid = somgrid(6, 4, "hexagonal"))
plot(sommap, type = "property", property = getCodes(sommap, 1)[,1],main = colnames(getCodes(sommap, 1))[1])

## Show the U matrix
Umat <- plot(sommap, type="dist.neighbours", main = "SOM neighbour distances")

## use hierarchical clustering to cluster the codebook vectors
som.hc <- cutree(hclust(object.distances(sommap, "codes")), 5)
add.cluster.boundaries(sommap, som.hc)

## and the same for rectangular maps
set.seed(7)
sommap <- som(scale(wines),grid = somgrid(6, 4, "rectangular"))
plot(sommap, type="dist.neighbours", main = "SOM neighbour distances")

## use hierarchical clustering to cluster the codebook vectors
som.hc <- cutree(hclust(object.distances(sommap, "codes")), 5)
add.cluster.boundaries(sommap, som.hc)
