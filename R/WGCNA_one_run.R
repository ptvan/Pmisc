WGCNA_one_run <- function(dat, netType="unsigned", pow=NULL, iter=1, defaultPow=3, showPlots=TRUE, ds=1) {
  ## modified from original code by Carl Murie
  require(WGCNA)

  ## calculate power estimate
  powers <- c(1:10, seq(from=12, to=20, by=2))
  sft <- pickSoftThreshold(dat, powerVector=powers, networkType=netType, verbose=5)

  if(showPlots) {
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="SoftThreshold(power)",
         ylab="ScaleFreeTopologyModelFit,signedR^2",type="n",main =paste("Scale independence"))
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=1,col="red")
  } ## end if showPlots

  ## set power parameter if not passed - if power estimate is NA then use defaultPow
  if(is.null(pow)) {
    if(!is.na(sft$powerEstimate)) {
      pow <- sft$powerEstimate
    } else {
      pow <- defaultPow
    }
  } ## end if is.null

  cat("Executing WGCNA with power=", pow, "\n")

  adjacency <- adjacency(dat, power=pow)                  ## calculate adjacency matrix of pearson correlations
  TOM <- TOMsimilarity(adjacency, TOMType=netType)        ## calculate similarity matrix
  dissTOM <- 1- TOM
  geneTree <- hclust(as.dist(dissTOM), method="average")  ## cluster on similarity matrix

  ## calculate tree cuts
  dynamicMods <- cutreeDynamic(dendro=geneTree, distM=dissTOM, deepSplit=ds, pamStage=FALSE,
                               pamRespectsDendro=FALSE, minClusterSize=20)
  dynamicColors <- labels2colors(dynamicMods)

  if(showPlots) {
    plotDendroAndColors(geneTree, dynamicColors, dendroLabels=FALSE, addGuide=TRUE)
  } ## end if showPlots

  MEList <- moduleEigengenes(dat, colors=dynamicColors)
  MEs <- MEList$eigengenes
  MEDiss <- 1 - cor(MEs)
  METree <- hclust(as.dist(MEDiss), method="average")
  MEDissThres <- 0.15

  if(showPlots) {
    plot(METree)
    abline(h=MEDissThres, col="red")
  } ## end showPlots

  merge <- mergeCloseModules(dat, dynamicColors, cutHeight=MEDissThres, verbose=3)
  mergedColors <- merge$colors
  mergedMEs <- merge$newMEs

  if(showPlots) {
    plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), dendroLabels=FALSE, addGuide=TRUE)
  } ## end if showPlots

  return(merge)
} ## end WGCNAOneRun
