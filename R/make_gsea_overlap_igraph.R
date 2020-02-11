#' this function takes a bioConductor EList (eg. voom-transformed RNASeq counts)
#' GSEA output matrix from limma's camera() and a user-specified geneSet
#' and create an igraph graph where the nodes are gene lists and the edge weight
#' is how many expressed genes are common between the gene lists
#' Phu T. Van, FHCRC 2017, w/ substantial help and input from C.Murie & V.Voillet
#'
#' @param expressionList a matrix of gene expression values
#' @param cameraMat matrix contain output from CAMERA
#' @param geneSets a list of list of geneIds
#' @param verbose Boolean whether to print out progress
#' @return table of genes' mean expression values for each group
#' @examples
#' \dontrun{
#' vDat <- voom(exprs(eDat), design=designMat, plot=FALSE, lib.size=libNorm)
#' res <- camera(vDat, setsIndices, design=designMat, contrast=cons[i], sort=TRUE)
#' geneSets <- geneIds(getGmt(gmtFile))
#' n <- make_gsea_igraph(vDat, res, geneSets)
#' plot(n, edge.width=E(n)$overlap*.1, vertex.color="white", vertex.label.cex=0.5)
#' }
#' @import igraph

make_gsea_overlap_igraph <- function(expressionList, cameraMat, geneSets, verbose=FALSE){
  # get the geneset categories, and the expressed genes
  geneSetCats <- rownames(cameraMat)
  expressedGenes <- rownames(expressionList$E)

  # create unconnected graph with geneSet categories as vertices
  net <- make_empty_graph(n = nrow(cameraMat), directed=FALSE)

  # add vertex attributes
  V(net)$label <- geneSetCats
  V(net)$totalGeneCount <- unlist(lapply(geneSets[rownames(cameraMat)], length))
  V(net)$sigGeneCount <- cameraMat$NGenes
  V(net)$direction <- cameraMat$Direction
  V(net)$FDR <- cameraMat$FDR
  V(net)$pval <- cameraMat$PValue

  # create list of edges
  tmp <- list()
  for (i in 1:length(geneSetCats)){
    if (verbose){
      cat("working on", geneSetCats[i], ":", length(unlist(geneSets[geneSetCats[i]])), "total genes,"
          ,length(which(expressedGenes%in%unlist(geneSets[geneSetCats[i]]))), "of which were expressed \n"
      )}
    tmp[[i]] <- expressedGenes[which(expressedGenes%in%unlist(geneSets[geneSetCats[i]]))]
    names(tmp)[i] <- geneSetCats[i]
  }
  nms <- combn( names(tmp), 2 , FUN = paste0 ,collapse = "xxx" ,simplify = FALSE )
  ll <- combn( tmp , 2 , simplify = FALSE )
  edges <- lapply( ll , function(x)  intersect( x[[1]] , x[[2]] )  )
  names(edges) <- nms
  edges <- edges[lapply(edges, length)>0]

  # add edges from edge list to  graph
  for (i in 1:length(edges)){
    nodes <- unlist(str_split(names(edges)[i],"xxx"))
    net <- add_edges(net, c(which(nodes[1] == V(net)$label), which(nodes[2] == V(net)$label)))
    net <- set_edge_attr(net, "overlap", index=i, value=length(edges[[i]]))
  }

  return(net)
}
