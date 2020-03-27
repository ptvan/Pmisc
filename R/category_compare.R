#' compares the mean expression between two groups of genes in a given
#' expression matrix, calculates the number and proportion of genes in the
#' input matrix that are up-regulated
#'
#' @param dat A matrix of gene expression
#' @param setsIndices A nested-list of gene indices for each GSEA category, most commonly made from limma::ids2indices()
#' @param grp1name String indicating name of the first group for the comparison
#' @param grp2name String indicating name of the second group for the comparison
#' @param grp1idx A vector of integer indices indicating columns belonging to group1
#' @param grp2idx A vector of integer indices indicating columns belonging to group1
#' @return table of genes' mean expression values for each group
#' @examples
#' \dontrun{
#' out <- category_compare(expressionMatrix,
#'                                   setIndices,
#'                                   grp1name = "control",
#'                                   grp2name = "experimental",
#'                                   which(grepl("control", colnames(expressionMatrix))),
#'                                   which(grepl("experimental", colnames(expressionMatrix)))
#'                                   )
#' }
#' @import data.table


category_compare <- function(dat, setsIndices, grp1name="group1", grp2name="group2", grp1idx, grp2idx){
  if (!is.null(names(setsIndices))){

    categories <- names(setsIndices)

    tab <- data.table(cbind(categories, rep(1, length(categories)), rep(0, length(categories))))
    setnames(tab, c("categories", "V2","V3")
             , c("category", "totalGenes", "genesUp"))
    tab$totalGenes <- as.numeric(tab$totalGenes)
    tab$genesUp <- as.numeric(tab$genesUp)

    for (i in 1:length(categories)){
      cat <- categories[i]
      cidxs <- unlist(setsIndices[cat], use.names = F)
      tab[tab$category==cat]$totalGenes <- length(cidxs)
      d <- data.table(dat[cidxs,])

      if (length(grp1idx) + length(grp2idx) == ncol(d) && !identical(grp1idx, grp2idx)){
        rows <- rownames(d)
        cols <- colnames(d)
        idxs <- cols
        idxs[grp1idx] <- grp1name
        idxs[grp2idx] <- grp2name
        idx <- data.frame(cbind(cols,idxs))
        d <- merge(melt(d, measure.vars=colnames(d)), idx, by.x="variable", by.y="cols")
        setnames(d, c("Var1","Var2","value", "idxs")
                 , c("geneName", "inputColumn","expr","group"))

        row <- merge(subset(d, group==grp1name)[,list(group1expr=mean(expr)),by=geneName],
                     subset(d, group==grp2name)[,list(group2expr=mean(expr)),by=geneName],
                     by="geneName")
        tab[tab$category==cat]$genesUp <- length(which(row$group2expr > row$group1expr))
      }

    }

    tab$proportionUp <- round(tab$genesUp / tab$totalGenes, 3)
    setnames(tab, c("genesUp"), c(paste0("genesUpIn", grp2name)))
    return(tab)
  } else { "setsIndices must be a named list-of-lists !!!" }

}
