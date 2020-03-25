#' compares the mean expression between two groups of genes in a given
#' expression matrix, optionally makes boxplots of the comparison, faceted by gene
#' used in a pipeline, so checks are pretty rudimentary
#'
#' @param dat A matrix of gene expression
#' @param camera2output output matrix from camera2()
#' @param genesetidx a list of lists of gene indices
#' @param grp1name String indicating name of the first group for the comparison
#' @param grp2name String indicating name of the second group for the comparison
#' @param grp1idx A vector of integer indices indicating columns belonging to group1
#' @param grp2idx A vector of integer indices indicating columns belonging to group1
#' @return table of genes' mean expression values for each group
#' @examples
#' \dontrun{
#' out <- gene_compare(expressionMatrix
#'                    ,c("LAP3","NCOA7","IFIT2","STAT1","TRIM21")
#'                    ,grp1name = "control"
#'                    ,grp2name = "experimental"
#'                    ,which(grepl("control", colnames(expressionMatrix)))
#'                    ,which(grepl("experimental", colnames(expressionMatrix)))
#'                                   )
#' }
#' @import ggplot2
#' @import data.table
#' @import reshape2

gene_compare_with_weights <- function (dat, camera2output, genesetidx, grp1name="group1", grp2name="group2", grp1idx, grp2idx){

  if(!is.list(camera2output) || length(camera2output) != 2 ) stop ("requires a list of length 2 ! ")
  if(length(genesetidx) != 1) stop ("requires exactly 1 genesetidx !")

  geneVals <- camera2output$geneVals
  summary <- camera2output$summary

  geneSet <- names(geneVals[[genesetidx]])
  geneSetName <- names(geneVals)[genesetidx]
  geneSetPValue <- summary[geneSetName,]$PValue
  geneSetFDR <- summary[geneSetName,]$FDR

  geneWeight <- geneVals[[genesetidx]]
  geneWeight <- data.frame(geneWeight)
  geneWeight$geneName <- rownames(geneWeight)

  # check that at least some of data has genes in the specified geneSet
  if (length(intersect(geneSet, rownames(dat))) > 0){

    dat <- subset(dat, rownames(dat) %in% geneSet)
    # check that all columns have group assignments and the groups don't overlap
    if (length(grp1idx) + length(grp2idx) == ncol(dat) && !identical(grp1idx, grp2idx)){
      rows <- rownames(dat)
      cols <- colnames(dat)
      idxs <- cols
      idxs[grp1idx] <- grp1name
      idxs[grp2idx] <- grp2name
      idx <- data.frame(cbind(cols,idxs))
      d <- data.table(merge(reshape2::melt(dat), idx, by.x="Var2", by.y="cols"))
      setnames(d, c("Var1","Var2","value", "idxs")
               , c("geneName", "inputColumn","expr","group"))

      tab <- merge(subset(d, group==grp1name)[,list(group1expr=mean(expr)),by=geneName],
                   subset(d, group==grp2name)[,list(group2expr=mean(expr)),by=geneName],
                   by="geneName")

      tab <- merge(tab, geneWeight, by="geneName")

      tab$higher <- tab$group2expr > tab$group1expr
      tab$higher <- gsub("FALSE", grp1name, tab$higher)
      tab$higher <- gsub("TRUE", grp2name, tab$higher)

      setnames(tab, c("group1expr", "group2expr")
               , c(grp1name,grp2name))

      attr(tab, "geneSetName") <- geneSetName
      attr(tab, "geneSetPValue") <- geneSetPValue
      attr(tab, "geneSetFDR") <- geneSetFDR

      return(tab)

    } else { stop ("must specify one unique index for each column of data !!!") }

    # colnames(d) <- c("gene", "group", "expression")
  } else { stop ("geneSet does not match ANY columns in data !!!") }
}
