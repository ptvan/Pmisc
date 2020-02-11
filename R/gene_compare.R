gene_compare <- function (dat, geneSet, grp1name="group1", grp2name="group2", grp1idx, grp2idx, plot=FALSE){
  # compares the mean expression between two groups of genes in a given
  # expression matrix, optionally makes boxplots of the comparison, faceted by gene
  # used in a pipeline, so checks are pretty rudimentary

  # example: out <- gene_compare(expressionMatrix
  #                                    ,c("LAP3","NCOA7","IFIT2","STAT1","TRIM21")
  #                                   ,grp1name = "control"
  #                                   ,grp2name = "experimental"
  #                                   ,which(grepl("control", colnames(expressionMatrix)))
  #                                   which(grepl("experimental", colnames(expressionMatrix)))
  #                                   )

  require(ggplot2)
  require(data.table)

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
      d <- data.table(merge(melt(dat), idx, by.x="Var2", by.y="cols"))
      setnames(d, c("Var1","Var2","value", "idxs")
               , c("geneName", "inputColumn","expr","group"))

      tab <- merge(subset(d, group==grp1name)[,list(group1expr=mean(expr)),by=geneName],
                   subset(d, group==grp2name)[,list(group2expr=mean(expr)),by=geneName],
                   by="geneName")

      tab$higher <- tab$group2expr > tab$group1expr
      tab$higher <- gsub("FALSE", grp1name, tab$higher)
      tab$higher <- gsub("TRUE", grp2name, tab$higher)

      setnames(tab, c("group1expr", "group2expr")
               , c(grp1name,grp2name))

      if(plot) {
        ggplot(d, aes(x=group, y=expr, col=group)) +
          geom_boxplot(outlier.shape = NA) +
          facet_wrap(~ geneName, scales="free") +
          theme(legend.position = "none") +
          labs(title=paste0(length(rows), " genes in list"))
      }
      return(tab)

    } else { stop ("must specify one unique index for each column of data !!!") }

    colnames(d) <- c("gene", "group", "expression")
  } else { stop ("geneSet does not match ANY columns in data !!!") }
}
