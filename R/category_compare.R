category_compare <- function(dat, setsIndices, grp1name="group1", grp2name="group2", grp1idx, grp2idx){
  # compares the mean expression between two groups of genes in a given
  # expression matrix, calculates the number and proportion of genes in the
  # input matrix that are up-regulated

  # requires a setsindices, a nested-list of gene indices for each GSEA category
  # setsindices are most commonly made from limma::ids2indices()

  # example: setsIndices <- ids2indices(geneIds, rownames(expressionMatrix))
  # categories <- names(setsIndices)
  # out <- category_compare(expressionMatrix,
  #                                   setIndices,
  #                                   grp1name = "control",
  #                                   grp2name = "experimental",
  #                                   which(grepl("control", colnames(expressionMatrix))),
  #                                   which(grepl("experimental", colnames(expressionMatrix)))
  #                                   )
  require(data.table)

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
      d <- dat[cidxs,]

      if (length(grp1idx) + length(grp2idx) == ncol(d) && !identical(grp1idx, grp2idx)){
        rows <- rownames(d)
        cols <- colnames(d)
        idxs <- cols
        idxs[grp1idx] <- grp1name
        idxs[grp2idx] <- grp2name
        idx <- data.frame(cbind(cols,idxs))
        d <- data.table(merge(melt(d), idx, by.x="Var2", by.y="cols"))
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
