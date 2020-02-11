camera2 <- function (y, index, design = NULL, contrast = ncol(design), weights = NULL,
                     use.ranks = FALSE, allow.neg.cor = FALSE, inter.gene.cor = 0.01,
                     trend.var = FALSE, sort = TRUE, ...)
{
  require(limma)
  dots <- names(list(...))
  if (length(dots))
    warning("Extra arguments disregarded: ", sQuote(dots))
  y <- getEAWP(y)
  G <- nrow(y$exprs)
  n <- ncol(y$exprs)
  ID <- rownames(y$exprs)
  if (G < 3)
    stop("Two few genes in dataset: need at least 3")
  if (!is.list(index))
    index <- list(set1 = index)
  nsets <- length(index)
  if (nsets == 0L)
    stop("index is empty")
  if (is.null(design))
    design <- y$design
  if (is.null(design))
    stop("design matrix not specified")
  else {
    design <- as.matrix(design)
    if (mode(design) != "numeric")
      stop("design must be a numeric matrix")
  }
  if (nrow(design) != n)
    stop("row dimension of design matrix must match column dimension of data")
  p <- ncol(design)
  df.residual <- n - p
  if (df.residual < 1)
    stop("No residual df: cannot compute t-tests")
  if (is.null(weights))
    weights <- y$weights
  fixed.cor <- !(is.na(inter.gene.cor) || is.null(inter.gene.cor))
  if (fixed.cor) {
    if (use.ranks)
      df.camera <- Inf
    else df.camera <- G - 2
  }
  else {
    df.camera <- min(df.residual, G - 2)
  }
  y <- y$exprs
  if (!is.null(weights)) {
    if (any(weights <= 0))
      stop("weights must be positive")
    if (length(weights) == n) {
      sw <- sqrt(weights)
      y <- t(t(y) * sw)
      design <- design * sw
      weights <- NULL
    }
  }
  if (!is.null(weights)) {
    if (length(weights) == G)
      weights <- matrix(weights, G, n)
    weights <- as.matrix(weights)
    if (any(dim(weights) != dim(y)))
      stop("weights not conformal with y")
  }
  if (is.character(contrast)) {
    contrast <- which(contrast == colnames(design))
    if (length(contrast) == 0)
      stop("coef ", contrast, " not found")
  }
  if (length(contrast) == 1) {
    j <- c((1:p)[-contrast], contrast)
    if (contrast < p)
      design <- design[, j]
  }
  else {
    QR <- qr(contrast)
    design <- t(qr.qty(QR, t(design)))
    if (sign(QR$qr[1, 1] < 0))
      design[, 1] <- -design[, 1]
    design <- design[, c(2:p, 1)]
  }
  if (is.null(weights)) {
    QR <- qr(design)
    if (QR$rank < p)
      stop("design matrix is not of full rank")
    effects <- qr.qty(QR, t(y))
    unscaledt <- effects[p, ]
    if (QR$qr[p, p] < 0)
      unscaledt <- -unscaledt
  }
  else {
    effects <- matrix(0, n, G)
    colnames(effects) <- ID
    unscaledt <- rep.int(0, G)
    names(unscaledt) <- ID
    sw <- sqrt(weights)
    yw <- y * sw
    for (g in 1:G) {
      xw <- design * sw[g, ]
      QR <- qr(xw)
      if (QR$rank < p)
        stop("weighted design matrix not of full rank for gene ",
             g)
      effects[, g] <- qr.qty(QR, yw[g, ])
      unscaledt[g] <- effects[p, g]
      if (QR$qr[p, p] < 0)
        unscaledt[g] <- -unscaledt[g]
    }
  }
  U <- effects[-(1:p), , drop = FALSE]
  sigma2 <- colMeans(U^2)
  U <- t(U)/sqrt(pmax(sigma2, 1e-08))
  if (trend.var)
    A <- rowMeans(y)
  else A <- NULL
  sv <- squeezeVar(sigma2, df = df.residual, covariate = A)
  modt <- unscaledt/sqrt(sv$var.post)
  if (use.ranks)
    Stat <- modt
  else {
    df.total <- min(df.residual + sv$df.prior, G * df.residual)
    Stat <- zscoreT(modt, df = df.total, approx = TRUE)
  }
  meanStat <- mean(Stat)
  varStat <- var(Stat)
  geneVals <- rep(list(list()), length(index))
  names(geneVals) <- names(index)
  tab <- matrix(0, nsets, 5)
  rownames(tab) <- names(index)
  colnames(tab) <- c("NGenes", "Correlation", "Down", "Up",
                     "TwoSided")
  for (i in 1:nsets) {
    iset <- index[[i]]
    if (is.character(iset))
      iset <- which(ID %in% iset)
    StatInSet <- Stat[iset]
    m <- length(StatInSet)
    m2 <- G - m
    if (fixed.cor) {
      correlation <- inter.gene.cor
      vif <- 1 + (m - 1) * correlation
    }
    else {
      if (m > 1) {
        Uset <- U[iset, , drop = FALSE]
        vif <- m * mean(colMeans(Uset)^2)
        correlation <- (vif - 1)/(m - 1)
      }
      else {
        vif <- 1
        correlation <- NA
      }
    }
    tab[i, 1] <- m
    tab[i, 2] <- correlation
    if (use.ranks) {
      if (!allow.neg.cor)
        correlation <- max(0, correlation)
      tab[i, 3:4] <- rankSumTestWithCorrelation(iset, statistics = Stat,
                                                correlation = correlation, df = df.camera)
    }
    else {
      if (!allow.neg.cor)
        vif <- max(1, vif)
      meanStatInSet <- mean(StatInSet)
      delta <- G/m2 * (meanStatInSet - meanStat)
      varStatPooled <- ((G - 1) * varStat - delta^2 * m *
                          m2/G)/(G - 2)
      two.sample.t <- delta/sqrt(varStatPooled * (vif/m +
                                                    1/m2))
      tab[i, 3] <- pt(two.sample.t, df = df.camera)
      tab[i, 4] <- pt(two.sample.t, df = df.camera, lower.tail = FALSE)
    }
    geneVals[[i]] <- StatInSet
  }
  tab[, 5] <- 2 * pmin(tab[, 3], tab[, 4])
  tab <- data.frame(tab, stringsAsFactors = FALSE)
  Direction <- rep.int("Up", nsets)
  Direction[tab$Down < tab$Up] <- "Down"
  tab$Direction <- Direction
  tab$PValue <- tab$TwoSided
  tab$Down <- tab$Up <- tab$TwoSided <- NULL
  if (fixed.cor)
    tab$Correlation <- NULL
  if (nsets > 1)
    tab$FDR <- p.adjust(tab$PValue, method = "BH")
  if (sort && nsets > 1) {
    o <- order(tab$PValue)
    tab <- tab[o, ]
  }

  foo <- list(tab, geneVals)
  names(foo) <- c("summary", "geneVals")
  foo
}

geneCompareWithWeights <- function (dat, camera2output, genesetidx, grp1name="group1", grp2name="group2", grp1idx, grp2idx){
  # compares the mean expression between two groups of genes in a given
  # expression matrix, optionally makes boxplots of the comparison, faceted by gene
  # used in a pipeline, so checks are pretty rudimentary

  # example: out <- geneCompare(expressionMatrix
  #                                    ,c("LAP3","NCOA7","IFIT2","STAT1","TRIM21")
  #                                   ,grp1name = "control"
  #                                   ,grp2name = "experimental"
  #                                   ,which(grepl("control", colnames(expressionMatrix)))
  #                                   which(grepl("experimental", colnames(expressionMatrix)))
  #                                   )

  require(data.table)

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
      d <- data.table(merge(melt(dat), idx, by.x="Var2", by.y="cols"))
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


      # if(plot) {
      #   ggplot(d, aes(x=group, y=expr, col=group)) +
      #     geom_boxplot(outlier.shape = NA) +
      #     facet_wrap(~ geneName, scales="free") +
      #     theme(legend.position = "none") +
      #     labs(title=paste0(length(rows), " genes in list"))
      # }

      return(tab)

    } else { stop ("must specify one unique index for each column of data !!!") }

    # colnames(d) <- c("gene", "group", "expression")
  } else { stop ("geneSet does not match ANY columns in data !!!") }
}
