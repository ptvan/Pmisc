#' this is a convenience wrapper to make plotting PCA biplots a bit nicer
#' it removes data columns that have zero variance, which do not 
#' contribute to the PCs but interferes with plotting
#' it also colors the groups with user-supplied labels

#' @param df data.frame of input data
#' @param title text to title the output plot
#' @param labels labels to label the groups
#' @param varname.adjust parameter to be passed to `ggbiplot`

#' @return a ggplot2 plot
#' @import ggbiplot

plotpca <- function (df, title=NULL, labels=NULL, varname.adjust=2, ...) {
  # remove data columns that have zero variance
  df <- df[,apply(df,2,var, na.rm=TRUE) != 0]
  # run PCA
  pcobj <- prcomp(df)
  if (length(labels) != nrow(df)){
    cat("length of `labels` did not match rows of input, groups will not be colored in biplot !\n")
    labels <- NULL
  }
  # plotting
  if (is.null(labels)){
    # no labels
    out <- ggbiplot(pcobj,...) +
      ggtitle(title)
    
  } else {
    
    out <- ggbiplot(pcobj, groups=labels, ...) +
      ggtitle(title) +
      geom_text(aes(label=labels))
  }
  print(out)  
  return(out)
}