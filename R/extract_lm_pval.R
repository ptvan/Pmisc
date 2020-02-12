#' extracts F-statistic from a linear model object
#' then calculates P-value
#' @param model an lm object
#' @return pvalues
#' @examples
#' \dontrun{
#' pvals <- extract_lm_pval(myModelObj)
#' }


extract_lm_pval <- function (model){
  if ( class(model) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(model)$fstatistic
  p <- pf(f[1],f[2],f[3], lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
