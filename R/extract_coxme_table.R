#' extracts model output from an lmekin object
#' @param model an lmekin object
#' @return data.frame containing beta's, SE's, Z-values and P for each contrast
#' @examples
#' \dontrun{
#' out <- extract_coxme(mylmekinObj)
#' }

# from https://stackoverflow.com/questions/43720260/how-to-extract-p-values-from-lmekin-objects-in-coxme-package
extract_coxme_table <- function (model){
  if ( class(model) != "lmekin") stop("Not an object of class 'lmekin' ")
  beta <- model$coefficients$fixed
  nvar <- length(beta)
  nfrail <- nrow(model$var) - nvar
  se <- sqrt(diag(model$var)[nfrail + 1:nvar])
  z <- round(beta/se, 2)
  p <- signif(1 - pchisq((beta/se)^2, 1), 2)
  table <- data.frame(cbind(beta,se,z,p))
  return(table)
}
