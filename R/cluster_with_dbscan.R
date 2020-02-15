#' cluster using DBscan 

#' @param dat data matrix
#' @param epsilon param to be passed to `dbscan()`
#' @param minpts param to be passed to `dbscan()`
#' 
#' @import dbscan
#' @import data.table

cluster_with_dbscan <- function(dat, epsilon, minpts = 5 ) {
  out <- dbscan(dat[,.(x,y)], eps = epsilon)
  dat[,"dbscan_cluster"] <- out$cluster
  return(dat)
}