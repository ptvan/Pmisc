#' plots the output of a tsne run on ICS flow data, highlights user-specified window, and return
#' cell events within that window. Also perform t.test() on extracted events

#' @param dat RTsne output matrix
#' @param x_min coordinate of bounding box to extract events
#' @param x_max coordinate of bounding box to extract events
#' @param y_min coordinate of bounding box to extract events
#' @param y_max coordinate of bounding box to extract events
#' @param antigen the antigen to subset, must exist in `dat`
#' @param marker the marker to facet on
#' @param polyDegree ICS polyfunctionality degree of data to extract
#' @param facet metadata to facet on, must exist in `dat`
#' @param title plot title
#'
#' @import ggplot2
#' @import gridExtra
#' @import RColorBrewer
#' @import sp
#' @import pheatmap


extract_tsne_rectangle <- function(dat, x_min, x_max, y_min, y_max, antigen, marker, polyDegree = 0 , facet, title) {
  facet <- as.symbol(facet)
  marker <- as.symbol(marker)
  # `marker` argument can be an actual marker (eg. "IFNg") or "degree"
  allData <- dat
  dat <- subset(dat, stim==antigen)

  m <- unique(dat[,.(ptid,neut_status,controller_status)])

  totals <- dat[, length(x), by=ptid]
  setnames(totals, "V1", "total")

  dat <- subset(dat, x>x_min&x<x_max&y>y_min&y<y_max)

  title <- paste(title, "," , antigen, sep=" ")

  if(polyDegree > 0){
    dat <- subset(dat, degree > polyDegree)
    title <- paste(title, ", degree", polyDegree+1, "and higher",sep=" ")
  }

  plotCall <- substitute(ggplot(allData, aes(x = x, y = y)) +
                           geom_point(aes(colour = mkr), size = 1, alpha = 0.8) +
                           geom_point(data=dat, colour="red", size =2 , alpha = 0.8) +
                           geom_rect(xmin=x_min, xmax=x_max, ymin=y_min, ymax=y_max,
                                     color="red", alpha=0) +
                           labs(title = title) +
                           facet_grid(fct ~ .) +
                           scale_color_gradientn(colours = rev(brewer.pal(11, "Spectral")))
                         ,
                         list(fct = facet, mkr = marker)
  )

  p1 <- eval(plotCall)


  non_h <- subset(dat, neut_status =="non")
  neut_h <- subset(dat, neut_status =="neut")

  non_s <- non_h[, length(x), by=ptid]
  setnames(non_s, "V1", "in_hull")
  neut_s <- neut_h[, length(x), by=ptid]
  setnames(neut_s, "V1", "in_hull")

  subs <- rbind(non_s, neut_s)

  m <- merge(m, subs, by="ptid")
  m <- merge(m, totals, by="ptid")
  m[,proportion := in_hull/total]
  m[,log_proportion := log10(proportion)]

  p <- t.test(log_proportion ~ neut_status, alternative="greater", data=m)

  p2 <- ggplot(m, aes(x=neut_status, y=log(in_hull))) +
    geom_boxplot(outlier.size = 3) +
    geom_point(aes(col=ptid), size=3) +
    theme(legend.position="none") +
    labs(title=paste0("selected proportion, p=", round(p$p.value, 3)), x="")

  p3 <- grid.rect(gp=gpar(col="white"))
  pl <- grid.arrange(p1, arrangeGrob(p2,p3, nrow=1), heights=c(1,0.5))
  print(pl)
  print(m)
  return(tmp)
}



