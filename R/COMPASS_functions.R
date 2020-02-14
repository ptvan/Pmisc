#' calculate cell proportions from the COMPASSResult container (CR[[antigen]]$data$count_s)
#' !!! IMPORTANT: assumes data comes from time-course experiment and therefore has
#' joint_id (<ptid>_<time>)
#'
#' @param CR A COMPASSResult object
#' @return table cell proportions, null category is included
#' @import COMPASS
#' @import data.table
#'
extract_COMPASS_cell_proportions <- function(CR){
 if(!is.list(CR) || class(CR[[1]]) != "COMPASSResult") {
    stop("this function requires a list containing at least 1 COMPASSResult !!!")
  }

  antigens <- names(CR)

  # get cytokines (last column is "Counts", exclude it)
  cytokines <- colnames(CR[[1]]$data$categories)
  cytokines <- cytokines[-length(cytokines)]

  # create the big data.frame to be returned later
  big <- data.frame(matrix(vector(),0, length(cytokines)+5))
  colnames(big) <- c("ptid","time",cytokines,"unstim_proportion" ,"stim_proportion","antigen")

  # for each antigen in the COMPASSResult structure...
  for (i in 1:length(antigens)){
    antigen <- antigens[i]
    cat(antigen, "\n")
    p <- as.data.frame(CR[[antigen]]$data$n_s)
    p_u <- as.data.frame(CR[[antigen]]$data$n_u)

    # get total cell counts
    totals <- CR[[antigen]]$data$counts_s
    totals_u <- CR[[antigen]]$data$counts_u

    p_new <- data.frame(matrix(vector(),0, length(cytokines)+5))
    colnames(p_new) <- c("ptid","time",cytokines,"unstim_proportion", "stim_proportion","antigen")

    k <- 0
    for (i in 1:nrow(p)){
      rown <- rownames(p)[i]
      # split jointID into ptid and time
      # eg. "1_0", ptid == 1, time == 0

      ptid <- unlist(strsplit(rown, "_"))[1]
      time <- unlist(strsplit(rown, "_"))[2]

      # calculate proportions, keeping the antigen combinations
      for (j in 1:ncol(p)){ # NOTE: this assumes count_s & count_u have same dimensions
        coln <- colnames(p)[j]
        x <- unlist(strsplit(coln, "&"))
        vec <- rep("+", length(cytokines))
        vec[grep("\\!", x)] = "-"
        vec <- c(ptid, time, vec, p[i,j] / totals[i], p_u[i,j] / totals_u[i], antigen)
        k <- k+1
        # cat("row", k, ":  ", vec, "\n")
        p_new[k,] <- vec
      }
    }
    # append to big matrix
    big <- rbind(big,p_new)
  }

  big$time <- as.factor(big$time)
  big$unstim_proportion <- as.numeric(big$unstim_proportion)
  big$stim_proportion <- as.numeric(big$stim_proportion)
  big <- as.data.table(big)
  return(big)
}


#' calculate cell proportions from the COMPASSResult container (CR[[antigen]]$data$count_s)
#' !!! IMPORTANT: assumes data comes from time-course experiment and therefore has
#' joint_id (<ptid>_<time>)
#'
#' This function calculates percentages EXCLUDING THE NULL CATEGORY
#' eg. you have 10000 cells
#' 5 cells express IL2 only, 3 cells express IFNg only, 2 cells express both (10 cells express any)
#' the percentages returned by this fx are IL2 = 5/10 = 0.5, IFNg = 3/10 = 0.3, IL2+IFNg = 2/10 = 0.2
#' this is different than when the calculation is over ALL cells (10000 in this example)
#'
#' @param CR A COMPASSResult object
#' @return table cell proportions, null category is NOT included
#'
#' @import COMPASS
extract_COMPASS_cell_proportions_nonull <- function(CR){
  if(!is.list(CR) || class(CR[[1]]) != "COMPASSResult") {
    stop("this function requires a list containing at least 1 COMPASSResult !!!")
  }

  antigens <- names(CR)

  # get cytokines (last column is "Counts", exclude it)
  cytokines <- colnames(CR[[1]]$data$categories)
  cytokines <- cytokines[-length(cytokines)]

  # create the big data.frame to be returned later
  big <- data.frame(matrix(vector(),0, length(cytokines)+6))
  colnames(big) <- c("ptid","time",cytokines,"unstim_proportion" ,"stim_proportion","antigen")

  for (i in 1:length(antigens)){
    antigen <- antigens[i]
    cat(antigen, "\n")

    # get cell counts, REMOVE NULL CATEGORY
    p <- as.data.frame(CR[[antigen]]$data$n_s)
    p <- p[,1:ncol(p)-1]
    p_u <- as.data.frame(CR[[antigen]]$data$n_u)
    p_u <- p_u[,1:ncol(p_u)-1]

    # RECALCULATE TOTAL
    totals <- apply(p, 1, sum)
    totals_u <- apply(p_u, 1, sum)

    p_new <- data.frame(matrix(vector(),0, length(cytokines)+5))
    colnames(p_new) <- c("ptid","time",cytokines,"unstim_proportion", "stim_proportion","antigen")
    k <- 0
    for (i in 1:nrow(p)){
      rown <- rownames(p)[i]
      # since this is time-series, need to split jointID into ptid and time
      # eg. "1_0", ptid == 1, time == 0

      ptid <- unlist(strsplit(rown, "_"))[1]
      time <- unlist(strsplit(rown, "_"))[2]
      # calculate proportions, keeping the antigen combinations
      for (j in 1:ncol(p)){ # NOTE: this assumes count_s & count_u have same dimensions
        coln <- colnames(p)[j]
        x <- unlist(strsplit(coln, "&"))
        vec <- rep("+", length(cytokines))
        vec[grep("\\!", x)] = "-"
        vec <- c(ptid, time, vec, p[i,j] / totals[i], p_u[i,j] / totals_u[i], antigen)
        k <- k+1
        # cat("row", k, ":  ", vec, "\n")
        p_new[k,] <- vec
      }
    }
    # append to big matrix
    big <- rbind(big,p_new)
  }

  big$time <- as.factor(big$time)
  big$unstim_proportion <- as.numeric(big$unstim_proportion)
  big$stim_proportion <- as.numeric(big$stim_proportion)
  big <- as.data.table(big)
  return(big)
}


#' extract cell counts from the COMPASSResult container (CR[[antigen]]$data$count_s)
#' !!! IMPORTANT: assumes data comes from time-course experiment and therefore has
#' joint_id (<ptid>_<time>)
#'
#' @param CR A COMPASSResult object
#' @return table of cell counts, null category is EXCLUDED
#' @import COMPASS
#'
COMPASS_cell_counts_nonull <- function(CR){

  if(!is.list(CR) || class(CR[[1]]) != "COMPASSResult") {
    stop("this function requires a list containing at least 1 COMPASSResult !!!")
  }

  antigens <- names(CR)

  # create the big data.frame to be returned later
  big <- data.frame(matrix(vector(),0, length(cytokines)+7))
  colnames(big) <- c("ptid","time",cytokines, "degree" ,"poly", "unstim_count" ,"stim_count","antigen")

  # get cytokines (last column is "Counts", exclude it)
  cytokines <- colnames(CR[[1]]$data$categories)
  cytokines <- cytokines[-length(cytokines)]

  for (i in 1:length(antigens)){
    antigen <- antigens[i]
    cat(antigen, "\n")

    # get cell counts, REMOVE NULL CATEGORY
    p <- as.data.frame(CR[[antigen]]$data$n_s)
    p <- p[,1:ncol(p)-1]
    p_u <- as.data.frame(CR[[antigen]]$data$n_u)
    p_u <- p_u[,1:ncol(p_u)-1]

    p_new <- data.frame(matrix(vector(),0, length(cytokines)+7))
    colnames(p_new) <- c("ptid","time",cytokines, "degree", "poly", "unstim_count", "stim_count","antigen")
    k <- 0

    for (i in 1:nrow(p)){
      rown <- rownames(p)[i]
      # since this is time-series, need to split jointID into ptid and time
      # eg. "1_0", ptid == 1, time == 0

      ptid <- unlist(strsplit(rown, "_"))[1]
      time <- unlist(strsplit(rown, "_"))[2]
      # calculate proportions, keeping the antigen combinations
      for (j in 1:ncol(p)){ # NOTE: this assumes count_s & count_u have same dimensions

        coln <- colnames(p)[j]
        x <- unlist(strsplit(coln, "&"))
        vec <- rep("+", length(cytokines))
        vec[grep("\\!", x)] = "-"
        poly <- paste(cytokines[which(vec=="+")], collapse=" ")
        degree <- length(which(vec=="+"))
        vec <- c(ptid, time, vec, degree, poly, p[i,j], p_u[i,j], antigen)
        k <- k+1
        # cat("row", k, ":  ", vec, "\n")
        p_new[k,] <- vec
      }
    }
    # append to big matrix
    big <- rbind(big,p_new)
  }

  big$time <- as.factor(big$time)
  big$unstim_count <- as.numeric(big$unstim_count)
  big$stim_count <- as.numeric(big$stim_count)
  big <- as.data.table(big)
  return(big)

}


#' extract cell counts from the COMPASSResult container (CR[[antigen]]$data$count_s)
#' !!! IMPORTANT: assumes data comes from time-course experiment and therefore has
#' joint_id (<ptid>_<time>)
#'
#' @param CR A COMPASSResult object
#' @return table of cell counts, including null category
#' @import COMPASS
#'
COMPASS_cell_counts <- function(CR){

  if(!is.list(CR) || class(CR[[1]]) != "COMPASSResult") {
    stop("this function requires a list containing at least 1 COMPASSResult !!!")
  }

  antigens <- names(CR)

  # create the big data.frame to be returned later
  big <- data.frame(matrix(vector(),0, length(cytokines)+7))
  colnames(big) <- c("ptid","time",cytokines, "degree" ,"poly", "unstim_count" ,"stim_count","antigen")

  # get cytokines (last column is "Counts", exclude it)
  cytokines <- colnames(CR[[1]]$data$categories)
  cytokines <- cytokines[-length(cytokines)]

  for (i in 1:length(antigens)){
    antigen <- antigens[i]
    cat(antigen, "\n")

    # get cell counts, **including** null category
    p <- as.data.frame(CR[[antigen]]$data$n_s)
    p_u <- as.data.frame(CR[[antigen]]$data$n_u)

    p_new <- data.frame(matrix(vector(),0, length(cytokines)+7))
    colnames(p_new) <- c("ptid","time",cytokines, "degree", "poly", "unstim_count", "stim_count","antigen")
    k <- 0

    for (i in 1:nrow(p)){
      rown <- rownames(p)[i]
      # since this is time-series, need to split jointID into ptid and time
      # eg. "1_0", ptid == 1, time == 0

      ptid <- unlist(strsplit(rown, "_"))[1]
      time <- unlist(strsplit(rown, "_"))[2]
      # calculate proportions, keeping the antigen combinations
      for (j in 1:ncol(p)){ # NOTE: this assumes count_s & count_u have same dimensions

        coln <- colnames(p)[j]
        x <- unlist(strsplit(coln, "&"))
        vec <- rep("+", length(cytokines))
        vec[grep("\\!", x)] = "-"
        poly <- paste(cytokines[which(vec=="+")], collapse=" ")
        degree <- length(which(vec=="+"))
        vec <- c(ptid, time, vec, degree, poly, p[i,j], p_u[i,j], antigen)
        k <- k+1
        # cat("row", k, ":  ", vec, "\n")
        p_new[k,] <- vec
      }
    }
    # append to big matrix
    big <- rbind(big,p_new)
  }

  big$time <- as.factor(big$time)
  big$unstim_count <- as.numeric(big$unstim_count)
  big$stim_count <- as.numeric(big$stim_count)
  big <- as.data.table(big)
  return(big)

}

#' extract polyfunctionality degree from a table of melted counts
#' column names must be Boolean expression of cytokines (eg. "IL2&!IFNg&TNFa"...)
#'
#' @param melted_counts A melted data.table of cell counts
#' @return table of cell counts
#' @import data.table
extract_degree_from_melted_counts <- function(melted_counts){
  melted_counts[,.(degree=sum(!grepl("!",str_split_fixed(poly,"&",4)))),.(1:nrow(melted_counts),jointID,poly,count)]}


#' extract polyfunctionality degree from a COMPASSResult object
#'
#' @param CR A COMPASSResult object
#' @return table of polyfunctionality degree
#' @import COMPASS
extract_degree_from_COMPASS <- function(CR){
  if(!is.list(CR) || class(CR[[1]]) != "COMPASSResult") {
    stop("this function requires a list containing at least 1 COMPASSResult !!!")
  }

  nu <- CR$data$n_u
  ns <- CR$data$n_s
  ns <- melt(ns)
  nu <- melt(nu)
  colnames(nu) <- c("jointID","poly","count")
  colnames(ns) <- c("jointID","poly","count")
  nu <- data.table(nu)
  ns <- data.table(ns)
  nu <- extractDegree(nu)
  ns <- extractDegree(ns)
  rbind(cbind(nu,isStim="nonStim"),cbind(ns,isStim="Stim"))
}
