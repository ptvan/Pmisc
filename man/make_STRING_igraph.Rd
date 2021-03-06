% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_STRING_igraph.R
\name{make_STRING_igraph}
\alias{make_STRING_igraph}
\title{this function takes an input matrix having one column named "gene" 
containing genes of interest and a stringDB object
it returns an igraph graph of the STRING network containing those genes
the node attributes are "name" (STRING id) and "geneName" (original gene names from input matrix)
genes that do not have STRING annotation are omitted
Phu T. Van, FHCRC 2017, w/ substantial help and input from C.Murie & V.Voillet}
\usage{
make_STRING_igraph(inputMatrix, STRINGdbObj)
}
\arguments{
\item{inputMatrix}{a matrix of gene expression values}

\item{STRINGdbObj}{a STRINGdbObj generated by connecting via STRINGdb}
}
\value{
table of genes' mean expression values for each group
}
\description{
this function takes an input matrix having one column named "gene" 
containing genes of interest and a stringDB object
it returns an igraph graph of the STRING network containing those genes
the node attributes are "name" (STRING id) and "geneName" (original gene names from input matrix)
genes that do not have STRING annotation are omitted
Phu T. Van, FHCRC 2017, w/ substantial help and input from C.Murie & V.Voillet
}
\examples{
\dontrun{
# DEGout <- topTable(fitBayes, number=nrow(vDat), coef="pttype", sort="P")
DEGout$gene <- rownames(DEGout)
string_db <- STRINGdb$new(version="10"
                          , species=9606
                          , score_threshold=0
                          , input_directory="/STRINGdb/HomoSapiens/")
n <- make_STRING_igraph(DEGout, string_db)
plot(n, vertex.label=V(n)$geneName)
}
}
