% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-pathway-network.R
\name{plotPathwayNetwork}
\alias{plotPathwayNetwork}
\title{Plot a pathway network}
\usage{
plotPathwayNetwork(
  PAResults,
  genesets,
  selectedPathways = NULL,
  statistic = "pFDR",
  mode = c("continuous", "discrete"),
  pThreshold = 0.05,
  useFDR = TRUE,
  edgeThreshold = 0.5,
  statLimit = 4,
  discreteColors = NULL,
  continuousScaleFunc = NULL,
  NAColor = "#ffffff",
  borderColor = "#333333",
  nodeSizeFnc = function(id) length(genesets[[id]])^0.75,
  borderWidthFnc = function(id) 1,
  edgeWidthFnc = function(from, to) 1,
  styleFile = system.file(package = "RCPA", "extdata", "pieStyle.js"),
  file = tempfile(fileext = ".html")
)
}
\arguments{
\item{PAResults}{A named list of data frame of Pathway analysis results.
The columns of each data frame should be at least ID, name, p.value and pFDR.
An optional column "color" can be used to specify the color of the nodes.
If the column "color" is not specified, the color of the nodes will be determined by the mode and the statistic.}

\item{genesets}{A genesets object that is obtained from getGeneSets function.}

\item{selectedPathways}{A vector of pathway IDs to be included in the plot.
If it is NULL, all pathways from genesets will be included.}

\item{statistic}{A character value of the statistic to use for the network.
The statistic should be one of the columns of the data frame in the results.}

\item{mode}{A character value of the mode to use to color the nodes.
The mode should be one of "discrete", "continuous".
If the mode is "discrete", the color of the nodes are determined by whether the p-value is significant or not.
If the mode is "continuous", the color of the nodes are proportional to the statistic.}

\item{pThreshold}{A numeric value of p-value threshold.}

\item{useFDR}{A logical value indicating whether to use FDR or not.
This parameter is independent of the pThreshold.}

\item{edgeThreshold}{A numeric from 0 to 1 indicating the threshold to draw edges.
edgeThreshold of 0.1 means that edges are drawn if the number of genes in common is greater than 1\% of the smaller gene set.}

\item{statLimit}{A numeric value of the maximum absolute value of the statistic.
If statistic is p.value or pFDR, this parameter is the limit of -log10(p-value).}

\item{discreteColors}{A character vector of colors to use for the discrete mode.
The length of the vector must be the same as the number of results.}

\item{continuousScaleFunc}{A function that takes a numeric value from -1 to 1 and returns a color.}

\item{NAColor}{A character value of the color to use for NA values.}

\item{borderColor}{A character value of the color to use for the border of the nodes.}

\item{nodeSizeFnc}{A function that takes a character value of the ID of the node and returns a numeric value of the size of the node.}

\item{borderWidthFnc}{A function that takes a character value of the ID of the node and returns a numeric value of the width of the border of the node.}

\item{edgeWidthFnc}{A function that takes a character value of the ID of the from node and a character value of the ID of the to node and returns a numeric value of the width of the edge.}

\item{styleFile}{A character value of the path to the style file.
If NULL, the default style file will be used, which is located at system.file(package="RCPA", "extdata", "pieStyle.js")}

\item{file}{A character value of the path to the html file to be created.}
}
\value{
A character value of the html content of the plot.
}
\description{
This function plots a pathway network.
This function needs an interactive environment with browser view support to work.
}
\details{
The function will plot a pathway network using the results of pathway analysis.
The nodes of the network are the pathways and the edges are the pathways that have at least a certain number of genes in common defined by the edgeThreshold.
The size of the nodes are proportional to the number of genes in the pathway.
The color of the nodes are proportional to the statistic used if the mode is "continuous".
If the mode is "discrete", the color of the nodes are determined by whether the p-value is significant or not.
The width of the edges are proportional to the number of genes in common.
}
\examples{
\donttest{

library(RCPA)

affyFgseaResult <- loadData("affyFgseaResult")
agilFgseaResult <- loadData("agilFgseaResult")
RNASeqFgseaResult <- loadData("RNASeqFgseaResult")
metaPAResult <- loadData("metaPAResult")
genesets <- loadData("genesets")

PAResults <- list(
    "Affymetrix - GSE5281" = affyFgseaResult,
    "Agilent - GSE61196" = agilFgseaResult,
    "RNASeq - GSE153873" = RNASeqFgseaResult,
    "Meta-analysis" = metaPAResult
)

genesetsToPlot <- metaPAResult$ID[order(metaPAResult$pFDR)][1:30]

pltHtml <- RCPA::plotPathwayNetwork(
    PAResults,
    genesets = genesets,
    selectedPathways = genesetsToPlot,
    edgeThreshold = 0.75,
    mode = "continuous",
    statistic = "normalizedScore"
)

}

}
