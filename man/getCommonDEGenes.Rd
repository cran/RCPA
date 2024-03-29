% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/common-de.R
\name{getCommonDEGenes}
\alias{getCommonDEGenes}
\title{Get common significant DE genes from multiple DE Analysis results}
\usage{
getCommonDEGenes(
  DEResults,
  pThreshold = 0.05,
  useFDR = TRUE,
  stat = "logFC",
  statThreshold = 0
)
}
\arguments{
\item{DEResults}{A list of data frames with the results of DE analysis.}

\item{pThreshold}{The p-value threshold to determine if a gene is differentially expressed.}

\item{useFDR}{Use the FDR adjusted p-value instead of the nominal p-value.}

\item{stat}{The additional statistics column to use for filtering differentially expressed genes.}

\item{statThreshold}{The absolute value of the statistic threshold to use for filtering differentially expressed genes.
Default is 0, which means no filtering.}
}
\value{
A data frame wtih three columns: ID (Entrez IDs), Symbol and Description
}
\description{
Get a list of common significant DE genes from multiple DE Analysis results.
}
\examples{
\donttest{
library(RCPA)
library(SummarizedExperiment)

affyDEExperiment <- loadData("affyDEExperiment")
agilDEExperiment <- loadData("agilDEExperiment")
RNASeqDEExperiment <- loadData("RNASeqDEExperiment")

DEResults <- list(
    "Affymetrix - GSE5281" = rowData(affyDEExperiment),
    "Agilent - GSE61196"   = rowData(agilDEExperiment),
    "RNASeq - GSE153873"   = rowData(RNASeqDEExperiment)
)

commonDEGenes <- RCPA::getCommonDEGenes(DEResults)

print(head(commonDEGenes))
}
}
