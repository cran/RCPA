% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotation.R
\name{getEntrezAnnotation}
\alias{getEntrezAnnotation}
\title{Get Entrez annotation}
\usage{
getEntrezAnnotation(entrezIds)
}
\arguments{
\item{entrezIds}{A vector of Entrez IDs.}
}
\value{
A data frame with Entrez annotation. The columns are ID (Entrez ID), Symbol, Description, OtherDesignations, OtherAliases and Chromosome.
}
\description{
This function gets Entrez annotation.
}
\examples{

library(RCPA)
geneAnno <- getEntrezAnnotation(c("77267466", "77267467"))

}
