\name{getPathwayStatistics.NGSk}
\alias{getPathwayStatistics.NGSk}
\title{Give the statistics for the probe sets in a pathway}
\description{
  Gives the statistics for the probe sets associated with a pathway.
}
\usage{
getPathwayStatistics.NGSk(statV, probeID, G, index,
                          keepUnknownProbes = FALSE, annotpkg = NULL)
}
\arguments{
  \item{statV}{a numeric vector of test statistic (not p-values) for
    each individual probe/gene}
  \item{probeID}{a character vector containing the names of probe sets
    associated with a matrix of expression values}
  \item{G}{a list containing the source, title, and probe sets
    associated with each curated pathway}
  \item{index}{an integer vector specifying the pathway(s) to summarize
    in \code{G}}
  \item{keepUnknownProbes}{a boolean indicating whether to keep the
    names of probe sets not represented in tab in the summary data
    frame}
  \item{annotpkg}{a character vector specifying the name of the
    BioConductor annotation package to use to fetch accession numbers,
    Entrez Gene IDs, gene name, and gene symbols}
}
\details{
  This function gives the test statistic for each probe in the pathway
  as indicated in \code{G[[index]]}.
}
\value{
  A list containing data frames (1 per pathway) with the probes' name
  and the corresponding test statistic.

  If a valid \code{annotpkg} is specified, the probes' accession
  numbers, Entrez Gene IDs, gene name, and gene symbols are also returned.
}
\note{See the help page for \code{calculate.NGSk} for example code that
  uses \code{getPathwayStatistics.NGSk}}
\author{Weil Lai}
\keyword{array}
\keyword{htest}
