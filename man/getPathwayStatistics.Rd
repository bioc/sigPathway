\name{getPathwayStatistics}
\alias{getPathwayStatistics}
\title{Give the statistics for the probe sets in a pathway}
\description{
  Gives the statistics for the probe sets associated with a pathway.
}
\usage{
getPathwayStatistics(tab, phenotype, G, index, ngroups = 2,
                     statList = NULL, keepUnknownProbes = FALSE,
                     annotpkg = NULL)
}
\arguments{
  \item{tab}{a numeric matrix of expression values, with the rows and
    columns representing probe sets and sample arrays, respectively}
  \item{phenotype}{a numeric (or character if \code{ngroups} >= 2)
    vector indicating the phenotype}
  \item{G}{a list containing the source, title, and probe sets
    associated with each curated pathway}
  \item{index}{an integer vector specifying the pathway(s) to summarize
    in \code{G}}
  \item{ngroups}{an integer indicating the number of groups in the expression
    matrix}
  \item{statList}{a list containing results from \code{calcTStatFast}}
  \item{keepUnknownProbes}{a boolean indicating whether to keep the
    names of probe sets not represented in tab in the summary data
    frame}
  \item{annotpkg}{a character vector specifying the name of the
    BioConductor annotation package to use to fetch accession numbers,
    Entrez Gene IDs, gene name, and gene symbols}
}
\details{
  This function gives the mean, standard deviation, and test statistic
  for each probe in the pathway as indicated in \code{G[[index]]}.
}
\value{
  A list containing data frames (1 per pathway) with the probes' name,
  mean, standard deviation, the test statistic (e.g., t-test), and the
  corresponding unadjusted p-value.

  If \code{ngroups} = 1, the Pearson correlation coefficient is also
  returned.

  If a valid \code{annotpkg} is specified, the probes' accession
  numbers, Entrez Gene IDs, gene name, and gene symbols are also
  returned.  This option only works if the probes in the gene set list
  \code{G} are manufacturer IDs corresponding to those used in making
  the BioConductor annotation package.
}
\note{See the help page of \code{calculate.NTk} or \code{calculate.NEk}
  for example code that uses \code{getPathwayStatistics}}
\author{Weil Lai}
\keyword{array}
\keyword{htest}
