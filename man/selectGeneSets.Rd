\name{selectGeneSets}
\alias{selectGeneSets}
\title{Select gene sets to be analyzed in pathway analysis}
\description{
  Selects gene sets to be analyzed in pathway analysis based on minimum
  and maximum number of probe sets to consider per pathway.
}
\usage{
selectGeneSets(G, probeID, minNPS = 20, maxNPS = 500)
}
\arguments{
  \item{G}{a list containing the source, title, and probe sets
    associated with each curated pathway}
  \item{probeID}{a character vector containing the names of probe sets
    associated with a matrix of expression values}
  \item{minNPS}{an integer specifying the minimum number of probe sets
    in \code{probeID} that should be in a gene set}
  \item{maxNPS}{an integer specifying the maximum number of probe sets
    in \code{probeID} that should be in a gene set}
}
\details{
  This function selects the appropriate pathways from a large, curated
  list based on the minimum and maximum number of probe sets that should
  be considered in a gene set.  It creates three vectors: \code{nprobesV}
  and \code{indexV} representing a sparse indicator matrix and
  \code{indGused} indicating which gene sets were selected from
  \code{G}.
}
\value{
  A list containing
  \item{nprobesV}{an integer vector indicating the number of probe sets
    in \code{probeID} that is in each selected gene set}
  \item{indexV}{an integer vector containing positions for each 1s in
    the sparse indicator matrix}
  \item{indGused}{an integer vector indicating which pathways in
    \code{G} were chosen}
}
\references{
  Tian L., Greenberg S.A., Kong S.W., Altschuler J., Kohane I.S., Park
  P.J. (2005)  Discovering statistically significant pathways in
  expression profiling studies.  \emph{Proceedings of the National
    Academy of Sciences of the USA}, \bold{102}, 13544-9.

  \url{http://www.pnas.org/cgi/doi/10.1073/pnas.0506577102}
}
\note{See the help page for \code{calculate.NTk} or \code{calculate.NEk}
  for example code that uses \code{getPathwayStatistics}}
\author{Lu Tian, Peter Park, and Weil Lai}
\keyword{array}
\keyword{htest}
