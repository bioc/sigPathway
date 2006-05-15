\name{rankPathways}
\alias{rankPathways}
\title{Summarizes Top Pathways from Pathway Analyses}
\description{
  Summarizes top pathways from pathway analyses.
}
\usage{
rankPathways(res.A, res.B, G, tab, phenotype, gsList, ngroups,
             methodNames = NULL, npath = 25, allpathways = FALSE)
}
\arguments{
  \item{res.A}{a list from the output of \code{calculate.NTk} or
    \code{calculate.NEk}}
  \item{res.B}{a list from the output of \code{calculate.NTk} or
    \code{calculate.NEk}}
  \item{G}{a list containing the source, title, and probe sets
    associated with each curated pathway}
  \item{tab}{a numeric matrix of expression values, with the rows and
    columns representing probe sets and sample arrays, respectively}
  \item{phenotype}{a numeric (or character if \code{ngroups} >= 2)
    vector indicating the phenotype}
  \item{gsList}{a list containing three vectors from the output of
    the \code{selectGeneSets} function}
  \item{ngroups}{an integer indicating the number of groups in the
    matrix}
  \item{methodNames}{a character vector of length 2 giving the names for
    \code{res.A} and \code{res.B}}
  \item{npath}{an integer indicating the number of top gene sets to
    consider from each statistic when ranking the top pathways}
  \item{allpathways}{a boolean to indicate whether to include the top
    \code{npath} pathways from each statistic or just consider the top
    \code{npath} pathways (sorted by the sum of ranks of both statistics) when
    generating the summary table}
}
\details{
  This function ranks together the statistics given in \code{res.A} and
  \code{res.B} and summarizes the top gene sets in a tabular format similar to
  Table 2 in Tian et al. (2005)
}
\value{
  A data frame showing the pathways' indices in \code{G}, gene set
  category, pathway title, set size, \code{res.A}'s statistics,
  \code{res.B}'s statistics, the corresponding q-values, and the ranks for the
  top gene sets.
}
\references{
  Tian L., Greenberg S.A., Kong S.W., Altschuler J., Kohane I.S., Park
  P.J. (2005)  Discovering statistically significant pathways in
  expression profiling studies.  \emph{Proceedings of the National
    Academy of Sciences of the USA}, \bold{102}, 13544-9.

  \url{http://www.pnas.org/cgi/doi/10.1073/pnas.0506577102}
}
\note{See the help page for \code{calculate.NTk} or \code{calculate.NEk}
  for example code that uses \code{rankPathways}}
\author{Lu Tian, Peter Park, and Weil Lai}
\keyword{array}
\keyword{htest}
