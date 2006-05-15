\name{rankPathways.NGSk}
\alias{rankPathways.NGSk}
\title{Summarizes Top Pathways from One of the Pathway Analyses}
\description{
  Summarizes top pathways from one of the pathway analyses (i.e.,
  \code{calculate.NTk}, \code{calculate.NEk}, \code{calculate.NGSk}, or
  \code{calculate.GSEA})
}
\usage{
rankPathways.NGSk(res.NGSk, G, gsList, methodName = "NGSk",
                  npath = 25)
}
\arguments{
  \item{res.NGSk}{a list from the output of \code{calculate.NGSk},
    \code{calculate.NTk}, \code{calculate.NEk}, or
    \code{calculate.GSEA}}
  \item{G}{a list containing the source, title, and probe sets
    associated with each curated pathway}
  \item{gsList}{a list containing three vectors from the output of
    the \code{selectGeneSets} function}
  \item{methodName}{a character vector of length 1 giving the name of
    the pathway analysis used in making \code{res.NGSk}}
  \item{npath}{an integer indicating the number of top gene sets to
    consider when ranking the top pathways}
}
\details{
  This function ranks the statistics given in \code{res.NGSk} and
  summarizes the top gene sets in a tabular format similar to Table 2 in
  Tian et al. (2005)
}
\value{
  A data frame showing the pathways' indices in \code{G}, gene set
  category, pathway title, set size, \code{res.NGSk}'s statistics,
  the corresponding q-values, and the numerical ranks for the top gene sets.
}
\references{
  Tian L., Greenberg S.A., Kong S.W., Altschuler J., Kohane I.S., Park
  P.J. (2005)  Discovering statistically significant pathways in
  expression profiling studies.  \emph{Proceedings of the National
    Academy of Sciences of the USA}, \bold{102}, 13544-9.

  \url{http://www.pnas.org/cgi/doi/10.1073/pnas.0506577102}
}
\note{See the help page for \code{calculate.NGSk} for example code that
  uses \code{rankPathways.NGSk}}
\author{Lu Tian, Peter Park, and Weil Lai}
\keyword{array}
\keyword{htest}
