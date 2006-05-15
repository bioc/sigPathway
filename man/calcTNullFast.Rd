\name{calcTNullFast}
\alias{calcTNullFast}
\title{Compute Null T Distribution for Each Gene}
\description{
  Computes a null t distribution for each gene by permuting the phenotypes.
}
\usage{
calcTNullFast(tab, phenotype, nsim, ngroups = 2, allphenotypes = FALSE)
}
\arguments{
  \item{tab}{a numeric matrix of expression values, with the rows and
    columns representing probe sets and sample arrays, respectively}
  \item{phenotype}{a numeric (or character if \code{ngroups} >= 2)
    vector indicating the phenotype}
  \item{nsim}{an integer indicating the number of permutations to use}
  \item{ngroups}{an integer indicating the number of groups in the expression
    matrix}
  \item{allphenotypes}{a boolean indicating whether the function should
    consider all possible permutations of the phenotype, including the
    original, non-permuted phenotype}
}
\details{
  Similar to \code{calcTStatFast} but calculates t-statistics over permuted
  phenotypes.  If \code{allphenotypes == FALSE}, then any permutation
  that has a permuted phenotype equal to the original phenotype will be
  repermuted.  For example, all the possible permutations for
  \code{phenotype == c(0,0,1,1)} are \code{c(0,0,1,1)},
  \code{c(0,1,0,1)}, \code{c(1,0,1,0)}, \code{c(1,0,0,1)},
  \code{c(0,1,1,0)}, and \code{c(1,1,0,0)}.  If \code{allphenotypes ==
    FALSE}, then the results will not include values from the
  \code{c(0,0,1,1)} case.
  
  The help file of \code{calcTStatFast} has more details on the
  different statistics one can calculate based on the
  value specified for \code{ngroups}.
}
\value{
  A matrix with \code{nsim} rows and \code{nrow(tab)} columns.
}
\author{Weil Lai}
\keyword{array}
\keyword{htest}
