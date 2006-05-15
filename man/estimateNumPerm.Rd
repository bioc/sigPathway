\name{estimateNumPerm}
\alias{estimateNumPerm}
\title{Compute the Number of Unique Permutations for a Phenotype Vector}
\description{
  Computes the number of unique permutations based on a vector of
  phenotypes and the number of groups.
}
\usage{
estimateNumPerm(phenotype, ngroups)
}
\arguments{
  \item{phenotype}{a numeric (or character if \code{ngroups} >= 2)
    vector indicating the phenotype}
  \item{ngroups}{an integer indicating the number of groups in the phenotype}
}
\details{
  This function calculates the number of unique permutations based on
  the given phenotype and the number of groups present in the
  phenotype.  This function is used internally in \code{sigPathway} and
  attempts to avoid numeric overflow associated with multiplying out
  large factorials.
}
\value{
  A numeric with length 1.
}
\author{Weil Lai}
\examples{
## One group: continuous observations
ptype1 <- c(24,25,17,26,25,16,14,17,12,15,19,20)
print(estimateNumPerm(ptype1, 1))

## Two groups
ptype2 <- c(0,1,1,0,1,0,1)
print(estimateNumPerm(ptype2, 2))

## Three groups
ptype3a <- c(2,0,1,2,0,1,2,0,0,1,1,2)
print(estimateNumPerm(ptype3a, 3))

ptype3b <- c("Banana","Apple","Lemon","Lemon","Lemon",
             "Apple","Lemon","Banana","Banana")
print(estimateNumPerm(ptype3b, 3))
}
\keyword{array}
\keyword{htest}
