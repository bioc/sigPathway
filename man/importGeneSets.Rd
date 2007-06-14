\name{importGeneSets}
\alias{importGeneSets}
\alias{gmtToG}
\alias{gmxToG}
\alias{grpToG}
\alias{xmlToG}
\title{Import gene sets stored in GMT, GMX, GRP, and XML file formats}
\description{
  Imports gene sets stored in GMT, GMX, GRP, and XML file formats and
  converts them to sigPathway's preferred format.
}
\usage{
importGeneSets(fileNames, verbose = TRUE)
gmtToG(fileNames, verbose = TRUE)
gmxToG(fileNames, verbose = TRUE)
grpToG(fileNames, verbose = TRUE)
xmlToG(fileNames, verbose = TRUE)
}
\arguments{
  \item{fileNames}{a character vector specifying the file(s) containing the gene
    sets of interest}
  \item{verbose}{a boolean to indicate whether to print debugging messages
    to the R console}
}
\details{
  These functions read in gene sets stored in GMT, GMX, GRP, and XML
  file formats and converts them to a list format that sigPathway can
  use.  Redundant gene IDs in each gene set are removed during
  conversion.  The \code{importGeneSets} function can read in GMT, GMX,
  GRP, and XML files in one pass.  The \code{gmtToG}, \code{gmxToG},
  \code{grpToG}, and \code{xmlToG} functions are specific to reading in
  their respective file formats.
}
\value{
  A list containing sublists representing each imported gene set.  The
  vignette contains more details about the list structure.
}
\note{These functions do not check whether the files are in the correct
  format and will give spurious output when given files in the wrong
  format.  The \code{xmlToG} function requires the \code{XML} package, which
  is available on CRAN.  The \code{xmlToG} function also requires XML
  files to be formatted based on the MSigDB Document Type Definition.
}
\references{
  \url{http://www.broad.mit.edu/cancer/software/gsea/wiki/index.php/Data_formats}
}
\author{Weil Lai}
\keyword{file}
\keyword{datagen}
