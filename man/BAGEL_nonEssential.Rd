\name{BAGEL_nonEssential}
\alias{BAGEL_nonEssential}
\docType{data}
\title{
Reference set of non essential genes
}
\description{
A list of reference non essential genes assembled from multiple RNAi studies used as classification template
by the BAGEL algorithm to call gene depletion significance [1].
}
\usage{data(BAGEL_nonEssential)}
\format{
 A vector of strings containing HGNC symbols of reference non essential genes.
}
\references{
[1] BAGEL: a computational framework for identifying essential genes from pooled library screens. Traver Hart and Jason Moffat. BMC Bioinformatics, 2016 vol. 17 p. 164.
}
\examples{
data(BAGEL_nonEssential)
head(BAGEL_nonEssential)
}
\seealso{
    \code{\link{BAGEL_essential}}
}
\keyword{datasets}
