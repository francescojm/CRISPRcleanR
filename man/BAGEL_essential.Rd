\name{BAGEL_essential}
\alias{BAGEL_essential}
\docType{data}
\title{
Reference Core fitness essential genes
}
\description{
A list of reference core fitness essential genes assembled from multiple RNAi studies used as classification template
by the BAGEL algorithm to call gene depletion significance [1].
}
\usage{data(BAGEL_essential)}
\format{
 A vector of strings containing HGNC symbols of reference core fitness essential genes.
}
\references{
[1] BAGEL: a computational framework for identifying essential genes from pooled library screens. Traver Hart and Jason Moffat. BMC Bioinformatics, 2016 vol. 17 p. 164.
}
\examples{
data(BAGEL_essential)
head(BAGEL_essential)
}
\seealso{
    \code{\link{BAGEL_nonEssential}}
}
\keyword{datasets}
