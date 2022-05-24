\name{ccr.genes2sgRNAs}
\alias{ccr.genes2sgRNAs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Targeting sgRNAs
}
\description{
    This function returns the set of sgRNAs targeting the set of genes provided in input, in a given pooled library.
}
\usage{
ccr.genes2sgRNAs(libraryAnnotation,genes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{libraryAnnotation}{A data frame with a named row for each sgRNA with the same format of \code{\link{KY_Library_v1.0}}}
  \item{genes}{A list of strings containing HGNC symbols}
}
\value{
A list of strings containing the identifiers of the sgRNAs targeting the inputted set of genes
}
\author{
Francesco Iorio (francesco.iorio@fht.org)
}
\examples{\dontshow{setwd(tempdir())}
    ## Loading an sgRNA pooled library annotation
    data(KY_Library_v1.0)
    ## Loading an example set of genes
    data(BAGEL_essential)
    
    ccr.genes2sgRNAs(KY_Library_v1.0,BAGEL_essential)
}
\seealso{
    \code{\link{KY_Library_v1.0}}
}
\keyword{utils}% use one of  RShowDoc("KEYWORDS")
