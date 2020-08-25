\name{AVANA_Library}
\alias{AVANA_Library}
\docType{data}
\title{
sgRNAs' genome-wide annotation for the AVANA sgRNA library
}
\description{
A data frame with a named row for each sgRNA of the AVANA sgRNA library [1] including annotations such as targeted genes, and genomic coordinates.
}
\usage{data("AVANA_Library")}
\format{
   A a row named data frame with 71482 observations (one for each sgRNA) of the following 7 variables.
  \describe{
    \item{\code{CODE}}{alphanumerical identifier of the sgRNAs;}
    \item{\code{GENES}}{targeted gene;}
    \item{\code{EXONE}}{exone of the targeted genomic region (string with 'ex' prefix followed by the exone number);}
    \item{\code{CHRM}}{chromosome of where the targeted region resides (string)}
    \item{\code{STRAND}}{targeted DNA strand ('+' or '-')}
    \item{\code{STARTpos}}{starting genomic coordinate of the targeted genomic region (numeric);}
    \item{\code{ENDpos}}{ending genomic coordinate of the targeted genomic region (numeric).}
  }
}
\references{
[1] Meyers RM, Bryan JG, McFarland JM, Weir BA. Computational correction of copy number effect improves specificity of CRISPR–Cas9 essentiality screens in cancer cells. Nature. 2017.
}
\examples{
data(AVANA_Library)
head(AVANA_Library)
}
\keyword{datasets}
