\name{CL.subset}
\alias{CL.subset}
\docType{data}
\title{
COSMIC identifiers of 15 immortalised human cancer cell lines
}
\description{
COSMIC identifiers [1] of 15 cell lines included in the GDSC1000 panel [2] that are used in [3] to assess CRISPRcleaneR results.
}
\usage{data(CL.subset)}
\format{
  A vector of strings.
  }
\references{
[1] Forbes SA, Beare D, Boutselakis H, et al. \emph{COSMIC: somatic cancer genetics at high-resolution}
Nucleic Acids Research, Volume 45, Issue D1, 4 January 2017, Pages D777-D783, \cr

[2] Iorio F, Knijnenburg TA, Vis DJ, Bignell GR, Menden MP, et al. \emph{A landscape of pharmacogenomic interactions in cancer}
Cell 2016 Jul 28;166(3):740-54

[3] Iorio, F., Behan, F. M., Goncalves, E., Beaver, C., Ansari, R., Pooley, R., et al. (n.d.). Unsupervised correction of gene-independent cell responses to CRISPR-Cas9 targeting. \cr
http://doi.org/10.1101/228189
}
\examples{
data(CL.subset)

## Loading annotation for the GDSC1000 cell lines
data(GDSC.CL_annotation)

## Visualising annotation
GDSC.CL_annotation[CL.subset,]
}
\keyword{datasets}
