\name{CCLE.gisticCNA}
\alias{CCLE.gisticCNA}
\docType{data}
\title{
Genome-wide copy number data for 13 human cancer cell lines.
}
\description{
Genome-wide Gistic [1] scores quantifying copy number status across a subset of the cell lines in
\code{\link{CL.subset}} that are used to assess CRISPRcleaneR results in [2].
}
\usage{data(CCLE.gisticCNA)}
\format{
  A data frame with one observations per gene across 13 variables (one per cell line).
  Row names indicate HGNC gene symbols and column names indicate cell line COSMIC identifiers [3].
}
\source{
This data frame has been derived from the tsv file downloadable at \cr
\url{http://www.cbioportal.org/study?id=cellline_ccle_broad#summary}.\cr
This has been obtained by processing Affymetrix SNP array data in the Cancer 
Cell Line Encyclopaedia [4] repository \cr 
(\url{https://portals.broadinstitute.org/ccle_legacy/data/}) 
}
\references{
[1] Mermel CH, Schumacher SE, Hill B, et al. \emph{GISTIC2.0 facilitates sensitive and confident localization of the targets of focal somatic copy-number alteration in human cancers.} Genome Biol. 2011;12(4):R41. doi: 10.1186/gb-2011-12-4-r41. \cr

[2] Iorio, F., Behan, F. M., Goncalves, E., Beaver, C., Ansari, R., Pooley, R., et al. (n.d.). Unsupervised correction of gene-independent cell responses to CRISPR-Cas9 targeting. \cr
http://doi.org/10.1101/228189

[2] Forbes SA, Beare D, Boutselakis H, et al. \emph{COSMIC: somatic cancer genetics at high-resolution}
Nucleic Acids Research, Volume 45, Issue D1, 4 January 2017, Pages D777-D783, \cr

[3] Barretina J, Caponigro G, Stransky N, et al. \emph{The Cancer Cell Line Encyclopedia enables predictive modelling of anticancer drug sensitivity.} Nature. 2012 Mar 28;483(7391):603-7. doi: 10.1038/nature11003. Erratum in: Nature. 2012 Dec 13;492(7428):290.
}
\examples{
data(CCLE.gisticCNA)
head(CCLE.gisticCNA)
}
\keyword{datasets}
