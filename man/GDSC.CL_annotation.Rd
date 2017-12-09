\name{GDSC.CL_annotation}
\alias{GDSC.CL_annotation}
\docType{data}
\title{
Tissue type and other annotations for 1,001 human cancer cell lines
}
\description{
Tissue type and other annotations for 1,001 human cancer cell lines
}
\usage{data(GDSC.CL_annotation)}
\format{
  A data frame with 1,001 observations of the following 7 variables.
  \describe{
    \item{\code{CL.name}}{Cell line name;}
    \item{\code{COSMIC.ID}}{Cosmic identifier of the cell line;}
    \item{\code{GDSC.description_1}}{Tissue descriptor (Genomics of Drug Sensitivity in Cancer - Level 1);}
    \item{\code{GDSC_description_2}}{Tissue descriptor (Genomics of Drug Sensitivity in Cancer - Level 2);}
    \item{\samp{TCGA type}}{Manaually curated matched TCGA cancer type;}
    \item{\code{MMR}}{Microsatellite instability status (MSI-S = Stable, MSI-L = Instable, MSI-H = highly-Instable).}
  }
}
\source{
This data frame has been derived from the xls table available at \url{http://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/suppData/TableS1E.xlsx}.
}
\references{
[1] Iorio F, Knijnenburg TA, Vis DJ, Bignell GR, Menden MP, et al. A landscape of pharmacogenomic interactions in cancer Cell 2016 Jul 28;166(3):740-54
}
\examples{
data(GDSC.CL_annotation)
head(GDSC.CL_annotation)
}
\keyword{datasets}
