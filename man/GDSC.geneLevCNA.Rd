\name{GDSC.geneLevCNA}
\alias{GDSC.geneLevCNA}
\docType{data}
\title{
Genome-wide copy number data for 15 human cancer cell lines.
}
\description{
Genome-wide copy number data derived from PICNIC analysis of Affymetrix SNP6 segmentation data (EGAS00001000978, part of the Genomics of Drug Sensitivity in 1,000 Cancer Cell Lines (GDSC1000) panel [1]) for 15 cell lines used in [2] to assess CRISPRcleaneR results. 
}
\usage{data(GDSC.geneLevCNA)}
\format{
  A data frame with HGNC gene symbols on the row cancer cell lines' cosmic identifiers on the columns. The entry in position
  \emph{i,j} indicates the copy number status of gene \emph{i} in cell line \emph{j}. 
}
\details{
Each entry of the data frame is a string made of four comma seperated peices of data (\code{n1,n2,n3,n4}), hyphen (\code{-}) is used when the corresponding data is unknown.

The four values indicate:
\itemize{
            \item \code{n1}: Maximum copy number of any genomic segment containing coding sequence of the gene (\code{-1} indicates a value could not be assigned).
            \item \code{n2}: Minimum copy number of any genomic segment containing coding sequence of the gene (\code{-1} indicates a value could not be assigned).
            \item \code{n3}: Zygosity - (\code{H}) if all segments containing gene sequence are heterozygous, (\code{L}) if any segment containing coding sequence has LOH, (\code{0}) if the complete coding sequence of the gene falls within a homozygous deletion.
            \item \code{n4}: Disruption (\code{D}) if the gene spans more than 1 genomic segment (\code{-}) if no disruption occures.
            }

}
\source{
This data frame has been derived from the xls table available at \url{ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/Gene_level_CN.xlsx}.\cr \cr
}
\references{
[1] Iorio F, Knijnenburg TA, Vis DJ, Bignell GR, Menden MP, et al. \emph{A landscape of pharmacogenomic interactions in cancer}
Cell 2016 Jul 28;166(3):740-54

[2] Iorio, F., Behan, F. M., Goncalves, E., Beaver, C., Ansari, R., Pooley, R., et al. (n.d.). Unsupervised correction of gene-independent cell responses to CRISPR-Cas9 targeting. \cr
http://doi.org/10.1101/228189
}
\examples{
data(GDSC.geneLevCNA)
GDSC.geneLevCNA[1:10,1:10]
}
\keyword{datasets}
