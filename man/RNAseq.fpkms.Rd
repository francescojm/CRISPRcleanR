\name{RNAseq.fpkms}
\alias{RNAseq.fpkms}
\docType{data}
\title{
    RNAseq derived genome-wide basal expression profiles for 15 cell lines.
}
\description{
    Genome-wide substitute reads with fragments per kilobase of exon per million reads mapped (FPKM) for the 15 cell lines specified in \code{\link{CL.subset}}, derived from a comprehensive
    collection of RNAseq profiles described in [1] and used in [2] to assess CRISPRcleaneR results.
}
\usage{data(RNAseq.fpkms)}
\format{
  A data frame with one bservations per gene and one variable per cell line. Row names indicates HGNC symbols and column names indicate cell line COSMIC identifiers [3].
}
\references{
[1] Garcia-Alonso L, Iorio F, Matchan A, et al. \emph{Transcription factor activities enhance markers of drug response in cancer}
doi: https://doi.org/10.1101/129478

[2] Iorio, F., Behan, F. M., Goncalves, E., Beaver, C., Ansari, R., Pooley, R., et al. (n.d.). Unsupervised correction of gene-independent cell responses to CRISPR-Cas9 targeting. \cr
http://doi.org/10.1101/228189

[3] Forbes SA, Beare D, Boutselakis H, et al. \emph{COSMIC: somatic cancer genetics at high-resolution}
Nucleic Acids Research, Volume 45, Issue D1, 4 January 2017, Pages D777-D783, \cr
}
\examples{
data(RNAseq.fpkms)
head(RNAseq.fpkms)
}
\seealso{\code{\link{CL.subset}}}
\keyword{datasets}
