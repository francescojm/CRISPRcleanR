\name{ccr.get.nonExpGenes}
\alias{ccr.get.nonExpGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Non expressed genes in a given cell line}
\description{
This function takes in input the name (or the COSMIC identifier [1]) of a cell line and it identifies genes that are not expressed (according to a user defined FPKM threshold) using a collection of RNAseq profile from [2].
}
\usage{
ccr.get.nonExpGenes(cellLine, th = 0.05,
                    amplified = FALSE, minCN = 8,
                    RNAseq.fpkms=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cellLine}{
 A string specifying the name of a cell line (or a COSMIC identifier [1]);
}
  \item{th}{
Minimum FPKM value for a gene to be considered as expressed;
}
  \item{amplified}{
A logic value specifying whether the selected not expressed genes should be also copy number amplified function;
}
  \item{minCN}{
If \code{amplified = TRUE}, this parameter defines a lower threshold for the minimum copy number of any genomic segment containing coding sequence of a gene in order for it to be considered as copy number amplified.}
  \item{RNAseq.fpkms}{
  Genome-wide substitute reads with fragments per kilobase of exon per million reads mapped (FPKM) across cell lines. These can be derived from a comprehensive collection of RNAseq profiles described in [2]. The format must be the same of the \code{\link{RNAseq.fpkms}} builtin data frame. If NULL then this function uses the \code{\link{RNAseq.fpkms}} builtin data fram containing data for 15 cell lines used in [3] to assess CRISPRcleaneR results.
}
}
\value{
A vector of string containing the HGNC symbols of non expressed (optionally copy number amplified) genes in the cell line under consideration.
}
\references{
[1] Forbes SA, Beare D, Boutselakis H, et al. \emph{COSMIC: somatic cancer genetics at high-resolution}
Nucleic Acids Research, Volume 45, Issue D1, 4 January 2017, Pages D777-D783.

[2] Garcia-Alonso L, Iorio F, Matchan A, et al. \emph{Transcription factor activities enhance markers of drug response in cancer}
doi: https://doi.org/10.1101/129478

[3] Iorio, F., Behan, F. M., Goncalves, E., Beaver, C., Ansari, R., Pooley, R., et al.\cr
Unsupervised correction of gene-independent cell responses to CRISPR-Cas9 targeting.\cr
http://doi.org/10.1101/228189
}
\author{
Francesco Iorio (francesco.iorio@fht.org)
}
\seealso{
\code{\link{ccr.get.gdsc1000.AMPgenes}}
}
\examples{
    ccr.get.nonExpGenes('HT-29',amplified = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utils}