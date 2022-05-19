\name{ccr.RecallCurves}
\alias{ccr.RecallCurves}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CRISPRcleanR correction assessment: Recall curve inspection
}
\description{
This function creates plots with Recall curve outcomes (a it computes areas under the Recall curves) resulting
from classifying defined sets of sgRNAs (respectively genes) based on their log fold change (respectively log fold changes averaged across targeting sgRNAs).
}
\usage{
ccr.RecallCurves(cellLine, correctedFCs, GDSC.geneLevCNA = NULL,
                 RNAseq.fpkms = NULL, minCN = 8, libraryAnnotation,
                 GeneLev = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cellLine}{
A string specifying the name of a cell line (or a COSMIC identifier [1]);
}
  \item{correctedFCs}{
sgRNAs log fold changes corrected for gene independent responses to CRISPR-Cas9 targeting, generated with the function \code{ccr.GWclean} (first data frame included in the list outputted by \code{ccr.GWclean}, i.e. \code{corrected_logFCs}).
}
  \item{GDSC.geneLevCNA}{
Genome-wide copy number data with the same format of \code{\link{GDSC.geneLevCNA}}. This can be assembled from the xls sheet specified in the source section [a] (containing data for the GDSC1000 cell lines). If NULL, then this function uses the built in \code{\link{GDSC.geneLevCNA}} data frame, containing data derived from [a] for 15 cell lines used in [2] to assess the performances of CRISPRcleanR.
}
  \item{RNAseq.fpkms}{
Genome-wide substitute reads with fragments per kilobase of exon per million reads mapped (FPKM) across cell lines. These can be derived from a comprehensive collection of RNAseq profiles described in [4]. The format must be the same of the \code{\link{RNAseq.fpkms}} builtin data frame. If NULL then this function uses the \code{\link{RNAseq.fpkms}} builtin data fram containing data for 15 cell lines used in [2] to assess CRISPRcleaneR results.

}
  \item{minCN}{
 A numerical value specifying the minimal copy number for a gene in order to be considered amplified
  based on the data in \code{GDSC.geneLevCNA}. This value can be 2, 4, 8 or 10.
}
  \item{libraryAnnotation}{
   The sgRNA library annotations formatted as specified in the reference manual entry of the \code{\link{KY_Library_v1.0}} built in library.
}
  \item{GeneLev}{
A logical value specifying if the Recall should be computed at level of genes. In this case average gene log fold changes are computed from the inputted corrected log fold changes across targeting sgRNAs. 
}
}
\details{
This function generates 2 plots, showing Recall curves resulting from classifying the following 4 sets of sgRNAs (or Genes, depending on the parameter \code{GeneLev}, based on their log fold changes (or log fold changes averaged across targeting guides):
\itemize{
    \item (i) Copy number amplified genes according to the data in \code{GDSC.geneLevCNA} based on the threshold value specified in \code{minCNs};
    \item (ii) Copy number amplified non expressed genes according to the data in \code{GDSC.geneLevCNA} based on the threshold value specified in \code{minCNs}, and the data in \code{RNAseq.fpkms} (FPKM < 0.05);
    \item (iv) reference sets of core-fitness-essential and non-essential genes assembled from multiple RNAi studies used as classification template by the BAGEL algorithm to call gene depletion significance [5] \cr (\code{\link{BAGEL_essential}}, \code{\link{BAGEL_nonEssential}}).
    }
}
\source{
[a] \url{ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/Gene_level_CN.xlsx}.\cr \cr
}
\references{
[1] Forbes SA, Beare D, Boutselakis H, et al. \emph{COSMIC: somatic cancer genetics at high-resolution}
Nucleic Acids Research, Volume 45, Issue D1, 4 January 2017, Pages D777-D783.

[2] Iorio, F., Behan, F. M., Goncalves, E., Beaver, C., Ansari, R., Pooley, R., et al. (n.d.). Unsupervised correction of gene-independent cell responses to CRISPR-Cas9 targeting. \cr
http://doi.org/10.1101/228189

[3] Mermel CH, Schumacher SE, Hill B, et al. \emph{GISTIC2.0 facilitates sensitive and confident localization of the targets of focal somatic copy-number alteration in human cancers.} Genome Biol. 2011;12(4):R41. doi: 10.1186/gb-2011-12-4-r41. \cr

[4] Garcia-Alonso L, Iorio F, Matchan A, et al. \emph{Transcription factor activities enhance markers of drug response in cancer}
doi: https://doi.org/10.1101/129478

[5] BAGEL: a computational framework for identifying essential genes from pooled library screens. Traver Hart and Jason Moffat. BMC Bioinformatics, 2016 vol. 17 p. 164.
}
\author{
Francesco Iorio (francesco.iorio@fht.org)
}
\seealso{
\code{\link{KY_Library_v1.0}}, \code{\link{ccr.GWclean}}, \cr
\code{\link{GDSC.geneLevCNA}}, \code{\link{RNAseq.fpkms}}, \cr
\code{\link{BAGEL_essential}}, \code{\link{BAGEL_nonEssential}}
}
\examples{
## loading corrected sgRNAs log fold-changes and segment annotations for an example
## cell line (EPLC-272H)
data(EPLC.272HcorrectedFCs)

## loading library annotation
data(KY_Library_v1.0)

## Creating recall curve plots and computing corresponding underlying area
## at the level of sgRNAs
ccr.RecallCurves('EPLC-272H',EPLC.272HcorrectedFCs$corrected_logFCs,
                  libraryAnnotation=KY_Library_v1.0)

## Creating recall curve plots and computing corresponding underlying area
## at the gene level
ccr.RecallCurves('EPLC-272H',EPLC.272HcorrectedFCs$corrected_logFCs,
                  libraryAnnotation=KY_Library_v1.0,GeneLev = TRUE)
}
\keyword{Assessment and Visualisation}% use one of  RShowDoc("KEYWORDS")
