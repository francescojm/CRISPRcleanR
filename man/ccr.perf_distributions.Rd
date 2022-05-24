\name{ccr.perf_distributions}
\alias{ccr.perf_distributions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CRISPRcleanR correction assessment: inspection of sgRNA log fold changes distributions
}
\description{
  This function creates distributions density plots of sgRNA log fold changes for defined sets of targeted genes
  prior/post CRISPRcleanR correction.
}
\usage{
ccr.perf_distributions(cellLine, correctedFCs,
                       GDSC.geneLevCNA = NULL,
                       CCLE.gisticCNA = NULL,
                       RNAseq.fpkms = NULL,
                       minCNs = c(8, 10),
                       libraryAnnotation,
                       GDSC.CL_annotation=NULL)
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
  \item{CCLE.gisticCNA}{
Genome-wide Gistic [3] scores quantifying copy number status across cell lines with the same format of \code{\link{CCLE.gisticCNA}}. If NULL then this function uses the \code{\link{CCLE.gisticCNA}} builtin data frame, containing data for 13 cell lines of the 15 used in [2] to assess the performances of CRISPRcleanR.

}
  \item{RNAseq.fpkms}{
Genome-wide substitute reads with fragments per kilobase of exon per million reads mapped (FPKM) across cell lines. These can be derived from a comprehensive collection of RNAseq profiles described in [4]. The format must be the same of the \code{\link{RNAseq.fpkms}} builtin data frame. If NULL then this function uses the \code{\link{RNAseq.fpkms}} builtin data fram containing data for 15 cell lines used in [2] to assess CRISPRcleaneR results.
}
  \item{minCNs}{
  A numerical vector with two entries specifying the minimal copy number for a gene in order to be considered amplified
  based on the data in \code{GDSC.geneLevCNA}. These two values can be 2, 4, 8 or 10.
}
\item{libraryAnnotation}{
    The sgRNA library annotations formatted as specified in the reference manual entry of the \code{\link{KY_Library_v1.0}} built in library.
}
  \item{GDSC.CL_annotation}{
    Cell lines annotation dataframe with the same structure of the \code{\link{GDSC.CL_annotation}}. If NULL then the \code{\link{GDSC.CL_annotation}} is used.
  }
}
\details{
This function generates 4 sets of plots. They contains log fold change distributions density plots prior/post
CRISPRcleanR correction respectively for
\itemize{
    \item (i) Copy number amplified genes according to the data in \code{GDSC.geneLevCNA} based on the two threshold values specified in \code{minCNs};
    \item (ii) Copy number amplified genes according to the data in \code{CCLE.gisticCNA} (gistic score = +2);
    \item (iii) Copy number amplified non expressed genes according to the data in \code{GDSC.geneLevCNA} based on the two threshold values specified in \code{minCNs}, and the data in \code{RNAseq.fpkms} (FPKM < 0.05);
    \item (iv) reference sets of core fitness essential genes from MSigDB [5] (included in the builtin vectors \code{EssGenes.DNA_REPLICATION_cons}, \code{EssGenes.KEGG_rna_polymerase}, \cr \code{EssGenes.PROTEASOME_cons}, \code{EssGenes.ribosomalProteins}, \cr \code{EssGenes.SPLICEOSOME_cons}, and reference core-fitness-essential and non-essential genes assembled from multiple RNAi studies used as classification template by the BAGEL algorithm to call gene depletion significance [6] \cr (\code{\link{BAGEL_essential}}, \code{\link{BAGEL_nonEssential}}).
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

[5] Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., et al. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proceedings of the National Academy of Sciences of the United States of America, 102(43), 15545-15550. http://doi.org/10.1073/pnas.0506580102

[6] BAGEL: a computational framework for identifying essential genes from pooled library screens. Traver Hart and Jason Moffat. BMC Bioinformatics, 2016 vol. 17 p. 164.
}
\author{
Francesco Iorio (francesco.iorio@fht.org)
}
\seealso{
\code{\link{KY_Library_v1.0}}, \code{\link{ccr.GWclean}}, \cr
\code{\link{GDSC.geneLevCNA}}, \code{\link{CCLE.gisticCNA}}, \code{\link{RNAseq.fpkms}}, \cr
\code{\link{EssGenes.DNA_REPLICATION_cons}}, \code{\link{EssGenes.KEGG_rna_polymerase}}, \code{\link{EssGenes.PROTEASOME_cons}}, \code{\link{EssGenes.ribosomalProteins}}, \code{\link{EssGenes.SPLICEOSOME_cons}} \cr
\code{\link{BAGEL_essential}}, \code{\link{BAGEL_nonEssential}}
}
\examples{\dontshow{setwd(tempdir())}
## loading corrected sgRNAs log fold-changes and segment annotations for an example
## cell line (HT-29)
data(HT.29correctedFCs)

## loading library annotation
data(KY_Library_v1.0)

## inpecting sgRNA log fold change distributions prior/post CRISPRcleanR correction
ccr.perf_distributions('HT-29',HT.29correctedFCs$corrected_logFCs,
                       libraryAnnotation = KY_Library_v1.0)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Assessment and Visualisation}% use one of  RShowDoc("KEYWORDS")
