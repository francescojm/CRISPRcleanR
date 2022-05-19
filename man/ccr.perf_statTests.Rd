\name{ccr.perf_statTests}
\alias{ccr.perf_statTests}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CRISPRcleanR correction assessment: Statistical tests
}
\description{
    This function tests the log fold changes of sgRNAs targeting different sets of genes for statistically significant differences with respect to background pre and post CRISPRcleanR correction, creating two sets of boxplots with outcomes and outputting statistical indicators.
}
\usage{
ccr.perf_statTests(cellLine, libraryAnnotation, correctedFCs,
                   outDir = "./",
                   GDSC.geneLevCNA = NULL,
                   CCLE.gisticCNA = NULL,
                   RNAseq.fpkms = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cellLine}{
    A string specifying the name of a cell line (or a COSMIC identifier [1]);
}
  \item{libraryAnnotation}{
    The sgRNA library annotations formatted as specified in the reference manual entry of the \code{\link{KY_Library_v1.0}} built in library.
}
  \item{correctedFCs}{
    sgRNAs log fold changes corrected for gene independent responses to CRISPR-Cas9 targeting, generated with the function \code{ccr.GWclean} (first data frame included in the list outputted by \code{ccr.GWclean}, i.e. \code{corrected_logFCs}).
}
  \item{outDir}{
    The path of the folder where the boxplot will be saved.
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
}
\details{
    This functions assess the statistical difference pre/post CRISPRcleanR correction of log fold changes for sgRNAs targeting respectively: 
    \itemize{
        \item copy number (CN) deleted genes according to the GDSC1000 repository
        \item CN deleted genes (gistic score = -2) according to the CCLE repository
        \item non expressed genes (FPKM < 0.05)
        \item genes with gistic score = 1
        \item genes with gistic score = 2
        \item non espressed genes (FPKM < 0.05) with gistic score = 1
        \item non espressed genes (FPKM < 0.05) with gistic score = 2
        \item genes with minimal CN = 2, according to the GDSC1000
        \item genes with minimal CN = 4, according to the GDSC1000
        \item genes with minimal CN = 8, according to the GDSC1000
        \item genes with minimal CN = 10, according to the GDSC1000
        \item non expressed genes (FPKM < 0.05) with minimal CN = 2, according to the GDSC1000
        \item non expressed genes (FPKM < 0.05) with minimal CN = 4, according to the GDSC1000
        \item non expressed genes (FPKM < 0.05) with minimal CN = 8, according to the GDSC1000
        \item non expressed genes (FPKM < 0.05) with minimal CN = 10, according to the GDSC1000
        \item core fitness essential genes, assembling signatures from MsigDB [5], included in the builtin vectors \code{EssGenes.DNA_REPLICATION_cons}, \code{EssGenes.KEGG_rna_polymerase}, \cr \code{EssGenes.PROTEASOME_cons}, \code{EssGenes.ribosomalProteins}, \cr \code{EssGenes.SPLICEOSOME_cons}
        \item Reference core fitness essential genes assembled from multiple RNAi studies used as classification template by the BAGEL algorithm to call gene depletion significance [6] \cr (\code{\link{BAGEL_essential}})
        \item Reference core fitness essential genes assembled from multiple RNAi studies used as classification template by the BAGEL algorithm to call gene depletion significance [6] after the removal core fitness essential genes from MsigDB [5]
        \item Reference non essential genes assembled from multiple RNAi studies used as classification template by the BAGEL algorithm to call gene depletion significance [6] \cr (\code{\link{BAGEL_nonEssential}})
}}
\value{
A list of three named 2x19 matrices, with one entry per statistical test, rows indicating pre/post CRISPRcleanR correction sgRNAs' log fold changes and one column per each tested gene set. In each matrix the entries contains, respectively
  \item{PVALS}{Pvalue resulting from a Student's t-test assessing the differences between sgRNAs log fold changes pre (first row) and post (second row) CRISPRcleanR correction with respect to background}
  \item{SIGNS}{The sign of the difference (1 = mean log fold change of the tested set larger that the mean of the background population, -1 = mean log fold change of the tested set smaller than the mean of the background population)}
  \item{EFFsizes}{Effect size (computing via the Cohen's D): difference of the means / pooled standard deviation.}
}
\source{
[a] \url{ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/Gene_level_CN.xlsx}.\cr
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
\examples{

## loading corrected sgRNAs log fold-changes and segment annotations for an example
## cell line (EPLC-272H)
data(EPLC.272HcorrectedFCs)

## loading library annotation
data(KY_Library_v1.0)

## Evaluate correction effects. Boxplots will be saved in EPLC-272H.pdf
## in the current directory
RES<-ccr.perf_statTests('EPLC-272H',libraryAnnotation = KY_Library_v1.0,
                   correctedFCs = EPLC.272HcorrectedFCs$corrected_logFCs)
RES$PVALS
RES$EFFsizes
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Assessment and Visualisation}% use one of  RShowDoc("KEYWORDS")
