\name{ccr.impactOnPhenotype}
\alias{ccr.impactOnPhenotype}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Assessing the impact and potential distortion introduced by the CRISPRcleanR correction on the genes showing loss/gain-of-fitness effect.
}
\description{
    This function compares two MAGeCK [1] gene summaries (obtained from sgRNA count files pre/post CRISPRcleanR correction) and it computes the percentages of genes whose loss/gain-of-fitness effect is attenuated post CRISPRcleanR correction or potentially distorted (i.e. loss-of-fitness genes are detected post CRISPRcleanR correction as gain-of-fitness genes, and viceversa).
    Results are returned in output and optionally plotted as bar/pie charts.
}
\usage{
ccr.impactOnPhenotype(MO_uncorrectedFile,
                      MO_correctedFile,
                      sigFDR = 0.05,
                      expName = "expName",
                      display = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MO_uncorrectedFile}{
    String specifying the path to a MAGeCK gene summary file produced by MAGeCK from non corrected sgRNA counts.
}
  \item{MO_correctedFile}{
    String specifying the path to a MAGeCK gene summary file produced by MAGeCK from CRISPRcleanR corrected sgRNA counts.
}
  \item{sigFDR}{
    A numerical value in [0,1] False discovery rate threshold at which genes are called as significantly exerting a loss/gain-of-fitness effect.
}
  \item{expName}{
    A string specifying the experiment name, used as main title in the figures (ignored if the \code{display} argument is set to \code{FALSE}).    
}
  \item{display}{
    Boolean value specifying whether figures sumarising the comparison results should be plotted.
}
}
\value{
 A list containing the following four numerical values and two data frames:
        \itemize{
            \item \code{GW_impact \%}: Percentage of genes impacted by the CRISPRcleanR correction, i.e. showing a gain/loss-of-fitness genes effect in the MAGeCK gene summary obtained from uncorrected sgRNA counts, over the total number of screened genes;
            \item \code{Phenotype_G_impact \%}: Percentage of genes impacted by the CRISPRcleanR correction, i.e. showing a gain/loss-of-fitness genes effect in the MAGeCK gene summary obtained from uncorrected sgRNA counts, over the total number of genes showing a gain/loss of fitness effect in the MAGeCK gene summary obtained from uncorrected sgRNA counts;
            \item \code{GW_distortion \%}: Percentage of genes distorted by the CRISPRcleanR correction, i.e. showing a gain/loss-of-fitness effect in the MAGeCK gene summary obtained from corrected sgRNA counts that is opposite to the effect in that obtained from uncorrected sgRNA counts, over the total number of screened genes;
            \item \code{Phenotype_G_distortion \%}: Percentage of genes distorted by the CRISPRcleanR correction, i.e. showing a gain/loss-of-fitness effect in the MAGeCK gene summary obtained from corrected sgRNA counts that is opposite to the effect in that obtained from uncorrected sgRNA counts, over the total number of screened genes,over the total number of genes showing a gain/loss of fitness effect in the MAGeCK gene summary obtained from uncorrected sgRNA countsl;
            \item \code{geneCounts}: A contingency table with gene counts as entries, with data referring to the original (uncorrected) sgRNA counts on the columns, and to the corrected sgRNA counts on the rows. There are three vectors for each dimensions, respectively for number of genes showing a significant loss of fitness effect (dep.), number of genes not showing any fitness effect (or with a not clear effect, i.e. showing both gain and loss of fitness effect, null), and number of genes showing a significant gain of fitness effect (enr.);
            \item \code{distortion}: a data frame showing genes whose fitness effect has been distorted by the CRISPRcleanR correction: one row per gene (as specified by the row names), with two column per condition (i.e. prior/post correction), indicating the loss of fitness effect fdr (neg.fdr and ccr.neg.fdr) and the gain of fitness effect fdr (pos.fdr and ccr.pos.fdr) as outputted by MAGeCK;
            \item \code{attenuation}: a data frame showing genes whose fitness effect has been attenuated by the CRISPRcleanR correction: one row per gene (as specified by the row names), with two column per condition (i.e. prior/post correction), indicating the loss of fitness effect fdr (neg.fdr and ccr.neg.fdr) and the gain of fitness effect fdr (pos.fdr and ccr.pos.fdr) as outputted by MAGeCK;
            }
}
\details{
For each of the considered MAGeCK gene summaries, this function calls loss/gain-of-fitness based on the MAGeCK negative/positive false discovery rate and the user defined threshold (as specified by the \code{sigFDR} argument). Particularly, are called as significant loss-of-fitness genes those with a negative fdr < \code{sigFDR} and a positive fdr >= \code{sigFDR}, and as significant gain-of-fitness genes those those with a positive fdr < \code{sigFDR} and a negative fdr >= \code{sigFDR}. All the other genes are deemed as not exerting any effect on cellular fitness.
}
\references{
[1] Li, W., Xu, H., Xiao, T., Cong, L., Love, M. I., Zhang, F., et al. (2014). MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biology, 15(12), 554. [2] Hart, T., & Moffat, J. (2016). BAGEL: a computational framework for identifying essential genes from pooled library screens. BMC Bioinformatics, 17(1), 164.
}
\author{
Francesco Iorio (francesco.iorio@fht.org)
}
\seealso{
    \code{\link{ccr.ExecuteMageck}}
}
\examples{
## Loading sgRNA library annotation file
data(KY_Library_v1.0)

## Deriving the path of the file with the example dataset,
## from the mutagenesis of the EPLC-272H colorectal cancer cell line
fn<-paste(system.file('extdata', package = 'CRISPRcleanR'),
          '/EPLC-272H_counts.tsv',sep='')

## Loading, median-normalizing and computing fold-changes for the example dataset
normANDfcs<-ccr.NormfoldChanges(fn,min_reads=30,
                                EXPname='EPLC-272H',
                                libraryAnnotation = KY_Library_v1.0)

uncorrected_fn<-ccr.PlainTsvFile(sgRNA_count_object = normANDfcs$norm_counts,
                                 fprefix = 'EPLC-272H')

## execute MAGeCK on uncorrected normalised counts
uncorrected_gs_fn<-ccr.ExecuteMageck(mgckInputFile = uncorrected_fn,
                                     expName = 'EPLC-272H',
                                     normMethod = 'none')

## Genome-sorting of the fold changes
gwSortedFCs<-ccr.logFCs2chromPos(normANDfcs$logFCs,KY_Library_v1.0)

## Identifying and correcting biased sgRNAs' fold changes
correctedFCs<-ccr.GWclean(gwSortedFCs,display=FALSE,label='EPLC-272H')

## correcting individual sgRNA treatment counts
correctedCounts<-ccr.correctCounts('EPLC-272H',normANDfcs$norm_counts,
                                   correctedFCs,
                                   KY_Library_v1.0,
                                   minTargetedGenes=3,
                                   OutDir='./')

## saving corrected/uncorrected sgRNA count files as plain tsv files
corrected_fn<-ccr.PlainTsvFile(sgRNA_count_object = correctedCounts,
                               fprefix = 'EPLC-272H_ccleaned')

## execute MAGeCK on corrected normalised counts
## - it requires MAGeCK to be pre-installed - 
corrected_gs_fn<-ccr.ExecuteMageck(mgckInputFile = corrected_fn,
                                   expName = 'EPLC-272H_ccleaned')

## If MAGeCK is installed and correctly executed then
## Assessing the impact of CRISPcleanR correction on gain/loss-of-fitness genes

if (corrected_gs_fn!=''){

  RES<-ccr.impactOnPhenotype(MO_uncorrectedFile = uncorrected_gs_fn,
                        MO_correctedFile = corrected_gs_fn,
                        expName = 'EPLC-272H')

  ## Percentage of genes whose gain/loss-of fitness effect is impacted by CRISPRcleanR 
  ## over the total number of screened genes
  RES[1]  

  ## Percentage of genes whose gain/loss-of fitness effect is impacted by CRISPRcleanR 
  ## over the total number of genes with a significant gain/loss-of fitness effect when 
  ## using uncorrected sgRNA counts
  RES[2]

  ## Percentage of genes whose gain/loss-of fitness effect is distorted by CRISPRcleanR 
  ## over the total number of screened genes
  RES[3]

  ## Percentage of genes whose gain/loss-of fitness effect is distorted by CRISPRcleanR 
  ## over the total number of genes with a significant gain/loss-of fitness effect when 
  ## using uncorrected sgRNA counts
  RES[4]

  ## Contingency table showing the impact of the CRISPRcleanR correction on the phenotype
  RES$geneCounts

  ## Genes whose gain/loss-of-fitness effect has been distorted by the CRISPRcleanR correction
  RES$distortion
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Assessment and Visualisation}% use one of  RShowDoc("KEYWORDS")
