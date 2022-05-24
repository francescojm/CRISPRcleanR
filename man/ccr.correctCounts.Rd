\name{ccr.correctCounts}
\alias{ccr.correctCounts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Correction of sgRNA treatment counts for gene independent responses to CRISPR-Cas9 targeting
}
\description{
This function applies an inverse transformation (described in [1]) to CRISPRcleanR corrected sgRNAs' log fold changes and produces in output normalised corrected sgRNA counts (across treatments and control replicates), suitable for gene depletion/enrichment statistical testing via mean-variance modeling (for example through MAGeCK [2]*).
*MAGeCK should be executed excluding initial normalisation, as the corrected sgRNA counts outputted by this function are already normalised.
}
\usage{
ccr.correctCounts(CL,normalised_counts,
                  correctedFCs_and_segments,
                  libraryAnnotation,
                  minTargetedGenes=3,
                  OutDir='./',
                  ncontrols=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CL}{A string specifying the name of the experiment. This will be used to compose names of files and folde where results will be saved.
}
  \item{normalised_counts}{
    A data frame containing normalised sgRNAs' read counts, which can be computed using the \code{ccr.NormfoldChanges} function from raw sgRNAs' counts.
}
  \item{correctedFCs_and_segments}{
    sgRNAs log fold changes corrected for gene independent responses, generated with the function \code{ccr.GWclean}.
}
\item{libraryAnnotation}{
    A data frame containing the sgRNAs' genome-wide annotations with at least a named row for each of the sgRNAs included in the foldchanges data frame provided in input.
    The following columns/headers should be present in this data frame (additional columns will be ignored):
    \itemize{
            \item GENES: string vector containing the HGNC symbols of the genes targeted by the sgRNA under consideration;
            \item EXONE: string vector containing the gene exon targeted by the sgRNA under consideration (these should include the prefix "ex" followed by the exone number);
            \item CHRM: string vector the chromosome of the gene targeted by the sgRNA under consideration (X and Y chromosome should be specified as "X" and "Y");
            \item STRAND: string vector containing the strand targeted by the sgRNA under consideration ("+" or "-");
            \item STARTpop: numeric vector containing the genomic coordinate of the starting position of the region targeted by the sgRNA under consideration;
            \item ENDpos: numeric vector containing the genomic coordinate of the ending position of the region targeted by the sgRNA under consideration;
            }}
  \item{minTargetedGenes}{
    Minimanl number of different genes targeted by sgRNAs in a biased segment in order for the corresponding counts to be corrected (default = 3).
}
  \item{OutDir}{
    Path of the folder where results and plots will be saved.
}
  \item{ncontrols}{
  A numerical value indicating the number of control replicates (therefore columns to be considered as controls in the normalised counts).
}
}
\value{
    A data frame with one entry per sgRNA and individual columns for the control/treatment samples included in the normalised count data object specified by the \code{normalised_counts} parameter, and containing sgRNA counts corrected for gene independent responses to CRISPR-Cas9 targeting and median-ratio normalised.
}
\references{
[1] Iorio F, Behan FM, Goncalves E, Bhosle SG, Chen E, Shepherd R, Beaver C, Ansari R, Pooley R, Wilkinson P, Harper S, Butler AP, Stronach EA, Saez-Rodriguez J, Yusa K, Garnett MJ. Unsupervised correction of gene-independent cell responses to CRISPR-Cas9 targeting. BMC Genomics. 2018 Aug 13;19(1):604. doi: 10.1186/s12864-018-4989-y.

[2] Li, W., Xu, H., Xiao, T., Cong, L., Love, M. I., Zhang, F., et al. (2014). MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biology, 15(12), 554.
}
\author{
Francesco Iorio (francesco.iorio@fht.orgfht.org)
}
\seealso{
\code{\link{ccr.NormfoldChanges}}, \code{\link{ccr.GWclean}}
}
\examples{\dontrun{
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

head(correctedCounts)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{analysis}% use one of  RShowDoc("KEYWORDS")
