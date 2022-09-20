\name{ccr.NormfoldChanges}
\alias{ccr.NormfoldChanges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normalisation of sgRNA counts and fold change computation
}
\description{
    This function normalises sgRNAs' counts stored in a tsv file whose path is provided in input, to adjust for the effect of library size and read count distributions, scaling by the total number of reads per sample or using the gene wise median of ratios method [1]. It computes log fold changes of transfected library replicates versus controls (tipically the sgRNA counts in the plasmid). The output of this function is returned as a list, and it is also saved into two tsv files.
 }
\usage{
ccr.NormfoldChanges(filename,
                    Dframe=NULL, display=TRUE,
                    saveToFig=FALSE,
                    outdir='./', min_reads=30, EXPname='',
                    libraryAnnotation, ncontrols=1,
                    method='ScalingByTotalReads')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
        A string specifying the path of a tsv file containing the raw sgRNA counts. This must be a tab delimited file with one row per sgRNA and the following columns/headers: 
        \itemize{
            \item sgRNA: containing alphanumerical identifiers of the sgRNA under consideration;
            \item gene: containing HGNC symbols of the genes targeted by the sgRNA under consideration; 
            }
followed by the columns containing the sgRNAs' counts for the controls and columns for library trasfected samples.
The argument is ignored if \code{Dframe} is not NULL.
}
\item{Dframe}{
        A data frame containing the raw sgRNA counts (usable as alternative to providing the path to a tsv file, i.e. previous argument). This must have one row per sgRNA and the following columns/headers: 
        \itemize{
            \item sgRNA: containing alphanumerical identifiers of the sgRNA under consideration;
            \item gene: containing HGNC symbols of the genes targeted by the sgRNA under consideration; 
            }
followed by the columns containing the sgRNAs' counts for the controls and columns for library trasfected samples.
If set to its default NULL value, then the function will try to load and use the file specified in \code{filename}.
}
  \item{display}{
    A logic value specifying whether figures containing boxplots with the count values pre/post normalisation and log fold-changes should be visualised (TRUE, by default). 
}
  \item{saveToFig}{
    A logic value specifying whether figures containing boxplots with the count values pre/post normalisation and log fold-changes should be saved as pdf files (FALSE, by default). Setting this parameter to TRUE overrides the value of the display parameter.
}
  \item{outdir}{
        Path of the directory where the normalised sgRNAs' counts and the log fold changes, as well as the pdf files (if the parameter saveToFig is set to TRUE), must be saved.
}
  \item{min_reads}{
    This parameter defines a filter threshold value for sgRNAs, based on their average counts in the control sample. Specifically, it indicates the minimal number of counts that each individual sgRNA needs to have in the controls (on average) in order to be included in the output.
}
  \item{EXPname}{
  A string specifying the name of the experiment. This will be used to compose main title of the generated figures and file names.
}
  \item{libraryAnnotation}{
  A data frame containing the sgRNA annotations, with a named row for each sgRNA, and columns for targeted genes, genomic coordinates and possibly other informations. This should be formatted as the \code{\link{KY_Library_v1.0}} data object containing the annotation of the sgRNA library presented in [2].
}
  \item{ncontrols}{
  A numerical value indicating the number of control replicates (therefore columns to be considered as control counts after the first two, in the inputted tsv file).
}
  \item{method}{
  A string specifying the normalisation method: 'ScalingByTotalReads' for scaling samples by total numbers or reads, 'MedRatios' to use the median of ratios method [1], or a gene name for scaling samples by total number of reads of the guides targeting that gene.
}
}
\value{
    A list containing two data frames: for the normalised sgRNAs' counts (norm_counts) and the sgRNAs' log fold changes (logFCs) respectively. First two columns in these data frames contain sgRNAs' identifiers and HGNC symbols of targete gene, respectively.
}
\references{
    [1] Anders S, Huber W. \emph{Differential expression analysis for sequence count data}. Genome Biol. 2010, 11: R106 \cr
    [2] Tzelepis K, Koike-Yusa H, De Braekeleer E, et al \emph{A CRISPR dropout screen identifies genetic vulnerabilities and therapeutic targets in acute myeloid leukaemia}. Cell Reports 2016 Oct 18;17(4):1193-1205\cr
}
\author{
Francesco Iorio (francesco.iorio@fht.org)
}
\seealso{
\code{\link{KY_Library_v1.0}}
}
\examples{\dontrun{
    ## loading sgRNA library annotation
    data(KY_Library_v1.0)
    
    ## derive path for an example dataset
    fn<-paste(system.file('extdata', package = 'CRISPRcleanR'),'/HT-29_counts.tsv',sep='')
    
    ## sgRNAs' normalisation and computation of log fold-changes
    normANDfcs<-ccr.NormfoldChanges(fn,
                                    min_reads=30,
                                    EXPname='Example',
                                    libraryAnnotation=KY_Library_v1.0)
    
    ## inspecting first 5 entries of the data frames containing the
    ## normalised counts and the log fold-changes
    head(normANDfcs$norm_counts)
    head(normANDfcs$logFCs)
}}
\keyword{analysis}% use one of  RShowDoc("KEYWORDS")

