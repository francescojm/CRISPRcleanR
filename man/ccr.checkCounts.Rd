\name{ccr.checkCounts}
\alias{ccr.checkCounts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Check consistency between library annotation and count files.
}
\description{
 This function takes as input the sgRNA library annotation and the counts matrix to ensure that they are consistent.
  This a utility function that runs automatically as part of the \link{ccr.AnalysisPipeline}.
 \cr
}
\usage{
  ccr.checkCounts(
    counts,
    libraryAnnotation,
    ncontrols = 1,
    min_reads = 30
  ) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{counts}{
    A data frame containing the raw sgRNA counts (usable as an alternative to providing the path to a tsv file, i.e. previous argument). This must have one row per sgRNA and the following columns/headers: 
    \itemize{
        \item sgRNA: containing alphanumerical identifiers of the sgRNA under consideration;
        \item gene: containing HGNC symbols of the genes targeted by the sgRNA under consideration; 
        }
    followed by the columns containing the sgRNAs' counts for the controls and columns for library trasfected samples.
  }
  \item{libraryAnnotation}{
    A data frame containing a sgRNAs library. This data frame must include one named row per each sgRNA and the at least following mandatory columns/headers:
    \itemize{
    \item CODE: the unique ID of the sgRNA;
    \item GENES: the gene symbol related to the sgRNA;
    \item seq: the nucleotidic sequnce of the sgRNA without PAM
    }
    All the built-in libraries included in the package are already compliant with this structure.
  }
  \item{ncontrols}{
    A numerical value used by the \link{ccr.NormfoldChanges} indicating the number of control replicates (therefore columns to be considered as control counts after the first two, in the inputted tsv file). 1 by default.
    The parameter will not be considered when the input are FASTQ / BAM files. In this case, the counts obtained by the files listed in files_FASTQ_controls / files_BAM_controls parameters will be used as controls.
  }
  \item{min_reads}{
   A numerical value used by the \link{ccr.NormfoldChanges} to define a filter threshold value for sgRNAs, based on their average counts in the control sample. Specifically, it indicates the minimal number of counts that each individual sgRNA needs to have in the controls (on average) in order to be included in the output. 30 by default.
  }
}
\value{
  A boolean value equal to TRUE if the function ends without errors.
}
\author{
Paolo Cremaschi (paolo.crmeaschi@fht.org)
}
\seealso{
\code{
  \link{ccr.AnalysisPipeline}
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{analysis}% use one of  RShowDoc("KEYWORDS")

