\name{ccr.getLibrary}
\alias{ccr.getLibrary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Convert library inputs in a standard library annotaion.
}
\description{
 This function takes as input a string identifying one of the sgRNA built-in library or a file name (comprehensive of the full path) to create a data frame with the standard format of the annotaion library requested by CRISPRcleanR.
  This is a utility function that runs automatically as part of the \link{ccr.AnalysisPipeline}.
 \cr
}
\usage{
  ccr.getLibrary(
    library_builtin,
    library_file,
    verbose = FALSE
  ) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{library_builtin}{
    A string containing the name of one of the libraries whose annotation is available in CRISPRcleanR.
  }
  \item{library_file}{
    A string specifying the path to a file containing the sgRNA annotations in TXT / TSV format, with columns for sgRNA ID ("CODE"), targeted genes ("GENES"), genomic coordinates, and possibly other information. This should be formatted as the \code{\link{KY_Library_v1.0}} data object containing the annotation of the sgRNA library presented in [1].
    If FASTQ files are used as input the annotation must include also a column "seq" with the nucleotidic sequence of the sgRNAs.
    This argument is ignored if a built-in library is specified.
  }
  \item{verbose}{
    Technical parameter undocumented. Used for the integration in web architecture.
  }
}

\value{
  A data frame compliant with the annotation library format requested by CRIPRcleanR.
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

