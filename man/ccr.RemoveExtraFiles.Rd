\name{ccr.RemoveExtraFiles}
\alias{ccr.RemoveExtraFiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Clean intermediate files after pipeline execution.
}
\description{
 This function removes from the output data folder intermediate files created by the analysis pipeline run.
  This is a utility function that runs automatically as part of the \link{ccr.AnalysisPipeline}.
 \cr
}
\usage{
  ccr.RemoveExtraFiles(
    is_web = FALSE,
    file_counts = NULL,
    files_FASTQ_controls = NULL,
    files_FASTQ_samples = NULL,
    files_BAM_controls = NULL,
    files_BAM_samples = NULL,
    outdir_data = NULL
  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{is_web}{
    Boolen flag that indicates if the piepline run was part of a web application backend process.
  }
  \item{file_counts}{
    A string specifying the path of a tsv file containing the raw sgRNA counts.
  }
  \item{files_FASTQ_controls}{
    List of FASTQ files used to generate the counts for the control samples. Each file name should include the path. 
    The argument must be NULL if counts / BAM files are specified as input.
  }
  \item{files_FASTQ_samples}{
    List of FASTQ files used to generate the counts for the samples. Each file name should include the path. 
    The argument must be NULL if counts / BAM files are specified as input.
  }
  \item{files_BAM_controls}{
    List of BAM files used to generate the counts for the control samples. Each file name should include the path. 
    The argument must be NULL if counts / FASTQ files are specified as input.
  }
  \item{files_BAM_samples}{
    List of BAM files used to generate the counts for the samples. Each file name should include the path. 
    The argument must be NULL if counts / FASTQ files are specified as input.
  }
  \item{outdir_data}{
    A string specifying folder where all the results data files are stored. 
  }
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

