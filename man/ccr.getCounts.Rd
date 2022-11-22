\name{ccr.getCounts}
\alias{ccr.getCounts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Convert count inputs in a standard count matrix format.
}
\description{
 The pipeline takes as input raw count/matrix or lists of FASTQ/BAM files to return a data frame containing all samples counts.
  This is a utility function that runs automatically as part of the \link{ccr.AnalysisPipeline}.
 \cr
}
\usage{
  ccr.getCounts(
    file_counts,
    files_FASTQ_controls,
    files_FASTQ_samples,
    files_BAM_controls,
    files_BAM_samples,
    libraryAnnotation,
    maxMismatches,
    nTrim5,
    nTrim3,
    nthreads,
    nBestLocations,
    duplicatedSeq,
    indexMemory,
    strand,
    EXPname,
    outdir_data,
    outdir_bam,
    aligner,
    fastqc_plots,
    verbose
  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_counts}{
    A string specifying the path to a tsv file containing the raw sgRNA counts. The file must be tab delimited with one row per sgRNA and the following columns/headers: 
    \itemize{
      \item sgRNA: containing alphanumerical identifiers of the sgRNA under consideration;
      \item gene: containing HGNC symbols of the genes targeted by the sgRNA under consideration; 
    }
    followed by the columns containing the sgRNAs' counts for the controls and columns for library transfected samples.
    The argument must be NULL if FASTQ / BAM files are specified as input.
  }
  \item{files_FASTQ_controls}{
    List of FASTQ files used to generate the counts for the control samples. Each file name should include the path. 
    If present, the name of each element of the list will be used as sample name for the BAM files and in the count matrix.
    The argument must be NULL if counts/BAM files are specified as input.
  }
  \item{files_FASTQ_samples}{
    List of FASTQ files used to generate the counts for the samples. Each file name should include the path. 
    If present, the name of each element of the list will be used as sample name for the BAM files and in the count matrix.
    The argument must be NULL if counts/BAM files are specified as input.
  }
  \item{files_BAM_controls}{
    List of BAM files used to generate the counts for the control samples. Each file name should include the path. 
    If present, the name of each element of the list will be used as sample name in the count matrix.
    The argument must be NULL if counts/FASTQ files are specified as input.
  }
  \item{files_BAM_samples}{
    List of BAM files used to generate the counts for the samples. Each file name should include the path. 
    If present, the name of each element of the list will be used as sample name in the count matrix.
    The argument must be NULL if counts/FASTQ files are specified as input.
  }
  \item{libraryAnnotation}{
    A data frame containing a sgRNAs library. This data frame must include one named row per each sgRNA and the at least following mandatory columns/headers:
    \itemize{
    \item CODE: the unique ID of the sgRNA;
    \item GENES: the gene symbol related to the sgRNA;
    \item seq: the nucleotidic sequence of the sgRNA without PAM
    }
    All the built-in libraries included in the package are compliant with this structure.
  }
  \item{maxMismatches}{
    Integer number containing the Ns allowed counting the reads. 
    The function will consider valid only the reads with a number of matched bases equal to or greater than the length of the sgRNA sequence, provided in the annotation library, minus the maxMismatches parameter.
    The parameter is ignored if the input are not FASTQ / BAM files.
  }
  \item{nTrim5}{
    Numeric value giving the number of bases trimmed off from 5' end of each read. 0 by default.
    If aligner = "Rsubreads" (the default ) the parameter accept only numeric values. 
    If MAGeCK is used the parameter can specify multiple trimming lengths, separated by comma (,); for example, "7,8" or can be set to "AUTO" to allow MAGeCK to determine the trimming length.
    The parameter is ignored if the input are not FASTQ files.
  }
  \item{nTrim3}{
    Numeric value giving the number of bases trimmed off from 3' end of each read. 0 by default.
    If aligner = "Rsubreads" (the default ) the parameter accept only numeric values. 
    The parameter is ignored if the input are not FASTQ files or if MAGeCK is used as aligner.
  }
  \item{nBestLocations}{
    Numeric value specifying the maximal number of equally-best mapping locations that will be reported for a multi-mapping read (max 16). 2 by default.
    Different tags will be included to the aligned sequences to specify the number of alignments reported. Please refer to the Rsubread [1] user guide for a complete description.
    The parameter is ignored if the input are not FASTQ files or if MAGeCK is used as aligner.
  }
  \item{duplicatedSeq}{
    A string defining the strategy to deal with the duplicated sequences in the library index creation. See ccr.CreateLibraryIndex for details.
    The possible options are "keep" (the default) or "exclude". The "keep" option will maintain the first occurrence of the duplicated sequences while the "exclude" option will remove all the sgRNA whose nucleotidic sequences occur more than once.
    The parameter is ignored if the input are not FASTQ files.
  }
  \item{indexMemory}{
    A numeric value specifying the amount of memory (in megabytes) used for storing the index during read mapping. 2000 MB by default.
    The parameter is ignored if the input are not FASTQ files.
  }
  \item{strand}{
    A string specifying the strand of the alinement used to count the reads. It accepts three different options: "F" to count only the reads equal to the sgRNA sequence (default); "R" the read only the reads complementary to the sgRNA sequence; "*" all the reads without any strand filter. See the function description for details.
    The parameter is ignored if the input are not FASTQ files.
  }
  \item{EXPname}{
    A string specifying the name of the experiment. This will be used as label in the output plots.
  }
  \item{outdir_data}{
    A string specifying the folder where all the results will be saved.  
  }
  \item{outdir_bam}{
    A string specifying the folder where all the BAM files created by the alignment process will be saved. 
    The parameter is not considered in the case of BAM/counts file input.
  }
  \item{aligner}{
    A string specifying the aligner used to map the reads to the library. 
    The possible options are "Rsubreads" (default) and "mageck" if the MAGeCK application [2] is installed.
    The parameter is ignored if the input are not FASTQ files.
  }
  \item{fastqc_plots}{
    A boolean value specifying if the QC plots for each FASTQ files will be generated during the alignment.
    The QC plots are created using the rqc function in the Rqc package.
    All the plots for each FASTQ file are collected in one HTML file named as the BAM file created in the alignment.
    The parameter is ignored if the input are not FASTQ files.
  }
  \item{verbose}{
    Technical parameter undocumented. Used for the integration in web architecture.
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

