\name{ccr.BAM2counts}
\alias{ccr.BAM2counts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Raw count file extraction from BAM file list.
}
\description{
 This function takes as input a list of BAM files and reads them to generate a raw count matrix with one column for each BAM file in the list.
  Based on the GenomicAlignments package [1], the function requires that the BAMs were created by aligning the FASTQ files to a genome file in which the seqnames match the sgRNA ID.
  The same IDs should be available in the CODE column from the library annotation data frame supplied as input.
  The function allows the user to select which strand to read the data from. 
  By default, the function counts only the reads with the same sequence of the sgRNA. However, depending on the strategy used to create the BAM files, the user can select to match only the reads with nsequnces complementary to the sgRNAs or to count the reads without any strand-related filter. \cr
}
\usage{
  ccr.BAM2counts(
    BAMfileList,
    libraryAnnotation,
    maxMismatches = 0,
    strand = c("F", "R", "*"),
    EXPname = "",
    outdir = "./",
    export_counts = TRUE,
    overwrite = TRUE
  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BAMfileList}{
    List of BAM files used to generate the count file. Each file should include the path to the BAM. 
    If present, the name of each element of the list will be used as sample name in the count matrix.
  }
  \item{libraryAnnotation}{
    A data frame containing a sgRNAs library. This data frame must include one named row per each sgRNA and the at least following mandatory columns/headers:
    \itemize{
    \item CODE: the unique ID of the sgRNA;
    \item GENES: the gene symbol related to the sgRNA;
    \item seq: the nucleotidic sequence of the sgRNA without PAM.
    }
    All the built-in libraries included in the package are already compliant with this structure.
  }
  \item{maxMismatches}{
    Integer number containing the Ns allowed to count the reads. 
    The function will consider as valid only the reads with a number of matched bases equal or greater than the length of the sgRNA sequence, provided in the annotation library, minus the maxMismatches parameter.
  }
  \item{strand}{
   A string specifying the strand of the alinement used to count the reads. It accepts three different options: "F" to count only the reads equal to the sgRNA sequence (default); "R" to count conly the reads complementary to the sgRNA sequence; "*" to count all the reads without any strand filter. See the function description for details.
  }
  \item{EXPname}{
    A string specifying the name of the experiment. This will be used to create the raw count file if the export_counts option is set to TRUE.
  }
  \item{outdir}{
    A string specifying folder where the raw count file will be created if the export_counts option is set to TRUE.
  }
  \item{export_counts}{
    A boolean value specifying if the raw count matrix should also be exported in a TXT file. TRUE by default.
  }
  \item{overwrite}{
    A boolean value specifying if the raw count file will overwrite any file with the same name already present in the outdir path (only when export_counts option is set to TRUE).
  }
}

\value{
  A dataframe with the raw counts related to each sample. The first two columns in these data frame contain sgRNAs' identifiers and HGNC symbols of target gene, respectively.
}
\references{
 
  [1] Lawrence, M., Huber, W., Pagès, H., Aboyoun, P., Carlson, M., Gentleman, R., et al. (2013). Software for Computing and Annotating Genomic Ranges. PLoS Computational Biology, 10.1371/journal.pcbi.1003118
    
}
\author{
Paolo Cremaschi (paolo.crmeaschi@fht.org)
}
\seealso{
\code{\link{ccr.FASTQ2counts}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{analysis}% use one of  RShowDoc("KEYWORDS")

