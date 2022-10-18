\name{ccr.CreateLibraryIndex}
\alias{ccr.CreateLibraryIndex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Library index generation based on the sgRNA library to allow FASTQ files alignment.
}
\description{
 This function takes in input a data frame with the structure of a sgRNA library (e.g. KY) to create a library index file used for the alignment of the FASTQ files. It uses the "buildindex" function of the Rsubread [1].
  The data frame must include a "seq" field containing the nucleotidic sequences of the sgRNA that will be used to create the library.
  The function creates in the current folder a set of files related to the library index. All the files will include the EXPname label and will be prefixed with the "Library_" tag.
  Moreover, the function will create an additional FA file reporting all the sgRNA sequences in the FASTA format and a TXT file that can be used for the alignment using MAGeCK [2]. 
  The sgRNA library includes sequences that are common to more than one sgRNA to avoid multiple alignments of the same read targeting those sequences. By default, the function keeps only the first of the instances.
  Since the presence of a sgRNA targeting multiple genes might represent an undesirable source of noise, the function allows the user to completely exclude all the sgRNAs whose sequences appear more than once in the library. \cr
}
\usage{
  ccr.CreateLibraryIndex(
    libraryAnnotation,
    duplicatedSeq = c("keep", "exclude"),
    EXPname = "",
    indexMemory = 2000,
    overwrite = FALSE
  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{libraryAnnotation}{
    A data frame containing a sgRNAs library. This data frame must include one named row per each sgRNA and the at least following mandatory columns/headers:
    \itemize{
    \item CODE: the unique ID of the sgRNA;
    \item GENES: the gene symbol related to the sgRNA;
    \item seq: the nucleotidic sequence of the sgRNA without PAM.
    }
    All the built-in libraries included in the package are already compliant with this structure.
  }
  \item{duplicatedSeq}{
    A string defining the strategy to deal with the duplicated sequences that might be present in the library. 
    The possible options are "keep" (the default) or "exclude". The "keep" option will maintain the first occurrence of the duplicated sequences, while the "exclude" option will remove all the sgRNA whose nucleotidic sequences occur more than once.
  }
  \item{EXPname}{
    A string specifying the name of the experiment. This will be included in all the library files created by the function.
  }
  \item{indexMemory}{
    A numeric value specifying the amount of memory (in megabytes) used for storing the index during read mapping. 2000 MB by default.
  }
  \item{overwrite}{
    A boolean value specifying whether the index files should be overwritten, in case theyalready exist. FALSE by default.
  }
}

\value{
  The name of the library file.
}
\references{
 
  [1] Liao, Y., Smyth, G.K., Shi, W. (2019). The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads. Nucleic Acids Research, 47, e47 DOI:10.1093/nar/gkz114
    
  [2] Li, W., Xu, H., Xiao, T., Cong, L., Love, M. I., Zhang, F., et al. (2014). MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biology, 15(12), 554. [2] Hart, T., & Moffat, J. (2016). BAGEL: a computational framework for identifying essential genes from pooled library screens. BMC Bioinformatics, 17(1), 164.
    
}
\author{
Paolo Cremaschi (paolo.crmeaschi@fht.org)
}
\seealso{
\code{\link{ccr.FASTQ2counts}}
}
\examples{\dontrun{
## Loading sgRNA library annotation file
data(KY_Library_v1.0)

## Create the index file based on the KY 1.0 libary
fn <- ccr.CreateLibraryIndex(KY_Library_v1.0)

}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{analysis}% use one of  RShowDoc("KEYWORDS")

