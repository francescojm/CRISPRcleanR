\name{ccr.FASTQ2counts}
\alias{ccr.FASTQ2counts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  FASTQ files alignment and raw count file extraction.
}
\description{
 This function take as input a list of FASTQ files and align them the sequences of the sgRNA library using by default the Rsubread [1] align function. The process generates a list of BAM files that are than used to extract a raw count file with one column for each FASTQ file in the list.
  If Rsubread is used for the alignment, the sorted BAM files will be generated in the outdir path.
  If MAGeCK [2] is installed, it can be specified as alternative aligner through its count function. In this case, some of the parameters might be unavailable and no BAM files will be created.
  The annotation data frame must include a "seq" field containing the nucleotidic sequences of the sgRNA that will be used to create the index library files for the alignment.
  By default the function allows the alignment of the reads in two different locations (nBestLocations = 2), while the duplications are removed in a second step based on the strand (strand = "F" of strand = "R" if the sequences in the library are complementary to the sgRNA sequences). This stragy allows the the identification of correct sgRNA in case of guides based on complementary sequnces.
  This alignment strtegy works optimally when no mismatchs (maxMismatches = 0) are allowed and the reads are trimmed to exact length of the sgRNA length.
  If the Ns are allowed in the alignment (maxMismatches > 0) or the reads after trimming are longer than the sgRNA sequences a more robust approach is based on the selection of only one locations (nBestLocations = 1) removing any filter based on the strand (strand = "*").
  As part of the alignment process the function will generate by default quality reports based on the Rqc package [3].
 \cr
}
\usage{
  ccr.FASTQ2counts(
    FASTQfileList,
    libraryAnnotation,
    maxMismatches = 0,
    nTrim5 = "0",
    nTrim3 = "0",
    nthreads = 1,
    nBestLocations = 2,
    strand = c("F", "R", "*"),
    indexMemory = 2000,
    duplicatedSeq = c("keep", "exclude"),
    EXPname = "",
    outdir = "./",
    aligner = c("Rsubreads", "mageck"),
    fastqc_plots = TRUE,
    export_counts = TRUE,
    overwrite = FALSE
  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FASTQfileList}{
    List of BAM files used to generate the count file. Each file should include the path to the FASTQ. 
    If present, the name of each element of the list will be used as sample name for the BAM files and in the count matrix.
  }
  \item{libraryAnnotation}{
    A data frame containing a sgRNAs library. This data frame must include one named row per each sgRNA and the at least following mandatory columns/headers:
    \itemize{
    \item CODE: the unique ID of the sgRNA;
    \item GENES: the gene symbol related to the sgRNA;
    \item seq: the nucleotic sequnce of the sgRNA without PAM
    }
    All the built-in libraries included in the package are already compliant with this structure.
  }
  \item{maxMismatches}{
    Integer number containing the Ns allowed to count the reads. 
    The function will consider as valid only the reads with a number of matched bases equal or greater than the length of the sgRNA sequence, provided in the annotation library, minus the maxMismatches parameter.
  }
  \item{nTrim5}{
    Numeric value giving the number of bases trimmed off from 5' end of each read. 0 by default.
    If aligner = "Rsubreads" (default), the parameter accept only numeric values.
    If MAGeCK is used, the parameter can specify multiple trimming lengths separated by comma (for example "7,8"), or can be set to "AUTO" to allow MAGeCK determine the trimming length.
  }
  \item{nTrim3}{
    Numeric value giving the number of bases trimmed off from 3' end of each read. 0 by default.
    If aligner = "Rsubreads" (default), the parameter accepts only numeric values. 
    If MAGeCK is used the parameter is ignored.
  }
  \item{nthreads}{
    Numeric value giving the number of threads used for mapping. 1 by default.
  }
  \item{nBestLocations}{
    Numeric value specifying the maximal number of equally-best mapping locations that will be reported for a multi-mapping read (max 16). 2 by default.
    Different tabs will be included to the aligned sequences to specify the number of alignments reported. Please refer to the Rsubread [1] user guide for a complelte description.
    If MAGeCK is used the parameter is ignored.
  }
  \item{strand}{
    A string specifying the strand of the alinement used to count the reads. It accepts three different options: "F" to count only the reads equal to the sgRNA sequence (default); "R" to count only the reads complementary to the sgRNA sequence; "*" to count all the reads without any strand filter. See the function description for details.
  }
  \item{indexMemory}{
    A numeric value specifying the amount of memory (in megabytes) used for storing the index during read mapping. 2000 MB by default.
  }
  \item{duplicatedSeq}{
    A string defining the strategy to deal with the duplicated sequences in the library index creation. See ccr.CreateLibraryIndex for details.
    The possible options are "keep" (the default) or "exclude". The "keep" option will maintain the first occurrence of the duplicates sequences while the "exclude" option will remove all the sgRNA whose nucleotidic sequences occur more than once.
    The parameter is ignored if the input are not FASTQ files.
  }
  \item{EXPname}{
    A string specifying the name of the experiment. This will be used to create the raw count file if the export_counts option is set to TRUE.
  }
  \item{outdir}{
    A string specifying the folder where the raw count file will be created if the export_counts option is set to TRUE.
  }
  \item{aligner}{
    A string specifying the aligner used to map the reads to the libary. 
    The possible options are "Rsubreads" (default) and "mageck" (if the MAGeCK application [2] is installed).
  }
  \item{fastqc_plots}{
    A boolean value specifying if the QC plots for each FASTQ files will be generated during the alignment.
    The QC plots are created using the rqc function in the Rqc package.
    All the plots for each FASTQ file are collected in one HTML file named as the BAM file created in the alignment.
  }
  \item{export_counts}{
    A boolean value specifying if the raw count matrix should also be exported in a TXT file. TRUE by default.
  }
  \item{overwrite}{
    A boolean value specifying if the files generated by the alignment will overwrite any file with the same name already present in the outdir path.
  }
}

\value{
  A data.frame containing the raw counts related to each sample. The first two columns in these data frame contain sgRNAs' identifiers and HGNC symbols of target gene, respectively.
}
\references{
 
  [1] Liao, Y., Smyth, G.K., Shi, W. (2019). The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads. Nucleic Acids Research, 47, e47 DOI:10.1093/nar/gkz114
    
  [2] Li, W., Xu, H., Xiao, T., Cong, L., Love, M. I., Zhang, F., et al. (2014). MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biology, 15(12), 554.

  [3] de Souza, W., Carvalho, B.S., Lopes-Cendes, I. (2018). Rqc: A Bioconductor Package for Quality Control of High-Throughput Sequencing Data. Journal of Statistical Software, Code Snippets, 87(2), 1–14. DIO:10.18637/jss.v087.c02

}
\author{
Paolo Cremaschi (paolo.crmeaschi@fht.org)
}
\seealso{
\code{
  \link{ccr.BAM2counts}
  \link{ccr.CreateLibraryIndex}
}}
\examples{\dontrun{
## Loading sgRNA library annotation file
data(KY_Library_v1.0)

## Define the list of FASTQ files to used for the alignment
fileList <- file.path(
  system.file("extdata", package = "CRISPRcleanR"),
  c("test_plasmid.fq.gz", "test_sample1.fq.gz", "test_sample2.fq.gz")
)

## Run the alignment and extract the raw counts
fn <- ccr.FASTQ2counts(
  FASTQfileList = fileList,
  libraryAnnotation = KY_Library_v1.0
)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{analysis}% use one of  RShowDoc("KEYWORDS")

