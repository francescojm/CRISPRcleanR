\name{ccr.PlainTsvFile}
\alias{ccr.PlainTsvFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Saving a sgRNA counts' object in plain tsv file
}
\description{
This function takes in input a sgRNA counts' object, as outputted (for example) by the
\code{\link{ccr.NormfoldChanges}} function 
and saves it as plaing tab delimited text file (which can be processed by MAGeCK [1]).
}
\usage{
ccr.PlainTsvFile(sgRNA_count_object,
                 fprefix = "",
                 path = "./")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sgRNA_count_object}{
    sgRNA counts data object.  
}
  \item{fprefix}{
    A string specifying a name prefix of the tsv file which will contain the inputted sgRNA counts data object.
}
  \item{path}{
    A string specifying the location where the tsv file will be saved.
}
}
\value{
    A string specifying the complete path of the saves tsv file.
}
\references{
[1] Li, W., Xu, H., Xiao, T., Cong, L., Love, M. I., Zhang, F., et al. (2014). MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biology, 15(12), 554. [2] Hart, T., & Moffat, J. (2016). BAGEL: a computational framework for identifying essential genes from pooled library screens. BMC Bioinformatics, 17(1), 164.
}
\author{
Francesco Iorio (francesco.iorio@fht.org)
}
\seealso{
\code{\link{ccr.NormfoldChanges}}
}
\examples{\dontshow{setwd(tempdir())}
## Loading sgRNA library annotation file
data(KY_Library_v1.0)

## Deriving the path of the file with the example dataset,
## from the mutagenesis of the EPLC-272H colorectal cancer cell line
fn<-paste(system.file('extdata', package = 'CRISPRcleanR'),
          '/EPLC-272H_counts.tsv',sep='')

## Loading, median-normalizing and computing fold-changes for the example dataset
normANDfcs<-ccr.NormfoldChanges(fn,min_reads=30,
                                EXPname='EPLC-272H',
                                libraryAnnotation = KY_Library_v1.0,
                                display=FALSE)
                                
## saving median-normalised sgRNA counts' as a plain tsv file in ./EPLC-272H_sgRNA_count.tsv
uncorrected_fn<-ccr.PlainTsvFile(sgRNA_count_object = normANDfcs$norm_counts,fprefix = 'EPLC-272H')

uncorrected_fn
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utils}% use one of  RShowDoc("KEYWORDS")
