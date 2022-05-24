\name{ccr.get.gdsc1000.AMPgenes}
\alias{ccr.get.gdsc1000.AMPgenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Copy number amplified genes in a given cell line from the GDSC1000
}
\description{
This function takes in input the name (or the COSMIC identifier [1]) of a cell line included in the GDSC1000 project [2] and it identifies the genes that are copy number amplified (according to a user defined minimal copy number value) in that cell line, using gene level copy number data from the Genomics of Drug Sensitivity in 1,000 Cancer Cell lines (GDSC1000) [2].
}
\usage{
ccr.get.gdsc1000.AMPgenes(cellLine, minCN = 8, exact = FALSE,
                          GDSC.geneLevCNA=NULL,GDSC.CL_annotation=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cellLine}{
    A string specifying the name of a cell line (or a COSMIC identifier [1]);
}
  \item{minCN}{
Lower threshold for the minimum copy number of any genomic segment containing coding sequence of a gene in order for it to be considered as copy number amplified.
}
  \item{exact}{
If TRUE, then those genes for which any genomic segment containing coding sequence has a minimum copy number equal to \code{minCN} are considered as copy number amplified.
  }
  \item{GDSC.geneLevCNA}{
    Genome-wide copy number data with the same format of \code{\link{GDSC.geneLevCNA}}. This can be assembled
    from the xls sheet specified in the source section [a] (containing data for the GDSC1000 cell lines). If NULL, then this function uses the data in the built in \code{\link{GDSC.geneLevCNA}} data frame, containing data derived from [a] for 15 cell lines used in [3] to assess the performances of CRISPRcleanR.
  }
  \item{GDSC.CL_annotation}{
    Cell lines annotation dataframe with the same structure of the \code{\link{GDSC.CL_annotation}}. If NULL then the \code{\link{GDSC.CL_annotation}} is used.
  }
}
\value{
 A data frame, containing one row for each copy number amplified gene with the following columns:
    \item{Gene}{
        HGNC symbol of the gene;
    }
    \item{minCN}{
        Minimum copy number of any genomic segment containing coding sequence of the gene in the
        cell line under consideration.
    }
}
\references{
[1] Forbes SA, Beare D, Boutselakis H, et al. \emph{COSMIC: somatic cancer genetics at high-resolution}
Nucleic Acids Research, Volume 45, Issue D1, 4 January 2017, Pages D777-D783, \cr

[2] Iorio F, Knijnenburg TA, Vis DJ, Bignell GR, Menden MP, et al. \emph{A landscape of pharmacogenomic interactions in cancer} Cell 2016 Jul 28;166(3):740-54 \cr

[3] Iorio, F., Behan, F. M., Goncalves, E., Beaver, C., Ansari, R., Pooley, R., et al. (n.d.). Unsupervised correction of gene-independent cell responses to CRISPR-Cas9 targeting. \cr
http://doi.org/10.1101/228189
}
\source{
[a] \url{ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/Gene_level_CN.xlsx}.\cr
}
\author{
Francesco Iorio (francesco.iorio@fht.org)
}
\seealso{
    \code{\link{ccr.get.CCLEgisticSets}}
}
\examples{\dontshow{setwd(tempdir())}
CNAgenes<-
    ccr.get.gdsc1000.AMPgenes('HT-29')
head(CNAgenes)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utils}% use one of  RShowDoc("KEYWORDS")
