\name{ccr.get.CCLEgisticSets}
\alias{ccr.get.CCLEgisticSets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CCLE gistic score gene sets
}
\description{
    This function splits all the genes into 5 classes (-2, -1, 0, +1 and +2) based on the CNA Gistic [1] score observed in a given cell line.
}
\usage{
ccr.get.CCLEgisticSets(cellLine,CCLE.gisticCNA=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cellLine}{A string specifying the name of a cell line (or a COSMIC identifier [2]);
}
  \item{CCLE.gisticCNA}{
    Genome-wide Gistic [1] scores quantifying copy number status across cell lines with the same format of
    \code{\link{CCLE.gisticCNA}}. If NULL then this function uses the \code{\link{CCLE.gisticCNA}} builtin data frame, containing data for 13 cell lines of the 15 used in [3] to assess the performances of CRISPRcleanR.
}
}
\value{
    A named list of vectors with the following fields:
    \item{gm2}{
        A vector of strings containing identifiers of sgRNAs targeting genes whit a Gistic score = -2 in the cell line under consideration;
    }
    \item{gm1}{
        A vector of strings containing identifiers of sgRNAs targeting genes whit a Gistic score = -1 in the cell line under consideration;
    }
    \item{gz}{
        A vector of strings containing identifiers of sgRNAs targeting genes whit a Gistic score = 0 in the cell line under consideration;
    }
    \item{gp1}{
        A vector of strings containing identifiers of sgRNAs targeting genes whit a Gistic score = +1 in the cell line under consideration;
    }
    \item{gp2}{
        A vector of strings containing identifiers of sgRNAs targeting genes whit a Gistic score = +2 in the cell line under consideration;
    }
    
}
\references{
[1] Mermel CH, Schumacher SE, Hill B, et al. \emph{GISTIC2.0 facilitates sensitive and confident localization of the targets of focal somatic copy-number alteration in human cancers.} Genome Biol. 2011;12(4):R41. doi: 10.1186/gb-2011-12-4-r41. \cr

[2] Forbes SA, Beare D, Boutselakis H, et al. \emph{COSMIC: somatic cancer genetics at high-resolution}
Nucleic Acids Research, Volume 45, Issue D1, 4 January 2017, Pages D777-D783, \cr

[3] Iorio, F., Behan, F. M., Goncalves, E., Beaver, C., Ansari, R., Pooley, R., et al. (n.d.). Unsupervised correction of gene-independent cell responses to CRISPR-Cas9 targeting. \cr
http://doi.org/10.1101/228189
}
\author{
Francesco Iorio (iorio@ebi.ac.uk)
}
\seealso{
\code{\link{ccr.get.gdsc1000.AMPgenes}}
}
\examples{
GS<-ccr.get.CCLEgisticSets('HT-29')

head(GS$gm2)
head(GS$gm1)
head(GS$gz)
head(GS$gp1)
head(GS$gp2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utils}% use one of  RShowDoc("KEYWORDS")
