\name{Whitehead_Library}
\alias{Whitehead_Library}
\docType{data}
\title{
Genome-wide annotation for the Whitehead sgRNA library
}
\description{
A data frame with a named row for each sgRNA of the Whitehead sgRNA library [1] including annotations such as targeted genes, and genomic coordinates.
}
\usage{data(Whitehead_Library)}
\format{
  A a row named data frame with 181131 observations of the following variables (among others)
  \describe{
    \item{\code{CODE}}{alphanumerical identifier of the sgRNAs;}
    \item{\code{GENES}}{targeted gene;}
    \item{\code{STARTpos}}{starting genomic coordinate of the targeted genomic region (numeric);}
    \item{\code{STRAND}}{targeted DNA strand ('sense' or 'antisense')}
    \item{\code{EXONE}}{exone of the targeted genomic region (exone number);}
    \item{\code{CHRM}}{chromosome of where the targeted region resides (string)}
    \item{\code{ENDpos}}{ending genomic coordinate of the targeted genomic region (numeric).}
  }
}
\references{
[1] Wang T, Birsoy K, Hughes NW, et al. Identification and characterization of essential genes in the human genome. Science. 2015;350(6264):1096-1101. doi:10.1126/science.aac7041
}
\source{
Discontinued by Addgene
}
\examples{
## Loading sgRNA Whitehead library annotation file
data(Whitehead_Library)
## Visualising first entries
head(Whitehead_Library)

## Deriving the path of an example count file
## from screening the HT-29 cell line with the Whitehead library
## [2]
fn<-paste(system.file('extdata', package = 'CRISPRcleanR'),
          '/HT-29-Whitehead_counts.tsv',sep='')

expName<-'HT29-Whitehead'

## Loading, median-normalizing and computing fold-changes
normANDfcs<-
    ccr.NormfoldChanges(filename = fn,
                        display = TRUE,
                        min_reads = 30,
                        EXPname = expName,
                        libraryAnnotation = Whitehead_Library)

## Genome-sorting the fold changes
gwSortedFCs<-
    ccr.logFCs2chromPos(foldchanges = normANDfcs$logFCs,
                        libraryAnnotation = Whitehead_Library)

## Identifying and correcting biased sgRNAs' fold changes
correctedFCs_and_segments<-
    ccr.GWclean(gwSortedFCs = gwSortedFCs,
                display=TRUE,
                label=expName)
}
\keyword{datasets}
\keyword{Supported sgRNA libraries}
