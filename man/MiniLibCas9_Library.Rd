\name{MiniLibCas9_Library}
\alias{MiniLibCas9_Library}
\docType{data}
\title{
Genome-wide annotation for the MiniLibCas9 sgRNA library
}
\description{
A data frame with a named row for each sgRNA of the MiniLibCas9 sgRNA library [1] including annotations such as targeted genes, and genomic coordinates.
}
\usage{data("MiniLibCas9_Library")}
\format{
  A data frame with 37701 observations on the following variables (among others).
  \describe{
  \item{\code{CODE}}{alphanumerical identifier of the sgRNAs;}
  \item{\code{GENES}}{targeted gene;}
  \item{\code{STARTpos}}{starting genomic coordinate of the targeted genomic region (numeric);}
  \item{\code{STRAND}}{targeted DNA strand ('+' or '-')}
  \item{\code{CHRM}}{chromosome of where the targeted region resides (string)}
  \item{\code{ENDpos}}{ending genomic coordinate of the targeted genomic region (numeric).}
    }
}
\source{
https://github.com/EmanuelGoncalves/crispy/blob/master/notebooks/minlib/libraries/MinLibCas9.csv.gz
}
\references{
[1] Goncalves E, Thomas M, Behan FM, Picco G, Pacini C, Allen F, Parry-Smith D, et al. 2019.
Minimal Genome-Wide Human CRISPR-Cas9 Library. BioRxiv, January, 848895. https://doi.org/10.1101/848895
}
\examples{\dontrun{
## Loading sgRNA MiniLibCas9 library annotation file
data(MiniLibCas9_Library)
## Visualising first entries
head(MiniLibCas9_Library)

## Deriving the path of an example count file
## from screening the HT-29 cell line with the Brunello library
## [1]

fn<-paste(system.file('extdata', package = 'CRISPRcleanR'),
          '/HT29-MiniLibCas9_counts.tsv',sep='')

expName<-'HT29-MiniLibCas9'

## Loading, median-normalizing and computing fold-changes
normANDfcs<-
    ccr.NormfoldChanges(filename = fn,
                        display = TRUE,
                        min_reads = 30,
                        EXPname = expName,
                        libraryAnnotation = MiniLibCas9_Library)

## Genome-sorting the fold changes
gwSortedFCs<-
    ccr.logFCs2chromPos(foldchanges = normANDfcs$logFCs,
                        libraryAnnotation = MiniLibCas9_Library)

## Identifying and correcting biased sgRNAs' fold changes
correctedFCs_and_segments<-
    ccr.GWclean(gwSortedFCs = gwSortedFCs,
                display=TRUE,
                label=expName)
}}
\keyword{datasets}
