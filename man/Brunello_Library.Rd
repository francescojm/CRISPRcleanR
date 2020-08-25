\name{Brunello_Library}
\alias{Brunello_Library}
\docType{data}
\title{
sgRNAs' genome-wide annotation for the Brunello sgRNA library
}
\description{
A data frame with a named row for each sgRNA of the Brunello sgRNA library [1] including annotations such as targeted genes, and genomic coordinates.
}
\usage{data("Brunello_Library")}
\format{
  A data frame with 76379 observations of the following variables (among others)
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
[1] Doench JG, Fusi N, Sullender M, Hegde M, Vaimberg EW, Donovan KF, et al.
Optimized sgRNA design to maximize activity and minimize off-target
effects of CRISPR-Cas9. Nat Biotechnol. 2016;34:184–91.
}
\source{
Addgene website (catalog number: 73179; file:
broadgpp-brunello-library-contents.txt, url:
https://www.addgene.org/static/cms/filer_public/8b/4c/8b4c89d9-eac1-44b2-bb2f-8fea95672705/broadgpp-brunello-library-contents.txt)
}
\examples{

## Loading sgRNA Brunello library annotation file
data(Brunello_Library)
## Visualising first entries
head(Brunello_Library)

## Deriving the paht of an example count file
## from screening the HT-29 cell line with the Brunello library
fn<-paste(system.file('extdata', package = 'CRISPRcleanR'),
          '/HT29-Brunello_counts.tsv',sep='')

expName<-'HT29-Brunello'

## Loading, median-normalizing and computing fold-changes
normANDfcs<-
    ccr.NormfoldChanges(filename = fn,
                        display = TRUE,
                        min_reads = 30,
                        EXPname = expName,
                        libraryAnnotation = Brunello_Library)

## Genome-sorting the fold changes
gwSortedFCs<-
    ccr.logFCs2chromPos(foldchanges = normANDfcs$logFCs,
                        libraryAnnotation = Brunello_Library)

## Identifying and correcting biased sgRNAs' fold changes
correctedFCs_and_segments<-
    ccr.GWclean(gwSortedFCs = gwSortedFCs,
                display=TRUE,
                label=expName)
}
\keyword{datasets}
