\name{ccr.logFCs2chromPos}
\alias{ccr.logFCs2chromPos}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genomic sorting of sgRNAs' log fold changes.
}
\description{
This function maps genome-wide sgRNAs' log fold changes (averaged across replicates) on the genome and returns them sorted according to the position of their targeted region on the chromosomes.
}
\usage{
ccr.logFCs2chromPos(foldchanges, libraryAnnotation)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{foldchanges}{
    A data frame containing genome-wide sgRNAs' log fold changes, one column per library transfection replicate, with first and second column
    containing the sgRNAs' identifiers and the HGNC symbols of the targeted genes, respectively. This can be generated from raw count files using the \code{\link{ccr.NormfoldChanges}} function.
}
  \item{libraryAnnotation}{
    A data frame containing the sgRNAs' genome-wide annotations with at least a named row for each of the sgRNAs included in the foldchanges data frame provided in input.
    The following columns/headers should be present in this data frame (additional columns will be ignored):
    \itemize{
            \item GENES: string vector containing the HGNC symbols of the genes targeted by the sgRNA under consideration;
            \item EXONE: string vector containing the gene exon targeted by the sgRNA under consideration (these should include the prefix "ex" followed by the exone number);
            \item CHRM: string vector the chromosome of the gene targeted by the sgRNA under consideration (X and Y chromosome should be specified as "X" and "Y");
            \item STRAND: string vector containing the strand targeted by the sgRNA under consideration ("+" or "-");
            \item STARTpop: numeric vector containing the genomic coordinate of the starting position of the region targeted by the sgRNA under consideration;
            \item ENDpos: numeric vector containing the genomic coordinate of the ending position of the region targeted by the sgRNA under consideration;
            }
    Additiol columns can be optionally included and will be ignored by this function.
    The annation for the genome-wide sgRNA library presented in [1] is included in the \code{\link{KY_Library_v1.0}} data object, formatted as described above.
}
}
\value{
    A data frame with a named row per each sgRNA and the following columns/headers:
    \itemize{
            \item \code{CHR}: the chromosome where the gene targeted by the sgRNA under consideration resides;
            \item \code{startp}: the genomic coordinate of the starting position of the region targeted by the sgRNA under consideration;
            \item \code{endp}: the genomic coordinate of the ending position of the region targeted by the sgRNA under consideration;
            \item \code{avgFC}: the log fold change of the sgRNA averaged across replicates;
            \item \code{BP}: the genomic coordinate of the sgRNA defined as STARTpos+(ENDpos-STARTpos)/2.
            }
}
\references{
[1] Tzelepis K, Koike-Yusa H, De Braekeleer E, Li Y, Metzakopian E, Dovey OM, Mupo A, Grinkevich V, Li M, Mazan M, Gozdecka M, Onishi S, Cooper J, Patel M, McKerrell T, Chen B, Domingues AF, Gallipoli P, Teichmann S, Ponstingl H, McDermott U, Saez-Rodriguez J, Huntly BJP, Iorio F, Pina C, Vassiliou GS, Yusa K. \emph{A CRISPR dropout screen identifies genetic vulnerabilities and therapeutic targets in acute myeloid leukaemia}.
Cell Reports 2016 Oct 18;17(4):1193-1205
}
\author{
Francesco Iorio (francesco.iorio@fht.org)
}
\seealso{
    \code{\link{ccr.NormfoldChanges}}, \code{\link{KY_Library_v1.0}}
}
\examples{\dontrun{
    data(KY_Library_v1.0)
    fn<-paste(system.file('extdata', package = 'CRISPRcleanR'),'/A2058_counts.tsv',sep='')
    normANDfcs<-ccr.NormfoldChanges(fn,min_reads=30,
                                    EXPname='Example',
                                    libraryAnnotation=KY_Library_v1.0)
    mappedLogFCs<-ccr.logFCs2chromPos(normANDfcs$logFCs,KY_Library_v1.0)
    head(mappedLogFCs)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{analysis}% use one of  RShowDoc("KEYWORDS")

