\name{HT.29correctedFCs}
\alias{HT.29correctedFCs}
\docType{data}
\title{
    CRISPRcleanR corrected data for an example cell line
}
\description{
This list contains corrected sgRNAs log fold-changes and segment annotations for an example cell line (HT-29), obtained using the \code{ccr.GWclean} function, as detailed in its reference manual entry \code{\link{ccr.GWclean}}.
}
\usage{data("HT.29correctedFCs")}
\format{
   A list containing two data frames and a vector of strings. The first data frame (corrected_logFCs) contains a named row per each sgRNA and the following columns/header:
        \itemize{
            \item \code{CHR}: the chromosome of the gene targeted by the sgRNA under consideration;
            \item \code{startp}: the genomic coordinate of the starting position of the region targeted by the sgRNA under consideration;
            \item \code{endp}: the genomic coordinate of the ending position of the region targeted by the sgRNA under consideration;
            \item \code{genes}: the HGNC symbol of the gene targeted by the sgRNA under consideration;
            \item \code{avgFC}: the log fold change of the sgRNA averaged across replicates;
            \item \code{correction}: the type of correction: 1 = increased log fold change, -1 = decreased log fold change. 0 indicates no correction;
            \item \code{correctedFC}: the corrected log fold change of the sgRNA
            }.\cr
    The second data frame (segments) contains the identified region of estimated equal log fold changes (one region per row) and the following
    columns/headers:
        \itemize{
            \item \code{CHR}: the chromosome of the gene targeted by the sgRNA under consideration;
            \item \code{startp}: the genomic coordinate of the starting position of the region targeted by the sgRNA under consideration;
            \item \code{endp}: the genomic coordinate of the ending position of the region targeted by the sgRNA under consideration;
            \item \code{genes}: the HGNC symbol of the gene targeted by the sgRNA under consideration;
            \item \code{avgFC}: the log fold change of the sgRNA averaged across replicates;
            \item \code{correction}: the type of correction: 1 = increased log fold change, -1 = decreased log fold change. 0 indicates no correction;
            \item \code{correctedFC}: the corrected log fold change of the sgRNA}.\cr
    The second data frame (segments) contains the identified region of estimated equal log fold changes (one region per row) and the following
    columns/headers:
        \itemize{
            \item \code{CHR}: the chromosome of the region under consideration;
            \item \code{startp}: the genomic coordinate of the starting position of the region under consideration;
            \item \code{endp}: the genomic coordinate of the ending position of the region under consideration;
            \item \code{n.sgRNAs}: the number of sgRNAs targeting sequences in the region under consideration;
            \item \code{avg.logFC}: the average log fold change of the sgRNAs in the region;
            \item \code{guideIdx}: the indexes range of the sgRNAs targeting the region under consideration as they appear in the gwSortedF Cs provided in input.
            }
    The string of vectors (SORTED_sgRNAs) contains the sgRNAs' identifiers in the same order as they are reported in the gwSortedFCs data frame inputted to the \code{ccr.GWclean} function.

}
\examples{
data(HT.29correctedFCs)
head(HT.29correctedFCs$corrected_logFCs)
head(HT.29correctedFCs$segments)
head(HT.29correctedFCs$SORTED_sgRNAs)
}
\keyword{datasets}
