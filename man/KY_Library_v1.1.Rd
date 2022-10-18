\name{KY_Library_v1.1}
\alias{KY_Library_v1.1}
\docType{data}
\title{
    Genome-wide annotation for the Sanger sgRNA Library v1.1
}
\description{
A data frame with a named row for each sgRNA of the updated Sanger sgRNA library presented in [1] including annotations such as targeted genes, and genomic coordinates.
}
\usage{data(KY_Library_v1.1)}
\format{
  A a row named data frame with 90709 observations (one for each sgRNA) of the following 7 variables.
  \describe{
    \item{\code{CODE}}{alphanumerical identifier of the sgRNAs;}
    \item{\code{GENES}}{targeted gene;}
    \item{\code{EXONE}}{exone of the targeted genomic region (string with 'ex' prefix followed by the exone number);}
    \item{\code{CHRM}}{chromosome of where the targeted region resides (string)}
    \item{\code{STRAND}}{targeted DNA strand ('+' or '-')}
    \item{\code{STARTpos}}{starting genomic coordinate of the targeted genomic region (numeric);}
    \item{\code{ENDpos}}{ending genomic coordinate of the targeted genomic region (numeric).}
    \item{\code{seq}}{nucelotidic sequence of the sgRNAs without the PAM. (string).}
  }
}
\references{
[1] Tzelepis K, Koike-Yusa H, De Braekeleer E, Li Y, Metzakopian E, Dovey OM, Mupo A, Grinkevich V, Li M, Mazan M, Gozdecka M, Onishi S, Cooper J, Patel M, McKerrell T, Chen B, Domingues AF, Gallipoli P, Teichmann S, Ponstingl H, McDermott U, Saez-Rodriguez J, Huntly BJP, Iorio F, Pina C, Vassiliou GS, Yusa K. \emph{A CRISPR dropout screen identifies genetic vulnerabilities and therapeutic targets in acute myeloid leukaemia}.
Cell Reports 2016 Oct 18;17(4):1193-1205
}
\examples{
data(KY_Library_v1.1)
head(KY_Library_v1.1)
}
\keyword{datasets}
\keyword{Supported sgRNA libraries}
