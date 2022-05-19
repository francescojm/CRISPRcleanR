\name{ccr.geneMeanFCs}
\alias{ccr.geneMeanFCs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Gene level log fold changes
}
\description{
This functions computes gene level log fold changes based on average log fold changes of targeting sgRNAs 
}
\usage{
ccr.geneMeanFCs(sgRNA_FCprofile, libraryAnnotation)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sgRNA_FCprofile}{
    A named numerical vector containing the sgRNAs' log fold-changes, with names corresponding to sgRNAs identifiers.
}
  \item{libraryAnnotation}{
    A data frame containing the sgRNA library annotation (with same format of \code{\link{KY_Library_v1.0}}).
}
}
\value{
    A numerical vector containing gene average log fold-changes, with corresponding HGNC symbols as names.
}
\author{
    Francesco Iorio (francesco.iorio@fht.org)
}
\seealso{
\code{\link{KY_Library_v1.0}}
}
\examples{
## loading corrected sgRNAs log fold-changes and segment annotations for
## an example cell line (EPLC-272H)
data(EPLC.272HcorrectedFCs)


## loading sgRNA library annotation
data(KY_Library_v1.0)

## storing sgRNA log fold-changes in a named vector
FCs<-EPLC.272HcorrectedFCs$corrected_logFCs$avgFC
names(FCs)<-rownames(EPLC.272HcorrectedFCs$corrected_logFCs)

## computing gene level log fold-changes
geneFCs<-ccr.geneMeanFCs(FCs,KY_Library_v1.0)

head(geneFCs)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utils}% use one of  RShowDoc("KEYWORDS")
