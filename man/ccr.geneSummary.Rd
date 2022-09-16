\name{ccr.geneSummary}
\alias{ccr.geneSummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gene level depletion summary
}
\description{
This function collapses single-guide RNAs (sgRNAs) depletion log fold-changes (logFCs)
on a targeted gene basis, by averaging (using the \code{ccr.geneMeanFCs}). In addition it computes also a logFC threshold \eqn{T} such that when considering as significantly depleted all the genes with a depletion logFC \eqn{< T}, the false discover rate
(FDR) of prior known non-essential genes is below a given threshold (specificed in input). Finally it calls
significantly depleted genes according to the computed threshold.\\
The significant threshold is computed using the \code{ccr.PrRc_Curve} function, employing a reference list of core-fitness essential genes and a list of non-essential genes, assembled from multiple RNAi studies used as classification template by the BAGEL algorithm to call gene depletion significance (included as the built-in objects \code{BAGEL_essential} and \code{BAGEL_nonEssential}) [1].
}
\usage{
ccr.geneSummary(sgRNA_FCprofile,
                libraryAnnotation,
                FDRth=0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sgRNA_FCprofile}{
    A named numerical vector containing the sgRNAs' log fold-changes, with names corresponding to sgRNAs identifiers.
    }
  \item{libraryAnnotation}{
    A data frame containing the sgRNA library annotation (with same format of \code{\link{KY_Library_v1.0}})
    }
  \item{FDRth}{
    The FDR threshold to consider in order to derive the significance threshold (FDR 5\% by default)
    }
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
    A data frame with gene symbols as row names and two columns: the first one indicating the gene depletion logFC and the second one including a boolean value specifying if the gene under consideration is significantly depleted at the indicated FDR level.
}
\references{
[1] BAGEL: a computational framework for identifying essential genes from pooled library screens. Traver Hart and Jason Moffat. BMC Bioinformatics, 2016 vol. 17 p. 164.
}
\author{
Francesco Iorio (francesco.iorio@fht.org)
}
\note{
%%  ~~further notes~~
}
\seealso{
\code{\link{KY_Library_v1.0}}, \code{\link{ccr.geneMeanFCs}},
\code{\link{ccr.PrRc_Curve}}, \code{\link{ccr.PrRc_Curve}},
\code{\link{BAGEL_essential}}, \code{\link{BAGEL_nonEssential}} 

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
    geneFCs<-ccr.geneSummary(FCs,KY_Library_v1.0)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
