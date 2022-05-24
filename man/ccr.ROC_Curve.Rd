\name{ccr.ROC_Curve}
\alias{ccr.ROC_Curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Classification performances of reference sets of genes (or sgRNAs) based on depletion log fold-changes
}
\description{
This functions computes Specificity/Sensitivity (or TNR/TPR, or ROC) curve, area under the ROC curve and (optionally) Recall (i.e. TPR) at fixed false discovery rate (computed as 1 - Precision (or Positive Predicted Value)) and corresponding log fold change threshold) when classifying reference sets of genes (or sgRNAs) based on their depletion log fold-changes
}
\usage{
ccr.ROC_Curve(FCsprofile,
                         positives,
                         negatives,
                         display = TRUE,
                         FDRth = NULL,
                         expName = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FCsprofile}{
A numerical vector containing gene average depletion log fold changes (or sgRNAs' depletion log fold changes) with names
corresponding to HGNC symbols (or sgRNAs' identifiers).
}
  \item{positives}{
A vector of strings containing a reference set of positive cases: HGNC symbols of essential genes or identifiers of their targeting sgRNAs. This must be a subset of \code{FCsprofile} names, disjointed from \code{negatives}.
}
  \item{negatives}{
A vector of strings containing a reference set of negative cases: HGNC symbols of essential genes or identifiers of their targeting sgRNAs. This must be a subset of \code{FCsprofile} names, disjointed from \code{positives}.
}
  \item{display}{
A logical parameter specifying if a plot containing the computed ROC curve with ROC indicators should be plotted (default = TRUE).
}
  \item{FDRth}{
  If different from NULL, will be a numerical value >=0 and <=1 specifying the false discovery rate threshold at which
  fixed recall will be computed. In this case, if the \code{display} parameter is TRUE, an orizontal dashed line will be   added to the plot at the resulting recall and its value will be visualised in the legend.
}
  \item{expName}{
  If different from NULL and \code{display} parameter is TRUE this parameter should be a string specifying the title of the plot with the computed ROC curve.
}
}
\value{
    A list containint three numerical variable \code{AUC}, \code{Recall}, and \code{sigthreshold} indicating the area under ROC curve and (if \code{FDRth} is not NULL) the recall at the specifying false discovery rate and the corresponding log fold change threshold (both equal to NULL, if \code{FDRth} is NULL), respectively. 
}
\author{
    Francesco Iorio (francesco.iorio@fht.org)
}
\seealso{
\code{\link{BAGEL_essential}}, \code{\link{BAGEL_nonEssential}},

\code{\link{ccr.genes2sgRNAs}}, \code{\link{ccr.VisDepAndSig}},

\code{\link{ccr.PrRc_Curve}}
}
\examples{\dontrun{
## loading corrected sgRNAs log fold-changes and segment annotations for an example
## cell line (EPLC-272H)
data(EPLC.272HcorrectedFCs)

## loading reference sets of essential and non-essential genes
data(BAGEL_essential)
data(BAGEL_nonEssential)

## loading library annotation
data(KY_Library_v1.0)

## storing sgRNA log fold-changes in a named vector
FCs<-EPLC.272HcorrectedFCs$corrected_logFCs$avgFC
names(FCs)<-rownames(EPLC.272HcorrectedFCs$corrected_logFCs)

## deriving sgRNAs targeting essential and non-essential genes (respectively)
BAGEL_essential_sgRNAs<-ccr.genes2sgRNAs(KY_Library_v1.0,BAGEL_essential)
BAGEL_nonEssential_sgRNAs<-ccr.genes2sgRNAs(KY_Library_v1.0,BAGEL_nonEssential)

## computing classification performances at the sgRNA level
ccr.ROC_Curve(FCs,BAGEL_essential_sgRNAs,BAGEL_nonEssential_sgRNAs)

## computing gene level log fold-changes
geneFCs<-ccr.geneMeanFCs(FCs,KY_Library_v1.0)

## computing classification performances at the sgRNA level, with Recall at 5\% FDR
ccr.ROC_Curve(geneFCs,BAGEL_essential,BAGEL_nonEssential,FDRth = 0.05)
}}
\keyword{Assessment and Visualisation}% use one of  RShowDoc("KEYWORDS")
