\name{ccr.VisDepAndSig}
\alias{ccr.VisDepAndSig}
\title{
    Depletion profile visualisation with genes signatures superimposed and recall
}
\description{
This functions ranks the gene (or sgRNAs) log fold changes. Based on this it determines a log fold change threshold based on a user defined false discovery rate when classifying two gene (sgRNA) positive/negative references sets (tipically core-fitness-essential and non-essential genes), and it computes the Recall (or True Positive Rate) of genes in other user defined sets at the determined threshold. It produces a plot where the log fold changes are visualised alongside the rank positions of the genes included in the inputted sets and, their recall and the determined FDR threshold. 
}
\usage{
ccr.VisDepAndSig(FCsprofile,SIGNATURES,TITLE='',
                 pIs=NULL,nIs=NULL,
                 th=0.05,plotFCprofile=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FCsprofile}{
    A numerical vector containing gene average depletion log fold changes (or sgRNAs' depletion log fold changes) with names
corresponding to HGNC symbols (or sgRNAs' identifiers).
}
  \item{SIGNATURES}{
    A named list of vectors containing HGNC gene symbols. Two of these lists are used
    as classification template (respectively for positive and negative cases) to determine
    a log fold-change threshold providing a user defined classification false discovery rate.
}
  \item{TITLE}{
A string specifiying the title of the plot.
}  
  \item{pIs}{
The index position of the signature that contains the positive cases of the classification template.
}  
  \item{nIs}{
The index position of the signature that contains the negative cases of the classification template.
}  
  \item{th}{
A numerical value specifying the desired classification false discovery rate (this must be a real number between 0 and 1).
}  
  \item{plotFCprofile}{
A logic value specifying whether the log fold changes should be plotted.
}  
}
\value{
    A named numerical vector containing recall scores for all the inputted signatures at the computed false discovery rate threshold for log fold-changes. 
}
\author{
Francesco Iorio (iorio@gmail.com)
}
\seealso{
\code{\link{ccr.PrecisionRecallCurve}}
}
\examples{
## loading corrected sgRNAs log fold-changes and segment annotations
## for an example cell line (EPLC-272H)
data(EPLC.272HcorrectedFCs)

## loading reference sets of essential and non-essential genes
data(BAGEL_essential)
data(BAGEL_nonEssential)

## loading other sets of core fitness genes
data(EssGenes.ribosomalProteins)
data(EssGenes.DNA_REPLICATION_cons)
data(EssGenes.KEGG_rna_polymerase)
data(EssGenes.PROTEASOME_cons)
data(EssGenes.SPLICEOSOME_cons)

## storing the sgRNA log fold changes into a name vector
FCs<-EPLC.272HcorrectedFCs$corrected_logFCs$avgFC
names(FCs)<-rownames(EPLC.272HcorrectedFCs$corrected_logFCs)

## loading sgRNA library annotation
data(KY_Library_v1.0)

## computing gene average log fold changes
FCs<-ccr.geneMeanFCs(FCs,KY_Library_v1.0)


## Assembling a named list with all the considered gene sets
SIGNATURES<-list(Ribosomal_Proteins=EssGenes.ribosomalProteins,
                 DNA_Replication = EssGenes.DNA_REPLICATION_cons,
                 RNA_polymerase = EssGenes.KEGG_rna_polymerase,
                 Proteasome = EssGenes.PROTEASOME_cons,
                 Spliceosome = EssGenes.SPLICEOSOME_cons,
                 CFE=BAGEL_essential,
                 non_essential=BAGEL_nonEssential)

## Visualising log fold change profile with superimposed signatures specifying
## that the reference gene sets are in positions 6 and 7
Recall_scores<-ccr.VisDepAndSig(FCsprofile = FCs,
                            SIGNATURES = SIGNATURES,
                            TITLE = 'EPLC-272H',
                            pIs = 6,
                            nIs = 7)

Recall_scores
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Assessment and Visualisation}% use one of  RShowDoc("KEYWORDS")
