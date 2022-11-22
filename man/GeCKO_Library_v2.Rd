\name{GeCKO_Library_v2}
\alias{GeCKO_Library_v2}
\docType{data}
\title{
Genome-wide annotation for the GeCKO (v2) sgRNA library
}
\description{
A data frame with a named row for each sgRNA of the GeCKO sgRNA library [1] including annotations such as targeted genes, and genomic coordinates.
}
\usage{data(GeCKO_Library_v2)}
\format{
    A a row named data frame with 121327 observations of the following variables (among others)
  \describe{
   \item{\code{CODE}}{alphanumerical identifier of the sgRNAs;}
    \item{\code{GENES}}{targeted gene;}
    \item{\code{STARTpos}}{starting genomic coordinate of the targeted genomic region (numeric);}
    \item{\code{STRAND}}{targeted DNA strand ('+' or '-')}
    \item{\code{EXONE}}{exone of the targeted genomic region (exone number);}
    \item{\code{CHRM}}{chromosome of where the targeted region resides (string)}
    \item{\code{ENDpos}}{ending genomic coordinate of the targeted genomic region (numeric).}
    \item{\code{seq}}{nucelotidic sequence of the sgRNAs without the PAM. (string).}
  }
}
\details{
GeCKO v2 library was developed with the aim of targeting all genes with a uniform number of sgRNAs, and included 6 sgRNAs per gene distributed over 3-4 constitutively expressed exons. Minimization of off-target effects was based on a specificity analysis. In addition the library included a number of sgRNAs targeting microRNAs (miRNAs) and 2,000 non targeting sgRNAs, for a total number of 123411 sgRNAs. 

Genomic coordinates of the sgRNAs (required by CRISPRcleanR) of the GeCKO v2 library were not available on the annotation file avaialable on AddGene [2], although some partial mappings are provided.

We generated the locations of these mapping positions on the reference genome using the sequence content of the sgRNAs available in the library annotation, using the latest human reference genome (GRCh38), using multiple tools, as detailed in the following steps:

\describe{
    \item{Step 1}{The sgRNAs were mapped onto the human reference sequence using the bwa short read mapper. Only the reads that were mapped to the reference genome uniquely were selected and their positions of mappings (start/and end positions) were superimposed to those of the intended targeted gene in Ensembl gene annotation v100. From these mappings all the sgRNAs that were mapped to the correct corresponding target genes were identified and retained. Although bwa is an efficient mapper, due to multiple mapping locations and some small insertions and deletions, some sgRNAs were not mapped to the reference sequence.}
    
    \item{Step 2}{All the sgRNAs that were mapped to the reference genome in multiple locations were selected and overlapped with the intended targeted gene locations. The sgRNAs mapped onto at least one gene-matching location were selected and retained.}
    
    \item{Step 3}{All the sgRNAs that were not mapped to their intended/declared target gene were selected and the intended gene symbol/name checked for alternative/more-recent gene symbols/names. All possible alternative gene names were identified and checked for overlap. After correction some of the mappings were corrected and the corresponding sgRNA retained.}
    
    \item{Step 4}{All the remaining sgRNAs (missing or not mapped) were selected and mapped to the reference genome using the blast tool. Here the mapping is slower but more accurate. The results of the blast psl files including all possible mappings of sgRNAs were parsed. The positions were similarly compared to the reference gene annotations and corrected for most recent gene names/symbols. The sgRNAs correctly mapped were retained.}
    
    \item{Step 5}{All the remaining sgRNAs were compared against miRBase for non-coding RNAS and for the consistency of the naming of these miRNAs. The matching sgRNAs were identified and retained.}
    
    \item{Step 6}{The remaining sgRNAs, matching many locations in the human reference genome or with an intendend target name different from that in the annotation file, were mapped to their targeted region using the Waterman-Smith local alignment manually. All the remaining sgRNAs were manually curated retained.}
    
    \item{Step 7}{Some of the sgRNAs were not added to the final annotation data object. The main reason for this is that these genes were removed from the primary human reference in the GRCh38 version. Also, some miRNAs are retracted as well as some genes. Finally some sgRNAs did not map to the gene that they are intended to target. 

These removed sgRNAs were declared to target:

GTF2H2D,
LILRA3,
LOC391322,
LOC653486,
PRAMEF16,
SMCR9,
hsa-mir-1273a,
hsa-mir-1273d,
hsa-mir-1273g,
hsa-mir-1302-5,
hsa-mir-3118-5,
hsa-mir-3118-6,
hsa-mir-320d-1,
hsa-mir-3669,
hsa-mir-3673,
hsa-mir-3910-1,
hsa-mir-3910-2,
hsa-mir-4419a,
hsa-mir-4419b,
hsa-mir-4459,
hsa-mir-4472-2,
hsa-mir-5096,
hsa-mir-548aa-2,
hsa-mir-548d-2,
hsa-mir-6087,
GTF2H2D,
LILRA3,
LOC391322,
LOC653486,
PRAMEF16,
SMCR9
}
}

The final list of retained sgRNAs included 121320  (out of 123411). Note that 2000 of the excluded sgRNAs were Non-targeting. 

}
\source{
Source of the Library: AddGene, https://www.addgene.org/pooled-library/zhang-human-gecko-v2
Source of the annotation file used for the sgRNA remapping: 
https://sourceforge.net/projects/mageck/files/libraries/Human_GeCKOv2_Library_combine.csv.zip/download

Sources for the tools:

\describe{
    \item{blat}{
    Standalone BLAT v. 36x5 fast sequence search command line tool.
    The executable is downloadable from http://hgdownload.soe.ucsc.edu/admin/exe/
    }
    \item{blast}{
    blastn: 2.6.0+
    Package: blast 2.6.0, build Jan 15 2017 17:12:27
    https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch
    }
    \item{miRbase}{
    http://www.mirbase.org/
    From available databases, we used has.gff3 which is from human reference sequence hg38
    }
    \item{Smith-Waterman}{
    https://www.ebi.ac.uk/Tools/psa/emboss_water/
    The web interface is used for the local alignment as well as the command line via REST API.
    }
}
}
\references{
[1] Sanjana NE, Shalem O, Zhang F. Improved vectors and genome-wide libraries for CRISPR screening. Nat Methods. 2014;11(8):783-784. doi:10.1038/nmeth.3047

[2] Aguirre AJ, Meyers RM, Weir BA, et al. Genomic Copy Number Dictates a Gene-Independent Cell Response to CRISPR/Cas9 Targeting. Cancer Discov. 2016;6(8):914-929. doi:10.1158/2159-8290.CD-16-0154
}
\examples{\dontrun{
## Loading sgRNA GeCKO library annotation file
data(GeCKO_Library_v2)
## Visualising first entries
head(GeCKO_Library_v2)

## Deriving the path of an example count file
## from screening the HT-29 cell line with the GeCKO v2
## library [2]
fn<-paste(system.file('extdata', package = 'CRISPRcleanR'),
          '/HT-29-GeCKOv2_counts.tsv',sep='')

expName<-'HT29-GeCKOv2'

## Loading, median-normalizing and computing fold-changes
normANDfcs<-
    ccr.NormfoldChanges(filename = fn,
                        display = TRUE,
                        min_reads = 30,
                        EXPname = expName,
                        libraryAnnotation = GeCKO_Library_v2)

## Genome-sorting the fold changes
gwSortedFCs<-
    ccr.logFCs2chromPos(foldchanges = normANDfcs$logFCs,
                        libraryAnnotation = GeCKO_Library_v2)

## Identifying and correcting biased sgRNAs' fold changes
correctedFCs_and_segments<-
    ccr.GWclean(gwSortedFCs = gwSortedFCs,
                display=TRUE,
                label=expName)
}}
\keyword{datasets}
