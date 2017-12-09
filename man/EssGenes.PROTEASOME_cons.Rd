\name{EssGenes.PROTEASOME_cons}
\alias{EssGenes.PROTEASOME_cons}
\docType{data}
\title{
 Core Fitness essential proteasome genes
}
\description{
List of core fitness essential proteasome genes assembled by merging together multilpe DNA replication signatures from MSigDB [1] as detailed in [2].
}
\usage{data("EssGenes.PROTEASOME_cons")}
\format{
 A vector of strings containing HGNC symbols.
}
\references{
[1] Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., et al. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proceedings of the National Academy of Sciences of the United States of America, 102(43), 15545-15550. http://doi.org/10.1073/pnas.0506580102 \cr

[2] Iorio, F., Behan, F. M., Goncalves, E., Beaver, C., Ansari, R., Pooley, R., et al. (n.d.). Unsupervised correction of gene-independent cell responses to CRISPR-Cas9 targeting. \cr
http://doi.org/10.1101/228189
}
\examples{
data(EssGenes.PROTEASOME_cons)
head(EssGenes.PROTEASOME_cons)
}
\keyword{datasets}
