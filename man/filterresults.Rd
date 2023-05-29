\name{filterResults}
\alias{filterResults}
\title{filterResults}
\description{
Filters the the most variable gene sets from gsea$runChunked() results.
}
\usage{
    gsea$filterResults(nFilteredGeneSets)
}
\arguments{
  \item{nFilteredGeneSets}{number of gene sets in the filtered results file}
}
\examples{

expressionMatrix <- readCsv("expressionMatrix.csv")
geneSets <- readGeneSets("geneSets.csv")
geneIds <- colnames(expressionMatrix)
sampleIds <- rownames(expressionMatrix)

gsea <- new(Gsea, sampleIds, geneIds, geneSets, 0)

gsea$runChunked(expressionMatrix[1:100])
gsea$runChunked(expressionMatrix[101:200])

gsea$filterResults(50)
}
