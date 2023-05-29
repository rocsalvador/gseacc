\name{runChunked}
\alias{runChunked}
\title{runChunked}
\description{
Run Gsea for the given chunk of the expression matrix
}
\usage{
gsea$runChunked(expressionMatrix)
}
\arguments{
  \item{expressionMatrix}{Numeric matrix containing a chunk of the expression matrix with samples in the rows and genes in the columns}
}
\examples{

expressionMatrix <- readCsv("expressionMatrix.csv")
geneSets <- readGeneSets("geneSets.csv")
geneIds <- colnames(expressionMatrix)
sampleIds <- rownames(expressionMatrix)

gsea <- new(Gsea, sampleIds, geneIds, geneSets, 0)

gsea$runChunked(expressionMatrix[1:100])
}
