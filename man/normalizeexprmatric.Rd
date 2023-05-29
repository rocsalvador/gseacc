\name{normalizeExprMatrix}
\alias{normalizeExprMatrix}
\title{normalizeExprMatrix}
\description{
Normalize the expression matrix using rpm and then mean center
}
\usage{
    gsea$normalizeExprMatrix()
}
\examples{

expressionMatrix <- readCsv("expressionMatrix.csv")
geneSets <- readGeneSets("geneSets.csv")
geneIds <- colnames(expressionMatrix)
sampleIds <- rownames(expressionMatrix)

gsea <- new(Gsea, sampleIds, geneIds, geneSets, 0)

gsea$normalizeExprMatrix()

gsea$run("results.csv", 10)
}
