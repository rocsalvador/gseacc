\name{run}
\alias{run}
\title{run}
\description{
Run Gsea for the given expression matrix and gene sets in the Gsea creator function
}
\usage{
    gsea$run(outFileName, ioutput)
}
\arguments{
    \item{outFileName}{Name of the results file}
    \item{ioutput}{Number of gene sets between status output}
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
