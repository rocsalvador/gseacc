\name{Gsea}
\alias{Gsea}
\title{Gsea}
\description{
Gsea creator function.

Class methods:

- \code{?run}

- \code{?runChunked}

- \code{?filterResults}

- \code{?normalizeExprMatrix}
}
\usage{
    # Create Gsea class to run gsea$runChunked()
    new(Gsea, sampleIds, geneIds, geneSets, nThreads)

    # Create Gsea class to run gsea$run()
    new(Gsea, expressionMatrix, geneSets, nThreads)
}
\arguments{
  \item{sampleIds}{Sample ids as a character vector of the expression matrix}
  \item{geneIds}{Gene ids as a character vector of the expression matrix}
  \item{expressionMatrix}{Expression matrix as numeric matrix, it has no null rownames and colnames}
  \item{geneSets}{Gene sets as list}
  \item{nThreads}{Number of threads used to run GSEA, 0 if all CPU threads want to be used}
}
\examples{

expressionMatrix <- readCsv("expressionMatrix.csv")
geneSets <- readGeneSets("geneSets.csv")
geneIds <- colnames(expressionMatrix)
sampleIds <- rownames(expressionMatrix)

gsea <- new(Gsea, sampleIds, geneIds, geneSets, 0)

gsea$runChunked(expressionMatrix[1:100])

gsea <- new(Gsea, expressionMatrix, geneSets, 4)

gsea$run("results.csv", 10)
}
