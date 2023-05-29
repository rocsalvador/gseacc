\name{GseaRcpp}
\alias{Gsea}
\title{Gsea}
\description{
Gsea creator function
}
\details{
Gsea creator function
}
\usage{
    # Create Gsea class to run gsea$runChunked()
    new(Gsea, sampleIds, geneIds, geneSets, nThreads)

    # Create Gsea class to run gsea$run()
    new(Gsea, expressionMatrix, geneSets, nThreads)
}
\arguments{
  \item{sampleIds}{sample ids as a character vector of the expression matrix}
  \item{geneIds}{gene ids as a character vector of the expression matrix}
  \item{expressionMatrix}{expression matrix as numeric matrix, it has no null rownames and colnames}
  \item{geneSets}{gene ids as a character vector of the expression matrix}
  \item{nThreads}{number of threads used to run GSEA, 0 if all CPU threads want to be used}
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
