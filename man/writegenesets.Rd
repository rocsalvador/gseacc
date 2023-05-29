\name{writeGeneSets}
\alias{writeGeneSets}
\title{writeGeneSets}
\description{
Write gene sets to a gene sets file
}
\usage{
    writeGeneSets(geneSets, fileName)
}
\arguments{
    \item{geneSets}{List of gene sets}
    \item{fileName}{Name of the gene sets file}
}
\examples{
library("AnnotationDbi")
library("org.Hs.eg.db")

xx <- as.list(org.Hs.egGO2ALLEGS)

geneSets <- lapply(xx, function(x) {
    mapIds(org.Hs.eg.db, keys = x, keytype = "ENTREZID", column = "ENSEMBL")
})

writeGeneSets(geneSets, "geneSets.csv")

geneSets <- readGeneSets("geneSets.csv")
}
