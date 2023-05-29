\name{readGeneSets}
\alias{readGeneSets}
\title{readGeneSets}
\description{
Read gene sets from a gene sets file
}
\usage{
    readGeneSets(fileName)
}
\arguments{
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
