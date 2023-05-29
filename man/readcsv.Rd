\name{readCsv}
\alias{readCsv}
\title{readCsv}
\description{
Read a csv file as a numeric matrix
}
\usage{
    readCsv(fileName, sep = ",", hasRowNames = TRUE, hasColNames = TRUE)
}
\arguments{
    \item{fileName}{Name of the csv file}
    \item{sep}{Separator of csv elements}
    \item{hasRowNames}{True if it has row names, false otherwise}
    \item{hasColNames}{True if it has column names, false otherwise}
}
\examples{

expressionMatrix <- readCsv("expressionMatrix.csv")

expressionMatrix <- readCsv("expressionMatrix.csv", "\t", FALSE, TRUE)
}
