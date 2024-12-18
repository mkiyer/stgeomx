% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalize.R
\name{normalize}
\alias{normalize}
\title{normalization internal function}
\usage{
normalize(ds, method = c("qn", "rle", "cpm", "q3", "tmm", "none"))
}
\arguments{
\item{ds}{dataset}

\item{method}{string normalization method}
}
\value{
matrix of normalized counts
}
\description{
scale/quantile normalization methods
}
\examples{
data(example_ds, package = "stgeomx")
example_ds <- preprocess(example_ds)
example_ds <- merge_probes(example_ds)
normalize(example_ds, "qn")
normalize(example_ds, "rle")
normalize(example_ds, "cpm")
normalize(example_ds, "q3")
normalize(example_ds, "tmm")
normalize(example_ds, "none")

}
