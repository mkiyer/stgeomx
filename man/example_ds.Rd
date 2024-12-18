% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{example_ds}
\alias{example_ds}
\title{Example GeoMx DSP Dataset}
\format{
## `example_ds`
A list of length 4:
\describe{
  \item{samples}{tibble of 2 rows and 36 columns, each row representing a region of interest}
  \item{meta}{tibble of 2 rows and 3 columns}
   \describe{
     \item{probe}{list of probes}
     \item{gene}{list of genes that correspond to the probes}
     \item{bg}{list correlating to background noise}
   }
  \item{counts}{tibble of 2 rows and 3 columns where each row contains the probes counts for each region of interest}
  \item{x}{tibble of 2 rows and 3 columns where each row contains the probe counts for each region of interest}
  ...
}
}
\source{
<https://nanostring.com/products/geomx-digital-spatial-profiler/spatial-organ-atlas/>
}
\usage{
example_ds
}
\description{
A small subset of data from the Nanostring Spatial Organ Atlas Dataset
}
\keyword{datasets}
