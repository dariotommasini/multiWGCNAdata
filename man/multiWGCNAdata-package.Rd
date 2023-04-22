\name{multiWGCNAdata-package}

\alias{multiWGCNAdata-package}
\alias{multiWGCNAdata}

\title{multiWGCNAdata available as R list objects}

\description{
  Stores expression profiling data from experiments compatible with
  the multiWGCNA R package. This includes human postmortem microarray data
  from patients and controls (GSE28521), astrocyte Ribotag RNA-seq data from
  EAE and wildtype mice (GSE100329), and mouse RNA-seq data from
  tau pathology (rTg4510) and wildtype control mice (GSE125957).
  These data have been parsed into \code{list} objects and
  are available in \code{ExperimentHub}.
}

\details{
  See the vignette for examples of using these data in differential
  gene expression analysis.

  browseVignettes("multiWGCNAdata")

  Details of how these data were creates are in the scripts/ directory
  of the source package.
}

\examples{
library(ExperimentHub)
hub <- ExperimentHub()
x <- query(hub, c("multiWGCNAdata", "autism"))
x
y <- query(hub, c("multiWGCNAdata", "astrocyte"))
y
z <- query(hub, c("multiWGCNAdata", "tau"))
z
\dontrun{
## download resource
autism_data = x[[1]]
astrocyte_data = y[[1]]
tau_data = z[[1]]
}}

\keyword{utilities}
