\name{multiWGCNAdata-package}

\alias{multiWGCNAdata-package}
\alias{multiWGCNAdata}

\title{multiWGCNAdata available as R objects}

\description{
  Stores expression profiling data from experiments compatible with
  the multiWGCNA R package. This includes human postmortem microarray data
  from patients and controls (GSE28521), astrocyte Ribotag RNA-seq data from
  EAE and wildtype mice (GSE100329), and mouse RNA-seq data from
  tau pathology (rTg4510) and wildtype control mice (GSE125957).
  These data have been parsed into \code{list} objects and
  are available in the \code{ExperimentHub}.
}

\details{
  See the vignette for examples of using these data in analyses.

  browseVignettes("multiWGCNAdata")

  Details of how these data were created are in the scripts/ directory
  of the source package.
}

\value{
NA
}

\examples{
library(ExperimentHub)
eh = ExperimentHub()
eh_query = query(eh, c("multiWGCNAdata"))
astrocyte_se = eh_query[["EH8223"]]
astrocyte_se
}

\keyword{utilities}
