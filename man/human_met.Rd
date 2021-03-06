% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{human_met}
\alias{human_met}
\title{a example output of \code{read_methylation_report}}
\format{
A \code{GRanges} object.

  The GRanges object contains three additional metadata columns:
  \itemize{ \item{ \code{methylated_reads}: the number if methylated reads for each C sites.}  
  \item{ \code{total_reads}: Total number of reads for each C sites.} 
  \item{ \code{type}:the C-context type: CG,CHG, or CHH} 
  }
}
\usage{
human_met
}
\value{
methylation data
}
\description{
a samll subset methylation data of human chr22
}
\seealso{
\code{\link{read_methylation_report}}
}
\keyword{datasets}
