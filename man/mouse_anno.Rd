% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{mouse_anno}
\alias{mouse_anno}
\title{a example output of \code{read_annotation}}
\format{
A \code{GRanges} object.

The GRanges object contains three additional metadata columns:
     \itemize{ \item{ \code{id}: EnsemblID}  
     \item{ \code{center}: the center of the region,(strat+end)/2} 
     \item{ \code{annotation}:the annotation of the region, it can be promoter,downstream,gene,UTR...} }
}
\usage{
mouse_anno
}
\value{
annotation data
}
\description{
a samll subset mouse genomic annotation of chr15
}
\seealso{
\code{\link{read_annotation}}
}
\keyword{datasets}
