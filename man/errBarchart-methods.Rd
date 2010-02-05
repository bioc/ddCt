\name{errBarchart-methods}
\docType{methods}
\alias{errBarchart}
\alias{errBarchart-methods}
\alias{errBarchart,ddCtExpression-method}
\title{ Draw barchart of relative expression level with error-bars}
\description{
  Draw barchart (with error-bars) of relative expression level
  represented in \code{\linkS4class{ddCtExpression}} object. The
  barchart is implemented as grid plot by \code{lattice} package, where
  each panel represents one sample and the relative expression values of detectors (as
  well as their standard errors) are depicted as bars.

  Detectors which are not determined are marked by grey ND.

  Two types of figures are supported: either condition on samples
  (by="Sample") or on detectors (by="Detector").
}
\section{Methods}{
\describe{
\item{object = "ddCtExpression"}{ An object of \code{ddCtExpression},
  constructed with the method \code{\link{ddCtExpression}}}
}}
\keyword{methods}

