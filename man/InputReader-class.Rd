\name{InputReader-class}
\Rdversion{1.1}
\docType{class}
\alias{InputReader-class}
\alias{SDMReader-class}
\alias{TSVReader-class}
\alias{QuantStudioReader-class}
\alias{ColMap}
\alias{TSVFrame}


\title{Class "InputReader"}
\description{
  Abstract factory for data input
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{files}:}{Input files}
    \item{\code{colmap}:}{Column mapping}
  }
}
\author{
  Rudolf Biczok and Jitao David Zhang
}
\examples{
showClass("InputReader")
}
\keyword{classes}
