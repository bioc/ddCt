\name{InputFrame}
\alias{InputFrame}
\alias{InputFrame,data.frame-method}
\title{Build an InputFrame from a ReaderClass or a data frame}
\description{
  Generally an InputFrame is built from a ReaderClass
  (e.g. \code{\linkS4class{InputReader}}), or a data.frame. See the example below
  for building an object from a valid data.frame.
}
\usage{
  InputFrame(object)
}
\arguments{
  \item{object}{A data.frame with three columns: Sample, Detector, and Ct}
}

\value{A object of class \code{InputFrame}}
\author{Jitao David Zhang \url{mailto: jitao_david.zhang@roche.com}}
\examples{
testDf <- data.frame(Sample=rep(paste("Sample", 1:3), each=2),
Detector=rep(paste("Gene", 1:2), 3),
Ct=30+rnorm(6))
testInputFrame <- InputFrame(testDf)
}
