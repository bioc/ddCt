\name{elistWrite-methods}
\docType{methods}
\alias{elist}
\alias{summary}
\alias{elistWrite}
\title{Write ddCtExpression object into data frame or files}
\description{\code{ddCtExpression} object contains a list of matrices as
the results of \code{\link{ddCt}} method. \code{elist} combines these
lists into one data frame, and \code{elistWrite} writes the data frame
into file.

\code{summary} is a wrapper for the \code{elist} method
}
\usage{
elist(object,...)
summary(object,...)
elistWrite(object,file,...)
}
\arguments{
  \item{object}{an \code{\link[ddCt:ddCtExpression-class]{ExpressionSet}} object.}
  \item{file}{output file.}
  \item{...}{additional arguments passed to write.table.}
}

\value{A data frame or output file.}

\details{
  \code{elist} is a wrapper to \code{as(object, "data.frame")} function.
}

\author{Jitao David Zhang \email{jitao_david.zhang@roche.com}}

\examples{
## read a SDM file
sampdat <- SDMFrame(system.file("extdata", "Experiment1.txt", package="ddCt"))

## call ddCtExpression method from class SDMFrame
## to get a ddCt calculated expression
result <- ddCtExpression(sampdat,
                         calibrationSample="Sample1",
                         housekeepingGenes=c("Gene1","Gene2"))

## call elist
elistResult <- elist(result)
elistResult
}
