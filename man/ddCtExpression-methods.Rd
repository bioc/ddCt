\name{ddCtExpression-methods}
\docType{methods}
\alias{ddCtExpression}
\alias{ddCtExpression-methods}
\alias{ddCtExpression,SDMFrame-method}
\alias{ddCtExpression,InputFrame-method}
\title{Apply the ddCt algorithm for a given data set}
\description{Apply the ddCt algorithm for a given data set}
\section{usage}{
    ddCtExpression(object,
                   warningStream = "warning.output.txt",
                   algorithm="ddCt"
                   calibrationSample,
                   housekeepingGenes,
                   type="mean",
                   sampleInformation=NULL,
                   toZero=TRUE,
                   efficiencies = NULL,
                   efficiencies.error = NULL)
 }
\arguments{
  \item{object}{SDMFrame Data object which holds a data set containing columns with the
    following names: 'Ct','Sample','Detector','Platename'. The column
    'Ct' must contain numeric values.}
  \item{algorithm}{character. Name of the calibration samples.}
  \item{warningStream}{character of length 1. The name of the file the
    warnings should be stored in.}
  \item{calibrationSample}{character. Name of the calibration samples.}
  \item{housekeepingGenes}{character. Name of the housekeeping genes.}
  \item{type}{character of length 1. `mean` or `median`- which method should be used for the
    aggregation of the repicates}
  \item{sampleInformation}{if specified it must be an object of class \code{phenoData} with a column named 'Sample'.}
  \item{toZero}{boolean - if there is only one replication should the error be
    treated as zero ? (only if 'type' is mean)}
  \item{efficiencies}{n.V.}
  \item{efficiencies.error}{n.V.}
}
\value{A an object of class
  \code{\linkS4class{ddCtExpression}}.
}
\section{Methods}{
  \describe{
    \item{object = "InputFrame"}{ An object of \code{\linkS4class{InputFrame}},
      constructed with the method \code{\link{InputFrame}}
    }
  }
}
\references{
  Analysis of relative gene expression data using real-time quantitative
  PCR and the 2(-Delta -Delta C(T)) Method. KJ Livak and TD Schmittgen,
  Methods, Vol. 25, No. 4. (December 2001), pp. 402-408
}
\seealso{
  \code{\linkS4class{InputFrame}}: reader for SDM files
  \code{\linkS4class{ddCtExpression}}:  representation for ddCt
  calculated expressions
}

\author{ Rudolf Biczok \email{mailto:r.biczok@dkfz.de}}

\examples{
## read a SDM file
sampdat <- SDMFrame(system.file("extdata", "Experiment1.txt", package="ddCt"))

## call ddCtExpression method from class SDMFrame
## to get a ddCt calculated expression
result <- ddCtExpression(sampdat,
                         calibrationSample="Sample1",
                         housekeepingGenes=c("Gene1","Gene2"))
result
}
