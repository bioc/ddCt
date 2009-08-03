\name{ddCtAbsolute}
\alias{ddCtAbsolute}
\title{absolute quantification for Taqman data}
\description{
  absolute quantification for Taqman data
}
\usage{
ddCtAbsolute(raw.table, addData, type = "mean", ADD = -30.234, DIV = -1.6268, sampleInformation = NULL, toZero = FALSE, filename = "warning.output.txt")
}

\arguments{
  \item{raw.table}{data frame. It must contain columns with the
    following names:'Ct','Sample','Detector','Platename'. The column
    'Ct' must contain numeric values. }
  \item{addData}{}
  \item{type}{character of length 1. `mean` or `median`- which method should be used for the
    aggregation of the repicates }
  \item{ADD}{  }
  \item{DIV}{  }
  \item{sampleInformation}{ if specified it must be an object of class
    \code{phenoData} with a column named 'Sample'. }
  \item{toZero}{boolean - if there is only one replication should the error be
    treated as zero ? (only if 'type' is mean) }
  \item{filename}{ character of length 1. The name of the file the
    warnings should be stored in. }
}

\value{A an object of class
  \code{\link[Bioconductor:eSet-class]{eSet}}. The assayData
  has the following components: exprs, error, Ct, Ct.error, Difference, number\_NA, number, Plate.}
\references{ ~put references to the literature/web site here ~ }
\author{Markus Ruschhaupt \url{mailto:m.ruschhaupt@dkfz.de} }

 
