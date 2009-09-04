\name{errBarchartParameter-class}
\Rdversion{1.1}
\docType{class}
\alias{errBarchartParameter-class}
\alias{exprsUndeterminedLabel,errBarchartParameter-method}
\alias{show,errBarchartParameter-method}

\title{Class "errBarchartParameter"}
\description{ Parameter object for errBarchart }
\section{Objects from the Class}{
Objects can be created by calls of the form
\code{new("errBarchartParameter", ...)}. So far the object is only
internally used, but in the near future it will be exported.
}
\section{Slots}{
  \describe{
    \item{\code{exprsUndeterminedLabel}:}{Object of class
      \code{"character"}, specifying the text label when the expression
      level is \sQuote{Undetermined}}
  }
}
\section{Methods}{
  \describe{
    \item{exprsUndeterminedLabel}{\code{signature(object =
	"errBarchartParameter")}:  getting the text label when the
      expression level is \sQuote{Undetermined}}
    \item{show}{\code{signature(object = "errBarchartParameter")}:
      print method}
	 }
}
\author{ Jitao David Zhang <j.zhang@dkfz.de>}
\note{ So far it is only internallly used }
\examples{
## Internally used
## param <- new("errBarchartParameter")
## exprsUndeterminedLabel(param)
}
\keyword{classes}
