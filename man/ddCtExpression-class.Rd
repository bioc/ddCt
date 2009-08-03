\name{ddCtExpression-class}
\docType{class}
\alias{ddCtExpression-class}
\alias{Ct,ddCtExpression-method}
\alias{CtErr,ddCtExpression-method}
\alias{dCt,ddCtExpression-method}
\alias{dCtErr,ddCtExpression-method}
\alias{ddCt,ddCtExpression-method}
\alias{ddCtErr,ddCtExpression-method}
\alias{elist,ddCtExpression-method}
\alias{elistWrite,ddCtExpression,character-method}
\alias{level,ddCtExpression-method}
\alias{levelErr,ddCtExpression-method}
\alias{numberCt,ddCtExpression-method}
\alias{numberNA,ddCtExpression-method}
%% currently ddCtExpression only
\alias{Ct}
\alias{CtErr}
\alias{dCt}
\alias{dCtErr}
\alias{ddCt}
\alias{ddCtErr}
\alias{level}
\alias{levelErr}
\alias{numberCt}
\alias{numberNA}

\title{ddCt Expression}
\description{This class is a subclass of \code{ExpressionSet} and
  represents objects which are produced
  by the ddCt algorithm in the \code{ddCtExpression} method
}
\section{Extends}{
  Class \code{ExpressionSet}, directly.
  Class \code{eSet}, by class "ExpressionSet", distance 2.
  Class \code{VersionedBiobase}, by class "ExpressionSet", distance 3.
  Class \code{Versioned}, by class "ExpressionSet", distance 4.
}
\section{Methods}{
  \describe{
    \item{Ct}{\code{signature(object = "ddCtExpression")}: returns the
      Ct value of this ddCtExpressionobject
    }
    \item{CtErr}{\code{signature(object = "ddCtExpression")}: returns
      the error number of the Ct value of this ddCtExpressionobject
    }
    \item{dCt}{\code{signature(object = "ddCtExpression")}: returns the
      dCt value of this ddCtExpressionobject
    }
    \item{dCtErr}{\code{signature(object = "ddCtExpression")}: returns
      the error number of the dCt value of this ddCtExpressionobject
    }
    \item{ddCt}{\code{signature(object = "ddCtExpression")}:returns the
      ddCt value of this ddCtExpressionobject
    }
    \item{ddCtErr}{\code{signature(object = "ddCtExpression")}: returns
      the error number of the ddCt value of this ddCtExpressionobject
    }
    \item{level}{\code{signature(object = "ddCtExpression")}: returns
      the levels in this ddCtExpressionobject
    }
    \item{levelErr}{\code{signature(object = "ddCtExpression")}: returns
      the error number of the levens in this ddCtExpressionobject
    }
    \item{numberCt}{\code{signature(object = "ddCtExpression")}: returns
      the Ct number of this ddCtExpressionobject
    }
    \item{numberNA}{\code{signature(object = "ddCtExpression")}: returns
      the NA number of this ddCtExpressionobject
    }
    \item{elist}{\code{signature(object = "ddCtExpression")}: returns a
      data frame which represents this expression object}
    \item{elistWrite}{\code{signature(object = "ddCtExpression", file =
	"character")}: {writes ddCtExpression object into a file}
    }
  }
}
\author{ Rudolf Biczok \email{mailto:r.biczok@dkfz.de}}

\seealso{
  \code{\linkS4class{SDMFrame}}: reader for SDM files
  \code{\link{elist}}, \code{\link{elistWrite}}: utility functions for
  ddCtExpression objects
  \code{\link{ddCtExpression}}:  the method which invokes the ddCt algorithm
}
\examples{
## read a SDM file
sampdat <- SDMFrame(system.file("extdata", "Experiment1.txt",
                                package="ddCt"))

## call ddCtExpression method to get a ddCt calculated expression
result <- ddCtExpression(sampdat,
                         calibrationSample="Sample1",
                         housekeepingGenes=c("Gene1","Gene2"))

## use getter methods
ddCt(result)
ddCtErr(result)
}
\keyword{classes}