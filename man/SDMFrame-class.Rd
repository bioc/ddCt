\name{SDMFrame-class}
\docType{class}
\alias{SDMFrame-class}
\alias{[,SDMFrame-method}
\alias{[[,SDMFrame-method}
\alias{$,SDMFrame-method}
\alias{detectorNames<-,SDMFrame,character-method}
\alias{detectorNames,SDMFrame-method}
\alias{fileNames,SDMFrame-method}
\alias{names,SDMFrame-method}
\alias{sampleNames<-,SDMFrame,character-method}
\alias{sampleNames,SDMFrame-method}
\alias{uniqueDetectorNames<-,SDMFrame,character,character-method}
\alias{uniqueDetectorNames<-,SDMFrame,missing,character-method}
\alias{uniqueDetectorNames,SDMFrame-method}
\alias{uniqueSampleNames<-,SDMFrame,character,character-method}
\alias{uniqueSampleNames<-,SDMFrame,missing,character-method}
\alias{uniqueSampleNames,SDMFrame-method}
\alias{coreData,SDMFrame-method}
\alias{coreData<-,SDMFrame,data.frame-method}
\alias{Ct,SDMFrame-method}
\alias{show,SDMFrame-method}

%% currently SDMFrame only
\alias{fileNames}
\alias{detectorNames<-}
\alias{detectorNames}
\alias{sampleNames<-}
\alias{sampleNames}
\alias{uniqueDetectorNames<-}
\alias{uniqueDetectorNames}
\alias{uniqueSampleNames<-}
\alias{uniqueSampleNames}
\alias{coreData}
\alias{coreData<-}
\alias{rightCensoring}
\alias{rightCensoring,SDMFrame,numeric-method}

\title{Working with SDM files}
\description{ The class SDMFrame provides core functionalities to read
  gene and sample information from SDM files and calculate them with a
  ddCt algorithm.

  The function \code{SDMFrame} reads the data given in the colums
  'Detector','Sample' and 'Ct' of the specified SDM output files and
  stores them as a data.frame. An additional column including the
  respective filename is added.
}
\section{Slots}{
  \describe{
    \item{\code{coreData}:}{Object of class \code{"data.frame"}: Holds
      all the required data extracted from the SDM file
    }
    \item{\code{files}:}{Object of class \code{"character"} contains the
      source SDM files
    }
  }
}
\section{Methods}{
  \describe{
    \item{[,[[,\$}{\code{signature(x = "SDMFrame")}: primitive
  accessors. Returns an object of \code{SDMFrame-class} with the subset data.
    }
    \item{names}{\code{signature(x = "SDMFrame")}: returns the column
      names in this SDM object
    }
    \item{ddCtExpression}{\code{signature(object = "SDMFrame")}: runs
      a ddCt algorithm with this SDM object and returns a object of
      class \code{\linkS4class{ddCtExpression}}
    }
    \item{fileNames}{\code{signature(object="SDMFrame")}: returns the
      source SDM file names.
    }
    \item{detectorNames}{\code{signature(object = "SDMFrame")}: returns
      the detector names in this SDM object
    }
    \item{detectorNames<-}{\code{signature(object = "SDMFrame", value =
	"character")}: replaces the detector names in this SDM object
    }
    \item{sampleNames}{\code{signature(object = "SDMFrame")}: returns
      the sample names in this SDM object
    }
    \item{sampleNames<-}{\code{signature(object = "SDMFrame", value =
	"character")}:  replaces the sample names in this SDM object
    }
    \item{uniqueDetectorNames}{\code{signature(object = "SDMFrame")}:
      returns a vector of unique detector names in this SDM object
    }
    \item{uniqueDetectorNames<-}{\code{signature(object = "SDMFrame",
	 target = "missing", value = "character")}:  replaces all detector
      names given by the 'names' attribute in 'value' with new detector
      names
    }
    \item{uniqueDetectorNames<-}{\code{signature(object = "SDMFrame",
	target = "character", value = "character")}:  replaces all detector
      names given by 'target' with new detector names
    }
    \item{uniqueSampleNames<-}{\code{signature(object = "SDMFrame",
	 target = "missing", value = "character")}: replaces all sample
      names given by the 'names' attribute in 'value' with new sample
      names
    }
    \item{uniqueSampleNames<-}{\code{signature(object = "SDMFrame",
	target = "character", value = "character")}: replaces all sample
      names given by 'target' with new sample names
    }
    \item{uniqueSampleNames}{\code{signature(object = "SDMFrame")}:
      returns a vector of unique sample names in this SDM object
    }
    \item{show}{\code{signature(object="SDMFrame")}: pretty print of the
      SDMFrame instance
    }
    \item{rightCensoring}{\code{signature(object="SDMFrame",
	threshold="numeric")}: Right censoring the Ct value,  which
      targets the data points above a certain value
      (\code{threshold}). High Ct values (higher than 40 or 45 by the rule
      of thumb) are often not accurate and may indicate too weak
      expression. The function performs the right censoring on the data and
      set the value above the threshold as \code{NA} (by default) or a
      given value. See the example.
    }
    \item{coreData}{\code{signature(object="SDMFrame")}: returns the
      data frame read from SDM file.
    }
    \item{coreData<-}{\code{signature(object="SDMFrame")}: replace the
      data frame read from SDM file.
    }
    \item{Ct}{\code{signature(object="SDMFrame")}: returns the Ct value
  of the SDM file.
    }
  }
}
\author{ Rudolf Biczok \email{mailto:r.biczok@dkfz.de}, Jitao David
  Zhang \email{mailto:j.zhang@dkfz.de}}
\seealso{
  \code{\link{SDMFrame}} function reads in data from SDM
  files. Data from SDM files is used to construct
  \code{\linkS4class{ddCtExpression}} objects to analyze differetial expression.
}
\examples{
## read a SDM file
sampdat <- SDMFrame(system.file("extdata", "Experiment1.txt",
                                package="ddCt"))
## you can also write 
## sampdat <- new("SDMFrame",system.file("extdata", "Experiment1.txt",
##                                       package="ddCt"))

## use the getter methods
sampleNames(sampdat)

## or the overloaded primitive accessors
sampdat[1:3,"Sample"]

## see all unique samples
uniqueSampleNames(sampdat)

## replace all sample names 'Sample1' and 'Sample2' in sampdat
## with 'NewSample1' and 'NewSample2'
uniqueSampleNames(sampdat,c("Sample1","Sample2")) <- c("NewSample1","NewSample2")
uniqueSampleNames(sampdat)

## or use this syntax to replace the gene names
uniqueDetectorNames(sampdat) <- c(Gene1="NewGene1", Gene2="NewGene2")
uniqueDetectorNames(sampdat)

## right censoring the data
rightCensoring(sampdat, 35)
rightCensoring(sampdat, 35, 35)
}
\keyword{classes}