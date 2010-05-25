\name{InputFrame-class}
\docType{class}
\alias{SDMFrame-class}
\alias{InputFrame-class}
\alias{ColMap-class}
\alias{[,InputFrame-method}
\alias{[[,InputFrame-method}
\alias{$,InputFrame-method}
\alias{detectorNames<-,InputFrame,character-method}
\alias{detectorNames,InputFrame-method}
\alias{fileNames,InputFrame-method}
\alias{names,InputFrame-method}
\alias{sampleNames<-,InputFrame,character-method}
\alias{sampleNames,InputFrame-method}
\alias{uniqueDetectorNames<-,InputFrame,character,character-method}
\alias{uniqueDetectorNames<-,InputFrame,missing,character-method}
\alias{uniqueDetectorNames,InputFrame-method}
\alias{uniqueSampleNames<-,InputFrame,character,character-method}
\alias{uniqueSampleNames<-,InputFrame,missing,character-method}
\alias{uniqueSampleNames,InputFrame-method}
\alias{coreData,InputFrame-method}
\alias{coreData<-,InputFrame,data.frame-method}
\alias{Ct,InputFrame-method}
\alias{Ct<-}
\alias{Ct<-,InputFrame,numeric-method}
\alias{show,InputFrame-method}


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

%% replace and remove
\alias{removeDetector}
\alias{removeDetector-methods}
\alias{removeDetector,InputFrame,character-method}
\alias{removeSample}
\alias{removeSample-methods}
\alias{removeSample,InputFrame,character-method}
\alias{replaceDetector}
\alias{replaceDetector-methods}
\alias{replaceDetector,InputFrame,character,character-method}
\alias{replaceSample}
\alias{replaceSample-methods}
\alias{replaceSample,InputFrame,character,character-method}

\title{Working with SDM files}
\description{ The class InputFrame provides core functionalities to read
  gene and sample information from SDM files and calculate them with a
  ddCt algorithm.

  The function \code{InputFrame} reads the data given in the colums
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
    \item{[,[[,\$}{\code{signature(x = "InputFrame")}: primitive
  accessors. Returns an object of \code{InputFrame-class} with the subset data.
    }
    \item{names}{\code{signature(x = "InputFrame")}: returns the column
      names in this SDM object
    }
    \item{ddCtExpression}{\code{signature(object = "InputFrame")}: runs
      a ddCt algorithm with this SDM object and returns a object of
      class \code{\linkS4class{ddCtExpression}}
    }
    \item{fileNames}{\code{signature(object="InputFrame")}: returns the
      source SDM file names.
    }
    \item{detectorNames}{\code{signature(object = "InputFrame")}: returns
      the detector names in this SDM object
    }
    \item{detectorNames<-}{\code{signature(object = "InputFrame", value =
	"character")}: replaces the detector names in this SDM object
    }
    \item{sampleNames}{\code{signature(object = "InputFrame")}: returns
      the sample names in this SDM object
    }
    \item{sampleNames<-}{\code{signature(object = "InputFrame", value =
	"character")}:  replaces the sample names in this SDM object
    }
    \item{uniqueDetectorNames}{\code{signature(object = "InputFrame")}:
      returns a vector of unique detector names in this SDM object
    }
    \item{uniqueDetectorNames<-}{\code{signature(object = "InputFrame",
	 target = "missing", value = "character")}:  replaces all detector
      names given by the 'names' attribute in 'value' with new detector
      names
    }
    \item{uniqueDetectorNames<-}{\code{signature(object = "InputFrame",
	target = "character", value = "character")}:  replaces all detector
      names given by 'target' with new detector names
    }
    \item{uniqueSampleNames<-}{\code{signature(object = "InputFrame",
	 target = "missing", value = "character")}: replaces all sample
      names given by the 'names' attribute in 'value' with new sample
      names
    }
    \item{uniqueSampleNames<-}{\code{signature(object = "InputFrame",
	target = "character", value = "character")}: replaces all sample
      names given by 'target' with new sample names
    }
    \item{uniqueSampleNames}{\code{signature(object = "InputFrame")}:
      returns a vector of unique sample names in this SDM object
    }
    \item{removeSample}{\code{signature(object = "InputFrame", sample="character")}:
      removes the sample(s) specified from the InputFrame object
    }
    \item{replaceDetector}{\code{signature(object = "InputFrame",
	target="character", value="character")}:
      replace the detectors equal to the target with the value. Both
      \code{target} and \code{value} can be vectors of the same length,
      then the replace takes place iteratively.
    }
    \item{replaceSample}{\code{signature(object = "InputFrame",
	target="character", value="character")}:
      replace the samples equal to the target with the value. Both
      \code{target} and \code{value} can be vectors of the same length,
      then the replace takes place iteratively.
    }
    \item{show}{\code{signature(object="InputFrame")}: pretty print of the
      InputFrame instance.
    }
    \item{rightCensoring}{\code{signature(object="InputFrame",
	threshold="numeric")}: Right censoring the Ct value,  which
      targets the data points above a certain value
      (\code{threshold}). High Ct values (higher than 40 or 45 by the rule
      of thumb) are often not accurate and may indicate too weak
      expression. The function performs the right censoring on the data and
      set the value above the threshold as \code{NA} (by default) or a
      given value. See the example.
    }
    \item{coreData}{\code{signature(object="InputFrame")}: returns the
      data frame read from SDM file.
    }
    \item{coreData<-}{\code{signature(object="InputFrame")}: replace the
      data frame read from SDM file.
    }
    \item{Ct}{\code{signature(object="InputFrame")}: returns the Ct value
      of the SDM file.
    }
    \item{Ct}{\code{signature(object="InputFrame", value="numeric")}:
      replace the Ct value in the object with the new values, and return the object.}
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

## remove sample or detector
removeSample(sampdat, "Sample1")
removeDetector(sampdat, "Gene1")

## replace sample or detector
replaceSample(sampdat, "Sample1", "Sample0")
replaceDetector(sampdat, "Gene1", "PLCG1")

## right censoring the data
rightCensoring(sampdat, 35)
rightCensoring(sampdat, 35, 35)
}
\keyword{classes}