\name{SDMFrame}
\alias{SDMFrame}
\alias{readSDM}
\alias{InputFrame,InputReader-method}
\alias{rightCensoring,InputFrame,numeric-method}
\title{Read an SDM file}
\description{Read an SDM file: Data Output File for SDS, Version 2.1}
\usage{
  SDMFrame(file)
  readSDM(file)
}
\arguments{
  \item{file}{Character vector of filenames}
}

\details{This function reads the data given in the colums
  'Detector','Sample' and 'Ct' of the specified SDM output file(s) and
  stores them as a data.frame. An additional column including the
  respective filename is added.

  This function is a wrapper for the SDMFrame constructor
}

\value{A object of class \code{\linkS4class{SDMFrame}}}

\author{ Rudolf Biczok \url{mailto:r.biczok@dkfz.de}}

\examples{
## read a SDM file
sampdat <- SDMFrame(system.file("extdata", "Experiment1.txt",
                                package="ddCt"))
## you can also write 
## sampdat <- new("SDMFrame",system.file("extdata", "Experiment1.txt",
##                                       package="ddCt"))

## or with
## sampdat <- readSDM(system.file("extdata", "Experiment1.txt",
##                                package="ddCt"))

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

}


