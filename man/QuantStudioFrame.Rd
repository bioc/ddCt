\name{QuantStudioFrame}
\alias{readQuantStudio}
\alias{QuantStudioFrame}
\title{Read QuantStudio file(s)}
\description{Read QuantStudio file(s)}
\usage{
QuantStudioFrame(file)
readQuantStudio(file)
}
\arguments{
  \item{file}{Character vector of filenames}
}

\details{This function reads the data given in the QuantStudio output files}

\value{A object of class \code{\linkS4class{SDMFrame}}}

\author{ Jitao David Zhang \url{mailto:jitao_david.zhang@roche.com}}

\examples{
sampdat <- QuantStudioFrame(system.file("extdata", c("QuantStudio_File1.txt", "QuantStudio_File2.txt"),package="ddCt"))

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
uniqueDetectorNames(sampdat) <- c(Hs00559368_m1="Gene1", Hs02576168_g1="Gene2")
uniqueDetectorNames(sampdat)

}


