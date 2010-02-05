################################################################################
##
## This software is created by Molecular Genom Analysis Group
## Department of German Cancer Research Center in Heidelberg
##
##
##  AllGenerics.R
##  Created on: Oct 23, 2008
##      Author: Jitao David Zhang <j.zhang@dkfz.de>, Rudolf Biczok <r.biczok@dkfz-heidelberg.de>
##      Description: ddCt generic methods
##
################################################################################


################################################################################
## Class ddCtExpression
################################################################################

setGeneric("errBarchart", function(object,by,...)	standardGeneric("errBarchart"))
setGeneric("elist", function(object,...) standardGeneric("elist"))
setGeneric("elistWrite", function(object,file,...) standardGeneric("elistWrite"))

setGeneric("level", function(object,...) standardGeneric("level"))
setGeneric("levelErr", function(object,...) standardGeneric("levelErr"))
setGeneric("Ct", function(object,...) standardGeneric("Ct"))
setGeneric("CtErr", function(object,...) standardGeneric("CtErr"))
setGeneric("dCt", function(object,...) standardGeneric("dCt"))
setGeneric("dCtErr", function(object,...) standardGeneric("dCtErr"))
setGeneric("ddCt", function(object,...) standardGeneric("ddCt"))
setGeneric("ddCtErr", function(object,...) standardGeneric("ddCtErr"))
setGeneric("numberCt", function(object,...) standardGeneric("numberCt"))
setGeneric("numberNA", function(object,...) standardGeneric("numberNA"))

################################################################################
## Class InputFrame and InputReader
################################################################################

setGeneric("detectorNames",function(object) standardGeneric("detectorNames"))
setGeneric("detectorNames<-",function(object, value) standardGeneric("detectorNames<-"))

setGeneric("sampleNames",function(object) standardGeneric("sampleNames"))
setGeneric("sampleNames<-",function(object, value) standardGeneric("sampleNames<-"))

setGeneric("uniqueDetectorNames",function(object) standardGeneric("uniqueDetectorNames"))
setGeneric("uniqueDetectorNames<-",function(object, target, value) standardGeneric("uniqueDetectorNames<-"))

setGeneric("uniqueSampleNames",function(object) standardGeneric("uniqueSampleNames"))
setGeneric("uniqueSampleNames<-",function(object, target, value) standardGeneric("uniqueSampleNames<-"))

setGeneric("removeSample", function(object, sample) standardGeneric("removeSample"))
setGeneric("removeDetector", function(object, detector) standardGeneric("removeDetector"))

setGeneric("replaceSample", function(object, target, value) standardGeneric("replaceSample"))
setGeneric("replaceDetector", function(object, target, value) standardGeneric("replaceDetector"))

setGeneric("fileNames",function(object) standardGeneric("fileNames"))

setGeneric("ddCtExpression",function(object, warningStream, algorithm, calibrationSample, housekeepingGenes, type, sampleInformation, toZero, efficiencies, efficiencies.error) standardGeneric("ddCtExpression"))

setGeneric("rightCensoring", function(object, threshold, value) standardGeneric("rightCensoring"))
           
## private
setGeneric("ddCtExec",function(object, calibrationSample, housekeepingGenes, type, sampleInformation, toZero) standardGeneric("ddCtExec"))
setGeneric("ddCtWithEExec",function(object, calibrationSample, housekeepingGenes,  type, sampleInformation, efficiencies, efficiencies.error) standardGeneric("ddCtWithEExec"))
setGeneric("coreData", function(object) standardGeneric("coreData"))
setGeneric("coreData<-", function(object, value) standardGeneric("coreData<-"))
setGeneric("headtailPrint", function(object,...) standardGeneric("headtailPrint"))

##----------------------------------------------------------------------------##
## Visualization parameter
##----------------------------------------------------------------------------##
setGeneric("exprsUndeterminedLabel", function(object) standardGeneric("exprsUndeterminedLabel"))

################################################################################
## Class InputReader
################################################################################

setGeneric("InputFrame",function(object) standardGeneric("InputFrame"))
setGeneric("readRawData",function(object, file.name) standardGeneric("readRawData"))
setGeneric("validate",function(object, file.name) standardGeneric("validate"))
