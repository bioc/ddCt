################################################################################
##
## This software is created by Molecular Genom Analysis Group
## Department of German Cancer Research Center in Heidelberg
##
##
##  runit.SDMFrame.R
##  Created on: Oct 23, 2008
##      Author: Rudolf Biczok <r.biczok@dkfz-heidelberg.de>
##      Description: Unit for SDMFrame class
##
################################################################################

### --- Test setup ---

load(file.path(unit.path,"runit.SDMFrame.RData"))
load(file.path(unit.path,"runit.ddCtExpression.RData"))
 
### --- Test functions ---


test.SDMFrame.initialize <- function()
{
  sdm <- SDMFrame(system.file("extdata", "Experiment1.txt", package="ddCt"))
  checkEquals(sdm@coreData, test.SDMFrame.data.sdm@coreData)
}

test.SDMFrame.ddCtExpression <- function() {
  result <- ddCtExpression(test.SDMFrame.data.sdm,
                           calibrationSample="Sample1",
                           housekeepingGenes=c("Gene1","Gene2"))
  checkEquals(assayData(test.ddCtExpression.data.ddCt), assayData(result))
}

test.SDMFrame.sampleNames <- function()
{
  checkEquals(test.SDMFrame.data.sdm@coreData[,"Sample"], sampleNames(test.SDMFrame.data.sdm))
}

test.SDMFrame.detectorNames <- function()
{
  checkEquals(test.SDMFrame.data.sdm@coreData[,"Detector"], detectorNames(test.SDMFrame.data.sdm))
}

test.SDMFrame.uniqueSampleNames <- function()
{
  checkEquals(unique(test.SDMFrame.data.sdm@coreData[,"Sample"]), uniqueSampleNames(test.SDMFrame.data.sdm))
}

test.SDMFrame.uniqueDetectorNames <- function()
{
  checkEquals(unique(test.SDMFrame.data.sdm@coreData[,"Detector"]), uniqueDetectorNames(test.SDMFrame.data.sdm))
}

test.SDMFrame.fileNames <- function()
{
  checkEquals(unique(test.SDMFrame.data.sdm@file), fileNames(test.SDMFrame.data.sdm))
}

assign("test.SDMFrame.uniqueSampleNames<-", function()
{
  sdm <- test.SDMFrame.data.sdm
  newNames <- oldNames <- unique(sdm@coreData[,"Sample"])
  newNames[2] <- "NewSample"
  uniqueSampleNames(sdm,oldNames[2]) <- newNames[2]
  
  checkEquals(unique(sdm@coreData[,"Sample"]), newNames)
})

assign("test.SDMFrame.uniqueDetectorNames<-", function()
{
  sdm <- test.SDMFrame.data.sdm
  newNames <- oldNames <- unique(sdm@coreData[,"Detector"])
  newNames[1] <- "NewGene"
  uniqueDetectorNames(sdm,oldNames[1]) <- newNames[1]
  
  checkEquals(unique(sdm@coreData[,"Detector"]), newNames)
})


assign("test.SDMFrame.[", function()
{
  checkEquals(unique(test.SDMFrame.data.sdm@coreData[,"Ct"]), test.SDMFrame.data.sdm[,"Ct"])
})

## currently disabled
##assign("test.SDMFrame.[[", function()
##{
##  checkEquals(unique(test.SDMFrame.data.sdm@coreData[[1]]), test.SDMFrame.data.sdm[[1]])
##})


assign("test.SDMFrame.$", function()
{
  checkEquals(unique(test.SDMFrame.data.sdm@coreData$Ct), test.SDMFrame.data.sdm$Ct)
})
