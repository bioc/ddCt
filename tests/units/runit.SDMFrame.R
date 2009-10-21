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

test.replaceVectorByEquality <- function() {
  vector <- c("USA","China","Russia","UK", "Germany", "Germany")
  perfectTarget <- c("Germany", "China", "Russia")
  perfectValue <- c("Deutschland", "VR China", "Russland")

  checkEquals(c("USA", "VR China", "Russland", "UK", "Deutschland", "Deutschland"),
              replaceVectorByEquality(vector=vector, target=perfectTarget,
                                      value=perfectValue))

  imperfectTarget <- c("Indian", "China", "Russia")
  imperfectValue <- c("Indien", "VR China", "Russland")

  checkEquals(c("USA", "VR China", "Russland", "UK", "Germany", "Germany"),
              replaceVectorByEquality(vector=vector, target=imperfectTarget,
                                      value=imperfectValue))
}

test.replaceDetector <- function() {
  newSDMframe <- replaceDetector(object=test.SDMFrame.data.sdm,
                                 target="Gene3",
                                 value="GeneDrei")
  checkEquals(detectorNames(newSDMframe),
              rep(rep(c("GeneDrei", "Gene2", "Gene1"), each=3),2))

  multiReplaceDetectorSDMframe <- replaceDetector(object=test.SDMFrame.data.sdm,
                                                   target=c("Gene3","Gene2"),
                                                  value=c("GeneDrei","GeneZwei"))
  checkEquals(detectorNames( multiReplaceDetectorSDMframe ),
              rep(rep(c("GeneDrei", "GeneZwei", "Gene1"), each=3),2))

  checkException(targetValueNotEquallyLong <- replaceDetector(object=test.SDMFrame.data.sdm,
                                                              target=c("Gene3","Gene2"),
                                                              value=c("GeneDrei")), silent=TRUE)
}

test.replaceSample <- function() {
  newSDMframe <- replaceSample(object=test.SDMFrame.data.sdm,
                               target="Sample2",
                               value="SampleZwei")
  checkEquals(sampleNames(newSDMframe),
              rep(c("Sample1", "SampleZwei"), each=9))

  multiReplaceSampleSDMframe <- replaceSample(object=test.SDMFrame.data.sdm,
                                                target=c("Sample1","Sample2"),
                                                value=c("SampleEins","SampleZwei"))
  checkEquals(sampleNames( multiReplaceSampleSDMframe ),
              rep(c("SampleEins", "SampleZwei"), each=9))

  checkException(targetValueNotEquallyLong <- replaceSample(object=test.SDMFrame.data.sdm,
                                                              target=c("Sample3","Sample2"),
                                                              value=c("SampleDrei")), silent=TRUE)
}

test.removeDetector <- function() {
  newSDMframe <- removeDetector(object=test.SDMFrame.data.sdm,
                                detector="Gene1")
  checkEquals(detectorNames(newSDMframe),
              rep(rep(c("Gene3", "Gene2"), each=3), 2))

  foreignNameSDMframe <- removeDetector(object=test.SDMFrame.data.sdm,
                                detector=c("Gene1","Gene4"))
  checkEquals(detectorNames(foreignNameSDMframe),
              rep(rep(c("Gene3", "Gene2"), each=3), 2))
}

test.removeSample <- function() {
  newSDMframe <- removeSample(object=test.SDMFrame.data.sdm,
                                sample="Sample1")
  checkEquals(sampleNames(newSDMframe),
              rep("Sample2", 9))
  
  foreignNameSDMframe <- removeSample(object=test.SDMFrame.data.sdm,
                                sample=c("Sample1","Sample4"))
  checkEquals(sampleNames(foreignNameSDMframe),
              rep("Sample2", 9))
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
