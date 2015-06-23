################################################################################
##
## This software is created by Molecular Genom Analysis Group
## Department of German Cancer Research Center in Heidelberg
##
##
##  unit.R
##  Created on: Oct 23, 2008
##      Author: Rudolf Biczok <r.biczok@dkfz-heidelberg.de>
##      Description: RUnit test suit for ddCt classes
##
################################################################################

pkg <- "ddCt"
unit.path <- file.path(getwd(), "units")

unitMode <- function() {
  if(Sys.getenv("R_DEVELOP_MODE") == "TRUE")
    return("unit")
  else
    return("normal")
}

normalTest <- function() {
  testFile <- system.file("./extdata/Experiment2.txt", package=pkg)
  
  ## Basic SDMFrame
  sdm <- SDMFrame(testFile)
  
  ## ddCt
  x <- ddCtExpression(sdm,
                  calibrationSample="Sample3",
                  housekeepingGenes="Gene2")
  
  ## coerece as data frame
  y1 <- as(x, "data.frame")
  y2 <- elist(x)
  stopifnot(all.equal(y1,y2))
  
  ## visualization
  errBarchart(x)
}


## --- Setup ---

library(package=pkg, character.only=TRUE)

# put this in an enclosure so we can return early
(function() {
    if(unitMode() != "unit") {
        normalTest()
        return()
    }

    if(!require("RUnit", quietly=TRUE)) {
        stop("cannot run unit tests -- package RUnit is not available")
    }

    ## --- Testing ---
    cat("------------------- BEGIN UNIT TESTS ----------------------\n\n")

    ## --- Setup test suit ---
    testSuite <- defineTestSuite(name=paste(pkg, "unit testing"), dirs=unit.path)
    tests <- runTestSuite(testSuite)

    ## --- Setup report directory ---
    pathReport <- file.path(getwd(),"report")
    if (!file.exists(pathReport)) {
        dir.create(pathReport)
    }

    ## --- Reporting ---
    cat("------------------- UNIT TEST SUMMARY ---------------------\n\n")

    printTextProtocol(tests, showDetails=FALSE)
    printTextProtocol(tests, showDetails=FALSE,
                      fileName=file.path(pathReport, "summary.txt"))
    printTextProtocol(tests, showDetails=TRUE,
                      fileName=file.path(pathReport, "summary-detail.txt"))
    printHTMLProtocol(tests,
                      fileName=file.path(pathReport, "summary.html"))

    errors <- getErrors(tests)
    if(errors$nFail > 0 | errors$nErr > 0) {
        warning(paste("\n\nunit testing failed (#unit failures: ", errors$nFail,
                      ", #R errors: ",  errors$nErr, ")\n\n", sep=""))
    }

    cat("------------------- END OF UNIT TESTING -------------------\n\n")
})()
