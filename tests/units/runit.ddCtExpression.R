################################################################################
##
## This software is created by Molecular Genom Analysis Group
## Department of German Cancer Research Center in Heidelberg
##
##
##  runit.ddCtExpression.R
##  Created on: Oct 23, 2008
##      Author: Rudolf Biczok <r.biczok@dkfz-heidelberg.de>
##      Description: Unit for ddCtExpression class
##
################################################################################

### --- Test setup ---

load(file.path(unit.path,"runit.ddCtExpression.RData"))
 
### --- Test functions ---

test.ddCtExpression.level <- function() {
  checkEquals(level(test.ddCtExpression.data.ddCt),
              exprs(test.ddCtExpression.data.ddCt))
}

test.ddCtExpression.levelErr <- function() {
  checkEquals(levelErr(test.ddCtExpression.data.ddCt),
              assayData(test.ddCtExpression.data.ddCt)$level.err)
}

test.ddCtExpression.Ct <- function() {
  checkEquals(Ct(test.ddCtExpression.data.ddCt),
              assayData(test.ddCtExpression.data.ddCt)$Ct)
}

test.ddCtExpression.CtErr <- function() {
  checkEquals(CtErr(test.ddCtExpression.data.ddCt),
              assayData(test.ddCtExpression.data.ddCt)$Ct.err)
}

test.ddCtExpression.dCt <- function() {
  checkEquals(dCt(test.ddCtExpression.data.ddCt),
              assayData(test.ddCtExpression.data.ddCt)$dCt)
}

test.ddCtExpression.dCtErr <- function() {
  checkEquals(dCtErr(test.ddCtExpression.data.ddCt),
              assayData(test.ddCtExpression.data.ddCt)$dCt.err)
}

test.ddCtExpression.ddCt <- function() {
  checkEquals(ddCt(test.ddCtExpression.data.ddCt),
              assayData(test.ddCtExpression.data.ddCt)$ddCt)
}

test.ddCtExpression.ddCtErr <- function() {
  checkEquals(ddCtErr(test.ddCtExpression.data.ddCt),
              assayData(test.ddCtExpression.data.ddCt)$ddCt.err)
}

test.ddCtExpression.numberNA <- function() {
  checkEquals(numberNA(test.ddCtExpression.data.ddCt),
              assayData(test.ddCtExpression.data.ddCt)$numberNA)
}

test.ddCtExpression.numberCt <- function() {
  checkEquals(numberCt(test.ddCtExpression.data.ddCt),
              assayData(test.ddCtExpression.data.ddCt)$number)
}
