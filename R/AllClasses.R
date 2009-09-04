################################################################################
##
## This software is created by Molecular Genom Analysis Group
## Department of German Cancer Research Center in Heidelberg
##
##
##  AllClasses.R
##  Created on: Oct 23, 2008
##      Author: Rudolf Biczok <r.biczok@dkfz-heidelberg.de>
##              Jitao David Zhang <j.zhang@dkfz-heidelberg.de>
##      Description: ddCt classes
##
################################################################################

################################################################################
## Class ddCtExpression
################################################################################

setClass("ddCtExpression",
         contains="ExpressionSet")
 
################################################################################
## Class SDMFrame
## Description: read SDM file into a data.frame object with four mandary
##     columns: Sample, Detector, Ct and Platename. (case senstive)
################################################################################
 
setClass("SDMFrame",
         representation(coreData="data.frame",
                        files="character")
         )

################################################################################
## Class ddCtParam
################################################################################
 
setClass("ddCtParam",
		 representation(type="character",default="ANY"))

##----------------------------------------------------------------------------##
## Visualization paramater set
##----------------------------------------------------------------------------##
setClass("errBarchartParameter",
         representation(exprsUndeterminedLabel="character"),
         prototype(exprsUndeterminedLabel="ND"))
