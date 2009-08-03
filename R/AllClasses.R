################################################################################
##
## This software is created by Molecular Genom Analysis Group
## Department of German Cancer Research Center in Heidelberg
##
##
##  AllClasses.R
##  Created on: Oct 23, 2008
##      Author: Rudolf Biczok <r.biczok@dkfz-heidelberg.de>
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
 
setClass("SDMFrame",representation(coreData="data.frame",files="character",withoutPath="logical"))

################################################################################
## Class ddCtParam
################################################################################
 
setClass("ddCtParam",
		 representation(type="character",default="ANY"))
