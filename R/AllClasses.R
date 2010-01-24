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
## utility classes
################################################################################

##----------------------------------------------------------------------------##
## Class ColMap
##----------------------------------------------------------------------------##
 
setClass("ColMap",
         representation(colmap="list")
         )

##----------------------------------------------------------------------------##
## Class ddCtParam
##----------------------------------------------------------------------------##
 
setClass("ddCtParam",
		 representation(type="character",default="ANY"))

##----------------------------------------------------------------------------##
## Visualization paramater set
##----------------------------------------------------------------------------##

setClass("errBarchartParameter",
         representation(exprsUndeterminedLabel="character"),
         prototype(exprsUndeterminedLabel="ND"))

################################################################################
## container classes
################################################################################

##----------------------------------------------------------------------------##
## Class ddCtExpression
##----------------------------------------------------------------------------##

setClass("ddCtExpression",
         contains="ExpressionSet")

##----------------------------------------------------------------------------##
## Class InputFrame
## Description: container for usable data in ddCt
##     columns: see AllGlobals
##----------------------------------------------------------------------------##
 
setClass("InputFrame",
         representation(coreData="data.frame",
                        files="character")
         )

################################################################################
## reader classes
################################################################################

##----------------------------------------------------------------------------##
## Class InputReader
## Description: abstract layer for valid input files
##----------------------------------------------------------------------------##
 
setClass("InputReader",
         representation(files="character",colmap="ColMap","VIRTUAL")
         )

##----------------------------------------------------------------------------##
## Class SDMReader
## Description: read a SDM file
##----------------------------------------------------------------------------##
 
setClass("SDMReader",
         contains=c("InputReader")
         )

##----------------------------------------------------------------------------##
## Class CSVReader
## Description: read a CSV file
##----------------------------------------------------------------------------##
 
setClass("CSVReader",
         contains=c("InputReader")
         )
