################################################################################
##
## This software is created by Molecular Genom Analysis Group
## Department of German Cancer Research Center in Heidelberg
##
##
##  AllClasses.R
##      Author: Rudolf Biczok <r.biczok@dkfz-heidelberg.de>
##              Jitao David Zhang <jitao_david.zhang@roche.com>
##      Description: ddCt classes
##
################################################################################

################################################################################
## utility classes
################################################################################

setClass("ColMap",
         representation=list(sample="character",
           feature="character",
           ct="character"),
         prototype=list(sample=DEFAULT.SAMPLE.COLNAME,
           feature=DEFAULT.FEATURE.COLNAME,
           ct=DEFAULT.CT.COLNAME))

setClass("ddCtParam",
         representation(type="character",default="ANY"))

setClass("errBarchartParameter",
         representation=list(exprsUndeterminedLabel="character"),
         prototype=list(exprsUndeterminedLabel="ND"))

################################################################################
## container classes
################################################################################

setClass("ddCtExpression",
         contains="ExpressionSet")

##----------------------------------------------------------------------------##
## Class InputFrame
## Description: container for usable data in ddCt
##     columns: see AllGlobals
##----------------------------------------------------------------------------##
 
setClass("InputFrame",
         representation(coreData="data.frame",
                        files="character"))

################################################################################
## reader classes
################################################################################

setClass("InputReader",
         representation(files="character",
                        colmap="ColMap"),
         contains=c("VIRTUAL"))

setClass("SDMReader",contains=c("InputReader"))
setClass("TSVReader",contains=c("InputReader"))
setClass("QuantStudioReader",contains=c("InputReader"))
