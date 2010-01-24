################################################################################
##
## This software is created by Molecular Genom Analysis Group
## Department of German Cancer Research Center in Heidelberg
##
##
##  AllGlobals.R
##  Created on: Oct 23, 2008
##      Author: Rudolf Biczok <r.biczok@dkfz-heidelberg.de>
##      Description: Some global definitions
##
################################################################################

############################################
## syntax of ddCt parameters
############################################
PARAM.SYNTAX            <- "^-([^-=\\s]+)=([^\\s]+)$"
SUB.PARAM.SYNTAX        <- "^([^-\\s]+)=([^\\s]+)$"

############################################
## syntax of the R commandline 
## programm Rscript
############################################
SYS.PARAM.SYNTAX        <- "^--([^=\\s]+)=?([^\\s]*)$"

############################################
## used columns of data
############################################
PRIMARY.INPUT.COLS      <- c("Sample","Detector","Ct");

