################################################################################
##
## This software is created by Molecular Genom Analysis Group
## Department of German Cancer Research Center in Heidelberg
##
##
##  ddCtSplit.R
##  Created on: Oct 23, 2008
##      Author: Rudolf Biczok <r.biczok@dkfz-heidelberg.de>
##      Description: TODO
##
################################################################################

ddCtSplit <- function(raw.table,
                      calibrationSample,
                      housekeepingGenes,
                      type="mean",
                      sampleInformation = NULL,
                      toZero = FALSE,
                      filename = "warning.output.txt",
                      split.field = NULL,
                      subset = NULL) {

  if (is.null(sampleInformation)) {
    return(ddCt(raw.table,
                calibrationSample=calibrationSample,
                housekeepingGenes=housekeepingGenes,
                type = type,
                toZero = toZero,
                filename = filename))
  }
  
  ## colnames of pd all lower-cased
  pd <- pData(sampleInformation)
  if(!("Sample" %in% colnames(pd) && "Sample" %in% colnames(raw.table))) {
    stop("Both sampleInformation and raw table must contain column 'Sample'!\n")
  }
  if(any(duplicated(pd$Sample))) {
    warning("Sample information has duplicated rows, they shall be removed!\n ")
    pd <- subset(pd, !duplicated(pd$Sample))
  }

  pdcolnames <- tolower(colnames(pd))  
  pdi <- function(x) match(tolower(x), pdcolnames)
  pdic.RO <- function(x) pd[,match(tolower(x),pdcolnames)]

  subset <- eval(substitute(subset), pd)
  if(!is.null(subset)) {
    raw.table.subset <- raw.table$Sample %in% pd$Sample[subset]
    
    pd <- pd[subset,]
    raw.table <- raw.table[raw.table.subset,]
    stopifnot(nrow(pd)>0)
  }

  cip <- match(raw.table$Sample,pd$Sample)
  ## "geneid" (case insensitive) as conservative words
  ## in case there's a column named "geneid", its content shall replace the "Sample" column
  if("geneid" %in% pdcolnames) {  
    ns <- pdic.RO("geneid")
    raw.table$Sample <- ns[cip]
    pd[,pdi("Sample")] <- pdic.RO("geneid")
  }

  m <- raw.table$Sample %in% pd$Sample
  if(!all(m)) {
    noin <- which(!m)
    stop("The following samples are not annotated:", paste(raw.table$Sample[noin], collapse=","),"!\n")
  }
  
  ## split raw.table
  if (is.null(split.field)) {
    return(ddCt(raw.table,
                calibrationSample=calibrationSample,
                housekeepingGenes=housekeepingGenes,
                sampleInformation <- new("AnnotatedDataFrame", data=pd, varMetadata=varMetadata(sampleInformation)),
                type = type,
                toZero = toZero,
                filename = filename))
  } else {
    if(!all(tolower(split.field) %in% pdcolnames)) {
      notin <- which(!tolower(split.field) %in% pdcolnames)
      stop("The following columns are not in the sampleInformation:", paste(split.field[notin],collapse=","),"!\n")
    }
    
    splits <- pd[,pdi(split.field)]
    split.pd <- split(pd, splits)
    ## recall the cip is the mapping of raw.table to pd on sample.
    split.raw.table <- split(raw.table, splits[cip])
  }
  
  var <- varMetadata(sampleInformation)
  result <- mapply(function(x,y) {
    a <- new("AnnotatedDataFrame", data=y, varMetadata=var)
    ddCt(x,
         calibrationSample=calibrationSample,
         housekeepingGenes=housekeepingGenes,
         type = type,
         sampleInformation = a,
         toZero = toZero,
         filename = filename)
  }, split.raw.table, split.pd)
  return(result)

}
