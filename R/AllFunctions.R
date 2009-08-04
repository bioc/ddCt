################################################################################
##
## This software is created by Molecular Genom Analysis Group
## Department of German Cancer Research Center in Heidelberg
##
##
##  ddCt.R
##  Created on: Mar 31, 2009
##      Author: Rudolf Biczok <r.biczok@dkfz-heidelberg.de>, Jitao David Zhang <j.zhang@dkfz.de>
##      Description: several internal functions used in the ddCt package
##    
################################################################################

##----------------------------------------##
## replaceNames used in class SDMFrame
##----------------------------------------##
replaceNames <- function(oldName ,targetName, newName) {
  if(is.null(targetName))
    stop("too few targets to replace")
  
  if(any(is.na(targetName)) || any(is.na(newName))) {
    warning("NAs will be ingonred")
    targetName <- as.character(na.omit(targetName))
    newName <- as.character(na.omit(newName))
  }
  
  for(i in 1:length(targetName)) {
    targetFound <- oldName == targetName[i]
    if(any(targetFound)) {
      if(is.na(newName[i]))
        stop(gettextf("missing new value for '%s'",targetName[i]))
      else
        oldName[targetFound] = newName[i]
    }
    else
      warning(gettextf("'%s' is not present in SDMFrame",targetName[i]))
  }
  return(oldName)
}

################################################################################
### alias functions
################################################################################

SDMFrame <- function(files) {
  return(new("SDMFrame",file=files));
}

readSDM <- function(files) {
  .Deprecated("SDMFrame",package="ddCt");
  return(SDMFrame(file=files))
}

################################################################################
### aux functions
################################################################################

na.sd <- function(x,...) {
  if (all(is.na(x))) return(NA)
  return (sd(x,...))
}
na.mad <- function(x,...) {
  if (all(is.na(x))) return(NA)
  return(mad(x,...))
}
na.mean <- function(x,...) {
  if (all(is.na(x))) return(NA)
  return(mean(x,...))
}
na.median <- function(x,...) {
  if (all(is.na(x))) return(NA)
  return(median(x,...))
}
levelfkt <- function(x) 2^(-x)
getDiff <- function(x) {
  if (any(is.na(x))|length(x) < 2 )
    y <- NA
  else
    y <- max(diff(sort(x)))
  return (y)
}
uniquePlate   <- function(x) {
  if (length(unique(x))!=1)
    warning(paste("g-s comb. on more than one plate:",paste(unique(x),collapse=",")), call.=FALSE)
  return (unique(x)[1])
}

################################################################################
### utility functions
################################################################################

getDir <- function(dir) {
  if(!file.exists(dir)) {
    dir.create(dir)
  }
  return(dir)
}

ddCtReport <- function(eSet, path, eSetLabel="bioRep") {
  table.path <- getDir(file.path(path,"table"))
  html.path <- getDir(file.path(path,"HTML"))
  
  elistWrite(eSet,
             file=file.path(path, sprintf("allValues_%s.txt",eSetLabel)))
  ad <- assayData(eSet)
  EE1 <- ad$exprs
  FF1 <- ad$level.err

  Ct <- round(ad$Ct, 2)
  lv <- round(EE1, 2)

  write.table(cbind(t(EE1), t(FF1)),
              file=file.path(table.path, sprintf("LevelPlusError_%s.txt", eSetLabel)),
              sep="\t", col.names=NA)

  write.table(lv,
              file=file.path(table.path, sprintf("level_matrix_%s.txt", eSetLabel)),
              sep="\t", col.names=NA)

  write.table(Ct,
              file=file.path(table.path, sprintf("Ct_matrix_%s.txt", eSetLabel)),
              sep="\t", col.names=NA)

  write.htmltable(cbind(rownames(lv),lv),title="Level",file=file.path(html.path,"level"))
  write.htmltable(cbind(rownames(Ct),Ct),title="Ct",file=file.path(html.path,"Ct"))
 
  dCtValues  <- round(ad$dCt,2)
  ddCtValues <- round(ad$ddCt,2)
  write.table(dCtValues,
              file=file.path(table.path,sprintf("dCt_matrix_%s.txt", eSetLabel)),
              sep="\t",col.names=NA)
  write.table(ddCtValues,
              file=file.path(table.path,sprintf("ddCt_matrix_%s.txt",eSetLabel)),
              sep="\t",col.names=NA)
  
  write.htmltable(cbind(rownames(dCtValues),dCtValues)  ,
                  title="dCt",
                  file=file.path(html.path,"dCt"))
  write.htmltable(cbind(rownames(ddCtValues),ddCtValues),
                  title="ddCt",
                  file=file.path(html.path,"ddCt"))
}

isGrep <- function(pattern,x,...) {
  res <- rep(FALSE, length(x))
  res[grep(pattern,x,...)] <- TRUE
  return(res)
}

################################################################################
### parameter parsing
################################################################################

getSysParams <- function() {
  sys.params <- list()
  for(arg in commandArgs()) {
    if(isTRUE(as.logical(grep(SYS.PARAM.SYNTAX, arg,perl=TRUE)))) {
      param.name  <- gsub(SYS.PARAM.SYNTAX, "\\1", arg, perl=TRUE)
      param.value <- gsub(SYS.PARAM.SYNTAX, "\\2", arg, perl=TRUE)
      sys.params[[param.name]] <- param.value
      ## all arguments after --args ar normal ddCt parameters
      if(sys.params[[param.name]] == "args")
        break
    }
  }
  return(sys.params)
}

getCmdParans <- function() {
  params.cmd <- list()
  for(arg in commandArgs(trailingOnly=TRUE)) {
    if(isTRUE(as.logical(grep(PARAM.SYNTAX, arg,perl=TRUE)))) {
      param.name  <- gsub(PARAM.SYNTAX, "\\1", arg, perl=TRUE)
      
      param.values <- NULL
      param.values.names <- NULL
      for(value in unlist(strsplit(gsub(PARAM.SYNTAX, "\\2", arg, perl=TRUE), ","))) {
        if(isTRUE(as.logical(grep(SUB.PARAM.SYNTAX, value,perl=TRUE)))) {
          param.values.names <- c(param.values.names, gsub(SUB.PARAM.SYNTAX, "\\1", value, perl=TRUE))
          param.values <- c(param.values, gsub(SUB.PARAM.SYNTAX, "\\2", value, perl=TRUE))
        }
        else
					param.values <- c(param.values, value)
      }
      names(param.values) <- param.values.names
      params.cmd[[param.name]] <- param.values
    }
  }
  return(params.cmd)
}

getConfParams <- function(conf.file) {
  params.conf <- NULL

  if(!isTRUE(as.logical(grep("^/.*",conf.file))))
    conf.file <- file.path(getwd(), conf.file)
  ## try to load config file
  out <- try(source(conf.file, local=TRUE), silent=TRUE)
  if(class(out) == "try-error") 
    warning(gettextf("Could not open ddCt config file \"%s\"", conf.file))
  
  return(params.conf)
}

getDefaultParams <- function(param.list) {
  params <- param.list
  for(param.name in names(param.list)) {
    params[[param.name]] <- param.list[[param.name]]@default
  }
  return(params)
}

getParams <- function(param.list) {
  params <- getDefaultParams(param.list)
  params.cmd <- getCmdParans()
  params.conf<- list()
  
  if(!is.null(params.cmd$confFile))
    params$confFile <- params.cmd$confFile
  
  if(!is.null(params$confFile))
    params.conf <- getConfParams(params$confFile)
  
  for(param.name in names(params)) {
    
    ## overwrite with config file values
    if(!is.null(params.conf[[param.name]]))
      params[[param.name]] <- params.conf[[param.name]]
    
    ## overwrite with command line values
    if(!is.null(params.cmd[[param.name]]))
      params[[param.name]] <- as(params.cmd[[param.name]], param.list[[param.name]]@type)
  }
  
  return(params)
}

checkParams <- function(params) {
  for(param.name in names(params)) {	
    if(!is.null(params[[param.name]]) && is.na(params[[param.name]]))
      stop(gettextf("value for parameter \"%s\" is missing", param.name))
  }
}

getExecMode <- function() {
  sys.params <- getSysParams()
  
  if(is.null(sys.params$file))
    return("source")
  else
    return("script")
}
