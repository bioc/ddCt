################################################################################
##
## This software is created by Molecular Genom Analysis Group
## Department of German Cancer Research Center in Heidelberg
##
##
##  output.R
##  Created on: Oct 23, 2008
##      Author: Rudolf Biczok <r.biczok@dkfz-heidelberg.de>
##      Description: Output utility
##
################################################################################

## xtable method
## TODO: change S3 methods from "ddCt" to "ddCtExpression"
xtable.ddCt <- function(expSet,name="exprs",na.threshold,...) {
  ad <- assayData(expSet)
  if(missing(name))
    name <- "exprs"
  if(!name %in%  names(ad)) 
    stop("Data entry '", name,"' not found! It should be of the following:\n\t",
         paste(names(ad), collapse=","),"\n")

  pos <- match(name,names(ad))
  nas <- t(ad$numberNA)
  numbers <- t(ad$number)
  
  exp <- t(ad[[pos]])

  if(!missing(na.threshold)) {
    if(!is.numeric(na.threshold)) {stop("'na.threshold' must be numeric!\n")}
    if(na.threshold >100 | na.threshold < 0) {
      stop ("'na.threshold' must range between 0 and 100 - as percentage!\n")
    }
    highna <- nas >= numbers * na.threshold / 100
    exp[highna] <- NA
  }
  
  xtable(exp,...)
}

## write ddCt into file
writeLines.ddCt <- function(expSet, con=stdout(),sep="\n",type="html",name="exprs",...) {
  txt <- print(xtable.ddCt(expSet,
                           xname=name,...),type=type)
  writeLines(txt, con=con,sep=sep)
}

## write tab-delimited, not quoted and no-row-name TSV files
writeSimpleTabCsv <- function(x, file="",...) {
  Call <- match.call(expand.dots=TRUE)
  for(argname in c("row.names", "sep", "quote")) {
    if (!is.null(Call[[argname]]))
      warning(gettextf("attemp to set '%s' ignored", argname), domain=NA)
  }
  Call$row.names <- FALSE
  Call$sep <- "\t"
  Call$quote <- FALSE
  Call[[1L]] <- as.name("write.table")
  eval.parent(Call)
}

## write html table, contributed by Andreas Buness
##------------------------------------------------------------
## writes a html table 
## x:        data frame
## filename: character with the name of the output file
## sortby:   name of a column in data frame din
##------------------------------------------------------------
write.htmltable <- function (x, file, title="", sortby=NULL, decreasing=TRUE, open="wt") {

  if(!is.null(sortby)) {
    if(!sortby %in% colnames(x))
      stop(paste("Invalid argument \"sortby\": could not find a column in data frame x with name", sortby))
    soby = x[, sortby]
    if(!is.numeric(soby))
      stop("Invalid argument \"sortby\": column is not numeric")
    x = x[order(soby, decreasing=decreasing), ]
  }
  
  outfile <- file(paste(file, ".html", sep=""), open=open)
  cat("<html>", "<STYLE>", 
      "<!--TD { FONT-FAMILY: Helvetica,Arial; FONT-SIZE: 14px;}-->",
      "<!--H1 { FONT-FAMILY: Helvetica,Arial; FONT-SIZE: 22px;}-->",
      "</STYLE>", "<head>", paste("<TITLE>", title, "</TITLE>", sep=""),
      "</head>", "<body bgcolor=#ffffff>", file = outfile, sep = "\n")
  
  if (title!="") 
    cat("<CENTER><H1 ALIGN=\"CENTER\">", title, " </H1></CENTER>\n", 
        file = outfile, sep = "\n")
  cat("<CENTER> \n", file = outfile)

  cat("<TABLE BORDER=0>", file = outfile, sep = "\n")
  
  cat("<TR>", file = outfile)
  for (j in 1:ncol(x))
    cat("<TD BGCOLOR=\"", c("#e0e0ff", "#d0d0f0")[j%%2+1], "\"><B>", colnames(x)[j], "</B></TD>\n", sep="", file = outfile)
  cat("</TR>", file = outfile)
  
  for (i in 1:nrow(x)) {
    cat("<TR>", file = outfile)
    for (j in 1:ncol(x)) {
      txt <- as.character(x[i,j])
      if (length(grep("^http:", txt))>0) {
        txt <- sub(";$", "", txt)              ## remove trailing semicolon
        s <- unlist(strsplit(txt, "[/?=]"))    ## split out last part of URL
        txt <- paste("<A HREF=\"", txt, "\" TARGET=\"z\">", s[length(s)], "</A>", sep="")
      }
      cat("<TD BGCOLOR=\"", c("#e0e0ff", "#d0d0f0", "#f0f0ff", "#e0e0f0")[i%%2*2+j%%2+1], "\">",
          txt, "</TD>\n", sep="", file = outfile)
    }
    cat("</TR>", file = outfile)
  }

  cat("</TABLE></CENTER>", "</body>", "</html>", sep = "\n", file = outfile)
  close(outfile)
}

