################################################################################
##
## This software is created by Molecular Genom Analysis Group
## Department of German Cancer Research Center in Heidelberg
##
##
##  ddCt.R
##  Created on: Oct 23, 2008
##      Author errBarchart: Rudolf Biczok <r.biczok@dkfz-heidelberg.de>
##      Description: commandline based ddCt script
##
################################################################################

################################################################################
### Libraries
################################################################################

library(ddCt)

################################################################################
### valid parameters
################################################################################

param.list <- list(
                   ## locations
                   loadPath              = new("ddCtParam", "character", getwd()),
                   savePath		 = new("ddCtParam", "character", getwd()),
                   ## input files         
                   confFile		 = new("ddCtParam", "character", NULL),
                   inputFile		 = new("ddCtParam", "character", NA),
                   sampleAnnotationFile  = new("ddCtParam", "character", NULL),
                   columnForGrouping	 = new("ddCtParam", "character", NULL),
                   onlyFromAnnotation	 = new("ddCtParam", "logical", FALSE),
                   ## transform
                   geneAlias		 = new("ddCtParam", "character", NULL),
                   sampleAlias		 = new("ddCtParam", "character", NULL),
                   ## housekeeping
                   referenceSample       = new("ddCtParam", "character", NA),
                   referenceGene         = new("ddCtParam", "character", NA),
                   ## threshold
                   threshold		 = new("ddCtParam", "numeric", 40),
                   ## core settings 
                   mode			 = new("ddCtParam", "character", "median"),
                   plotMode		 = new("ddCtParam", "character", c("level","Ct")),
                   algorithm		 = new("ddCtParam", "character", "ddCt"),
                   ## efficiencies
                   efficiencies		 = new("ddCtParam", "character", NULL),
                   efficienciesError	 = new("ddCtParam", "character", NULL),
                   ## remaining
                   genesRemainInOutput	 = new("ddCtParam", "character", NULL),
                   samplesRemainInOutput = new("ddCtParam", "character", NULL),
                   genesNotInOutput      = new("ddCtParam", "character", NULL),
                   samplesNotInOutput    = new("ddCtParam", "character", NULL),
                   ## grouping
                   groupingBySamples     = new("ddCtParam", "logical", TRUE),
                   ## plot
                   plotPerObject	 = new("ddCtParam", "logical", TRUE),
                   groupingForPlot	 = new("ddCtParam", "character", NULL),
                   groupingColor	 = new("ddCtParam", "character", c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00")),
                   cutoff		 = new("ddCtParam", "numeric", NULL),
                   brewerColor		 = new("ddCtParam", "character", c("Set3","Set1","Accent","Dark2","Spectral","PuOr","BrBG")),
                   legend		 = new("ddCtParam", "logical", TRUE),
                   ## ttest
                   ttestPerform		 = new("ddCtParam", "logical", FALSE),
                   ttestCol		 = new("ddCtParam", "character", NULL),
                   pairsCol		 = new("ddCtParam", "character", NULL),
                   samplesRemainInTTest	 = new("ddCtParam", "character", NULL),
                   samplesNotInTTest	 = new("ddCtParam", "character", NULL),
                   samplesRemainInCor	 = new("ddCtParam", "character", NULL),
                   samplesNotInCor	 = new("ddCtParam", "character", NULL))

################################################################################
### Functions
################################################################################

##----------------------------------------##
## myFoldChange
##----------------------------------------##
myFoldChange <- function(x){
  x <- as.numeric(x)
  if(length(x) == 2 ){
    resu <- 2^(x[1]-x[2])
  }else if( length(x) == 1  ){
    resu <- 2^x
  }else{
    stop("unexpected situation in myFoldChange")
  }
  return(resu)
}

##----------------------------------------##
## readAll - reads all required input
##----------------------------------------##
prepareInput <- function(params, sdmframe, datadir) {
  input <- list()
  input$CtData = sdmframe
  inputCt <- Ct(input$CtData)
  
  if (!is.null(params$threshold)){
    A <- which(inputCt>= params$threshold)
    if (length(A) > 0)
      input$CtData[ which(inputCt>= params$threshold),"Ct"] <- NA
  }
  
  if (! is.null(params$geneAlias)) detectorNames(input$CtData) <- params$geneAlias[detectorNames(input$CtData)]
  if (! is.null(params$sampleAlias)) sampleNames(input$CtData) <- params$sampleAlias[sampleNames(input$CtData)]
  
  
  if (! is.null(params$sampleAnnotationFile)){
    info <- file.path(datadir, params$sampleAnnotationFile)
    input$sampleInformation <- read.AnnotatedDataFrame(info,header=TRUE,sep="\t", row.names=NULL)
  } else {
    input$sampleInformation <- new("AnnotatedDataFrame",
                             data=data.frame(Sample=uniqueSampleNames(input$CtData)),
                             varMetadata=data.frame(labelDescription=c("unique identifier"),row.names=c("Sample")))
  }
  
  if(params$onlyFromAnnotation && !is.null(params$sampleAnnotationFile)){
    A <- sampleNames(input$CtData) %in% pData(input$sampleInformation)[,"Sample"]
    warning( paste("This samples will be removed:",paste(unique(input$CtData[!A,"Sample"]),collapse=", ")))
    input$CtData <- input$CtData[A,]
  }

  return(input)
}

##----------------------------------------##
## doTables - HTML and Table output
##----------------------------------------##
doTables <- function(params, result, savedir) {
  htmlName   <- "HTML"
  tablesName <- "Tables"
  
  html.path <- ddCt:::getDir(file.path(savedir, htmlName))
  table.path <- ddCt:::getDir(file.path(savedir, tablesName))
  
  if(!is.null(params$sampleAnnotationFile) & !is.null(params$columnForGrouping))
    result <- result[,order(pData(result)[,params$columnForGrouping])]
  
  elistWrite(result,file=file.path(savedir, paste("allValues",".csv",sep="")),sep="\t",row.names=FALSE)
  EE1 <- level(result)
  FF1 <- levelErr(result)
  
  Ct  <- round(Ct(result),2)
  lv  <- round(EE1,2)
  
  write.table(cbind(t(EE1),t(FF1)),file=file.path(table.path,"levelPlusError.txt"),sep="\t",col.names=NA)
  write.table(lv,file=file.path(table.path,"level_matrix.txt"),sep = "\t", col.names = NA)
  write.table(Ct,file=file.path(table.path,"Ct_matrix.txt"),sep="\t", col.names = NA)
  
  write.htmltable(cbind(rownames(lv),lv),title="Level",file=file.path(html.path,"level"))
  write.htmltable(cbind(rownames(Ct),Ct),title="Ct",file=file.path(html.path,"Ct"))
  
  
  if(is.null(params$efficiencies)){
    dCtValues  <- round(dCt(result),2)
    ddCtValues <- round(ddCt(result),2)
    write.table(dCtValues,file=file.path(table.path,"dCt_matrix.txt"),  sep="\t", col.names=NA)
    write.table(ddCtValues,file=file.path(table.path,"ddCt_matrix.txt"),sep="\t",col.names=NA)
    
    write.htmltable(cbind(rownames(dCtValues),dCtValues)  ,title="dCt",file=file.path(html.path,"dCt"))
    write.htmltable(cbind(rownames(ddCtValues),ddCtValues),title="ddCt",file=file.path(html.path,"ddCt"))
  }
}

##----------------------------------------##
## doPlot
##----------------------------------------##
doPlot <- function(params, result, savedir) {
  if(!is.null(params$groupingForPlot)) 
    GFP1 <- pData(result)[,params$groupingForPlot]
  
  for(KindOfPlot in params$plotMode) {
    
    if(KindOfPlot=="level") {
      EE1 <- level(result)
      FF1 <- levelErr(result)
      theTitle <- "level"
    }
    
    if(KindOfPlot=="ddCt") {
      EE1 <- ddCt(result)
      FF1 <- ddCtErr(result)
      theTitle <- "ddCt"
    }
    
    if(KindOfPlot=="dCt") {
      EE1 <- dCt(result)
      FF1 <- dCtErr(result)
      theTitle <- "dCt"
    }
    
    if(KindOfPlot=="Ct") {
      EE1 <- Ct(result)
      FF1 <- CtErr(result)
      theTitle <- "Ct"
    }
    
#################
## order the set
#################
    if (!is.null(params$geneAlias)) {
      EE2  <- EE1[match(params$geneAlias,rownames(EE1)),,drop=FALSE]
      FF2  <- FF1[match(params$geneAlias,rownames(EE1)),,drop=FALSE]
    } else {
      EE2 <- EE1
      FF2 <- FF1
    }
    
    if (!is.null(params$sampleAlias)) {
      EE  <- EE2[,match(params$sampleAlias,colnames(EE2)),drop=FALSE]
      FF  <- FF2[,match(params$sampleAlias,colnames(EE2)),drop=FALSE]
      if(!is.null(params$groupingForPlot)) GFP<- GFP1[match(params$sampleAlias,colnames(EE2))]
    } else {
      EE <- EE2
      FF <- FF2
      if(!is.null(params$groupingForPlot)) GFP<- GFP1
    }
    
####################
## Reducing the set
####################
    if (!is.null(params$genesRemainInOutput)){
      Gred <- (rownames(EE) %in% params$genesRemainInOutput)
    } else {
      Gred <- !(rownames(EE) %in% params$genesNotInOutput)
    }

    
    if (!is.null(params$samplesRemainInOutput)){
      Sred <- (colnames(EE) %in% params$samplesRemainInOutput)
    } else {
      Sred <- !(colnames(EE) %in% params$samplesNotInOutput)
    }
    
    EEN <- EE[Gred,Sred,drop=FALSE]
    FFN <- FF[Gred,Sred,drop=FALSE]
    if(!is.null(params$groupingForPlot)) GFPN <- as.factor(as.character(GFP[Sred]))
    
#############
## the color
#############
    COLORS <- c()
    for (i in 1:length(params$brewerColor))
      COLORS <- c(COLORS,brewer.pal(brewer.pal.info[params$brewerColor[i],]$maxcolors,params$brewerColor[i]))
    if(params$groupingBySamples) {
      THECO  <- COLORS[1:sum(Sred)]
    } else {
      THECO  <- COLORS[1:sum(Gred)]
    }
    
############
## plotting
############    
    if(params$plotPerObjec){
      pdf(w=15,h=15,file=file.path(savedir, paste(theTitle,"Result",".pdf",sep="")))
      if(params$groupingBySamples){
        for(k in 1:dim(EEN)[1]){
          EENN <- EEN[k,,drop=FALSE]
          FFNN <- FFN[k,,drop=FALSE]
          suppressWarnings(barploterrbar(EENN,EENN-FFNN,EENN+FFNN,barcol=THECO,legend=params$legend,columnForDiffBars=params$groupingBySamples,theCut=params$cutoff,ylab=theTitle))
        }
      }else{
        for(k in 1:dim(EEN)[2]){
          EENN <- EEN[,k,drop=FALSE]
          FFNN <- FFN[,k,drop=FALSE]
          suppressWarnings(barploterrbar(EENN,EENN-FFNN,EENN+FFNN,barcol=THECO,legend=params$legend,columnForDiffBars=params$groupingBySamples,theCut=params$cutoff,ylab=theTitle))
        }
      }
      if(params$groupingBySamples & !is.null(params$groupingForPlot)){
        for(k in 1:dim(EEN)[1]){
          EENN <- EEN[k,,drop=FALSE]
          FFNN <- FFN[k,,drop=FALSE]
          suppressWarnings(barploterrbar(EENN,EENN-FFNN,EENN+FFNN,barcol=params$groupingColor,legend=params$legend,columnForDiffBars=params$groupingBySamples,theCut=params$cutoff,ylab=theTitle,las=2,names.arg=colnames(EENN),main=rownames(EENN),groups=GFPN))
        }
      }
      dev.off()
    }else{
      suppressWarnings(barploterrbar(EEN,EEN-FFN,EEN+FFN,barcol=THECO,legend=params$legend,las=2,columnForDiffBars=params$groupingBySamples,theCut=params$cutoff,ylab=theTitle))
      dev.copy(pdf,w=15,h=15,file=file.path(savedir, paste(theTitle,"Result",".pdf",sep="")))
      dev.off()
    }
  }
}

##----------------------------------------##
## doCorrelation
##----------------------------------------##
doCorrelation <- function(params, result, savedir) {
  A <- params$referenceGene 
  if (length(A) > 1) { 
    if (!is.null(params$samplesRemainInCor)) {
      corRed <- (rownames(pData(result)) %in% params$samplesRemainInCor)
    } else {
      corRed <- !(rownames(pData(result)) %in% params$samplesNotInCor)
    }
    result2 <- Ct(result[,corRed])
    
    K <- combn(1:length(A),2)
    U <- ncol(K)
    pdf(file=file.path(savedir, paste("correlationResult",".pdf",sep="")))
    par(mfrow=c(ceiling(sqrt(U)),ceiling(sqrt(U))))
    for (i in 1:U) {
      Gen1 <- A[K[1,i]] ## first Housekeepinggene
      Gen2 <- A[K[2,i]] ## sceond Housekeepinggene
      BART <- cor(result2[Gen1,],result2[Gen2,],use="pairwise.complete.obs")
      if(!is.na(BART)) {
        plot(result2[Gen1,],result2[Gen2,],xlab=Gen1,ylab=Gen2,pch="*",col="red", main=paste("correlation:",BART))
      } else {
        gen1.allna <- all(is.na(result2[Gen1,]))
        warn <- sprintf("No correlation available becaouse %s is undetermined in all sample'\n",ifelse(gen1.allna, Gen1, Gen2))

        plot.new()
        text(0.5,0.5, warn)
      }
    }
    dev.off()
  } else {
    if (!is.null(params$samplesRemainInCor)){
      corRed <- (rownames(pData(result)) %in% params$samplesRemainInCor)
    } else {
      corRed <- !(rownames(pData(result)) %in% params$samplesNotInCor)
    }
    result2 <- assayData(result[,corRed])$Ct

    pdf(file=file.path(savedir, paste("expressionHKGeneResult",".pdf",sep="")))
    plot(result2[A,],pch="*",col="red", main="Expression HK Gene")
    dev.off()
  }
  
  if(params$ttestPerform) {
    ttestName   <- "tTests"
    ttest.path <- ddCt:::getDir(file.path(params$savePath, ttestName))
    
    if (!is.null(params$samplesRemainInTTest)){
      ttestRed <- (rownames(pData(result)) %in% params$samplesRemainInTTest)
    }else{
      ttestRed <- !(rownames(pData(result)) %in% params$samplesNotInTTest)
    }
  
    result3 <- result[,ttestRed]
  
    daten  <- dCt(result3) ## TTest always with normal values
    if( ! params$ttestCol %in% colnames(pData(result3)) )
      stop(paste(" did not find :", params$ttestCol,": in pData ",sep=""))
    faktor <- as.character(pData(result3)[,params$ttestCol])
    mmm <- nlevels(as.factor(faktor))
    if( mmm == 1 ){
      stop( " found only a single group for t-test ")
    }else if( mmm == 2 ){
      aa <- matrix(c(1,2), ncol=1) ## aa = comparing by pair
    }else{
      aa <- combn(1:mmm,2)
    }
    
    
    for(k in 1:ncol(aa)) {
      Groups  <- levels(as.factor(faktor))[aa[,k]]
      subs    <- faktor %in%  Groups
      datenS  <- daten[,subs]
      faktorS <- as.factor(faktor[subs])  ## force factor with only two elements
      
      if( ! is.null(params$pairsCol) ) {
        if( ! params$pairsCol %in% colnames(pData(result)) )
          stop(paste(" did not find :", params$pairsCol,": in pData ",sep=""))
        paarungS <- as.character(pData(result3)[,params$pairsCol])
        paarungS <- paarungS[subs]
        wenigerRes <- 1
        optTest <- "paired "
      } else {
        wenigerRes <- 0
        optTest <- ""
      }
      
      res  <- matrix(NA,nrow=nrow(datenS),ncol=8-wenigerRes)
      res2 <- matrix(NA,nrow=nrow(datenS),ncol=2-wenigerRes)
      
      for (i in 1:nrow(datenS)){ 
        a  <- datenS[i,]
        b  <- is.na(a)
        cc <- a[!b]
        d  <- faktorS[!b]
        if( ! is.null(params$pairsCol) ){ 
          paar <- as.character(paarungS[!b])
          ## restrict to valid pair data only
          validPaarItems <- paar[duplicated(paar)]
          valid <- which( as.character(paar) %in% validPaarItems )
          paar <- paar[valid]
          cc <- cc[valid]
          d <- d[valid]
          if(all(table(d) >1)) {
            sel1 <- which(as.character(d) == Groups[1])
            sel2 <- which(as.character(d) == Groups[2])
            stopifnot( all( paar[sel1][order(paar[sel1])] ==paar[sel2][order(paar[sel2])] )   )
            group1 <- cc[sel1][order(paar[sel1])]
            group2 <- cc[sel2][order(paar[sel2])]
            ff <- t.test(x=group1, y=group2, paired=TRUE)
            ff2<- wilcox.test(x=group1, y=group2, paired=TRUE)
          }else{
            ff <- NULL
            ff2 <- NULL
          }
        }else{
          if(all(table(d)>1)) {
            ff <- t.test(cc~d)
            ff2<- wilcox.test(cc~d)
          }else{
            ff  <- NULL
            ff2 <- NULL
          }
        }
        if( !is.null(ff) ){
          res[i,]  <-c(signif(ff$statistic),
                       signif(ff$p.value),
                       ff2$statistic,
                       signif(ff2$p.value),
                       myFoldChange(ff$estimate),
                       ff$parameter["df"],
                       ff$estimate)          ## 1 or 2 Units (paired)
        
          res2[i,] <-c(names(ff$estimate))  ## 1 or 1 Units (paired)
        }
      }
      AllGe      <- rownames(level(result))
      theHKGenes <- rep("",length(AllGe))
      theHKGenes[AllGe %in% params$referenceGene] <- "X"
      gg <- cbind(AllGe,theHKGenes,res)
      
      myColnames <- c("Name",
                      "Housekeeping Gene",
                      paste("statistic(", optTest,"t.test)", sep=""),
                      paste("pvalue(",    optTest,"t.test)", sep=""),
                      paste("statistic(", optTest,"Wilcox)", sep=""),
                      paste("pvalue(",    optTest,"Wilcox)", sep=""),
                      "foldChange",
                      "degreeOfFreedom"
                      )
      
      if( ! is.null(params$pairsCol) ){
        Mr1 <- res2[,1]
        extraName <- unique(Mr1[!is.na(Mr1)])
        if( length(extraName) < 1 ){
          extraName <- "fehlenUnklar"
        }
        myColnames <- c(myColnames, extraName) 
      }else{
        Mr1 <- res2[,1]
        Mr2 <- res2[,2]
        stopifnot(length(unique(Mr1[!is.na(Mr1)]))==1 & length(unique(Mr2[!is.na(Mr2)]))==1   )
        myColnames <- c(myColnames,
                        unique(Mr1[!is.na(Mr1)]),
                        unique(Mr2[!is.na(Mr2)])
                        ) 
      }
      
      colnames(gg) <- myColnames
      pVSpalte <- paste("pvalue(",    optTest,"t.test)", sep="")
      gg <-gg[order(as.numeric(gg[,pVSpalte])),]
      
      SAVED <- paste("ttest",Groups[1],Groups[2],sep="")
      write.table(gg,file=file.path(ttest.path, paste(SAVED,"txt",sep=".")),sep="\t",row.names=FALSE)
      write.htmltable(gg,file=file.path(ttest.path,SAVED))
    }
  }
}

##----------------------------------------##
## doBoxplotsHousekkepingGenes
##----------------------------------------##
doBoxplotsHousekkepingGenes <- function(params, result, savedir) {
  if(params$ttestPerform){
    if (!is.null(params$samplesRemainInTTest)) {
      ttestRed <- (rownames(pData(result)) %in% params$samplesRemainInTTest)
    } else {
      ttestRed <- !(rownames(pData(result)) %in% params$samplesNotInTTest)
    }
    result3 <- result[,ttestRed]
    daten   <- Ct(result3)  ## its about the expression of 
    ## Housekeeping Gene before normalisation
    faktor	<- as.character(pData(result3)[,params$ttestCol])
    mmm     <- levels(as.factor((faktor)))
    N <- length(mmm)
    
    daten2 <- daten[params$referenceGene,,drop=FALSE] 
    BoxPl <- list()
    for (i in 1:N) {
      BoxPl[[i]] <- t(daten2[,mmm[i]==faktor])
    }
    
    res <- list()
    for(i in 1:length(params$referenceGene)){
      A <- lapply(BoxPl, function(x) x[,i])
      names(A) <- rep(params$referenceGene[i], N)
      res      <- c(res,A)
    }
    
    theColor <- 2 + (1:N)
    pdf(file=file.path(savedir, paste("HKGenesPerGroup",".pdf",sep="")),w=15,h=15)
    boxplot(res,col=theColor,main="Ct expression of housekeeping genes per group")
    dev.off()
  }
}

##----------------------------------------##
## ddCtExec - executes ddCt with the given
##             arguments
##----------------------------------------##
ddCtExec <- function(params.new, sdmframe) {
  params <- ddCt:::getDefaultParams(param.list)
  for(param.name in names(params)) {  
    ## overwrite with config file values
    if(!is.null(params.new[[param.name]]))
      params[[param.name]] <- params.new[[param.name]]
  }

  if(!missing(sdmframe))
    params$inputFile = fileNames(sdmframe)

  ddCt:::checkParams(params)
  
  if(is.null(params$params$samplesRemainInTTest))
    params$params$samplesRemainInTTest = params$samplesRemainInOutput
  if(is.null(params$samplesNotInTTest))
    params$samplesNotInTTest = params$samplesNotInOutput
  if(is.null(params$samplesRemainInCor))
    params$samplesRemainInCor = params$samplesRemainInOutput
  if(is.null(params$samplesNotInCor))
    params$samplesNotInCor = params$samplesNotInOutput

  
  sho <- paste(gsub("^([^.]+)?.*$", "\\1", basename(params$inputFile)), collapse="_")
  
  savedir <- ddCt:::getDir(file.path(params$savePath, paste("Result", sho, sep="_")))
  datadir <- params$loadPath

  if(missing(sdmframe))
    sdmframe <- SDMFrame(file.path(datadir,params$inputFile))
  
  warningFile <- file.path(savedir, paste("warning.output",".txt",sep=""))

#########################
## reading & calculation
#########################
  
  input <- prepareInput(params, sdmframe, datadir)
  result <- ddCtExpression(input$CtData,
                           algorithm=params$algorithm,
                           calibrationSample=params$referenceSample,
                           housekeepingGenes=params$referenceGene,  
                           type=params$mode,
                           sampleInformation=input$sampleInformation,
                           efficiencies=params$efficiencies,
                           efficiencies.error=params$efficienciesError)

###################
## output
###################
  doTables(params, result, savedir)
  doPlot(params, result, savedir)
  doCorrelation(params, result, savedir)
  doBoxplotsHousekkepingGenes(params, result, savedir)
  
###################
## error-handling
###################
  ## If warnings have been created during the calculation process they appear here:
  if(file.exists(warningFile)) {
    bart  <- unlist(read.delim(warningFile,as.is=TRUE,header=FALSE))
    error <- sapply(bart, function(y) gsub("simpleWarning in withCallingHandlers\\(\\{: ","",y))
    error
  }  
}

################################################################################
### Main
################################################################################

## auto execude ddCt.R if script is loaded with Rscript
if(ddCt:::getExecMode() == "script") {
  params <- ddCt:::getParams(param.list)
  ddCtExec(params)
}
