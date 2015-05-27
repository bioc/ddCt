################################################################################
##
## This software is created by Molecular Genom Analysis Group
## Department of German Cancer Research Center in Heidelberg
##
##
##  AllMethods.R
##  Created on: Oct 23, 2008
##      Author: Rudolf Biczok <r.biczok@dkfz-heidelberg.de>
##      Description: ddCt classes
##
################################################################################

################################################################################
## Class ddCtExpression
################################################################################

##----------------------------------------##
## errBarchart
##----------------------------------------##
setMethod("errBarchart", c("ddCtExpression", "missing"),
          definition = function(object, by,...) {
            res <- elist(object)
            ddCtErrBarchart(res,by="Sample",...)
          })

setMethod("errBarchart", c("ddCtExpression", "character"),
          definition = function(object, by, ...) {
            res <- elist(object)
            ddCtErrBarchart(res,by=by, ...)
          })

##----------------------------------------##
## elist, elistWrite, as(from ,"data.frame")
##----------------------------------------##
setAs(from="data.frame",to="ddCtExpression", 
      def=function(from) {
        to@coreData <- from
      })

##----------------------------------------##
## elist, elistWrite, as(from ,"data.frame")
##----------------------------------------##
setAs(from="ddCtExpression",to="data.frame", 
      def=function(from) {
        aD <- assayData(from)
        if(is.data.frame(aD)) {
          return(aD)
        } else {
          stopifnot(is.list(aD))
          a <- data.frame(expand.grid(Detector=featureNames(from), Sample=sampleNames(from)))
          b <- colnames(a)
          for (i in 1:length(aD)) {
            a <- cbind(a, as.vector(aD[[i]]))
          }
          colnames(a) <- c(b, names(aD))
          return(a)
        }
      })

setMethod("elist", "ddCtExpression", 
          definition = function(object,colnames=NULL) {
            df <- as(object, "data.frame")
            if(!is.null(colnames)) {
              colnames(df) <- as.character(colnames)
            }
            return(df)
          })

setMethod("summary", "ddCtExpression",
          function(object, colnames=NULL) {
            elist(object, colnames)
          })

setMethod("elistWrite",signature("ddCtExpression","character"), 
          definition = function(object,file, sep="\t", quote=FALSE, row.names=FALSE,colnames=NULL,...) {
            a <- elist(object, colnames)
            write.table(a,file=file,sep=sep, quote=quote, row.names=row.names,...)
          })

##----------------------------------------##
## getter methods
##----------------------------------------##
setMethod("level", "ddCtExpression", function(object) {
  return(exprs(object))
})

setMethod("levelErr", "ddCtExpression", function(object) {
  return(assayData(object)$level.err)
})

setMethod("Ct", "ddCtExpression", function(object) {
  return(assayData(object)$Ct)
})

setMethod("CtErr", "ddCtExpression", function(object) {
  return(assayData(object)$Ct.error)
})

setMethod("dCt", "ddCtExpression", function(object) {
  return(assayData(object)$dCt)
})

setMethod("dCtErr", "ddCtExpression", function(object) {
  return(assayData(object)$dCt.error)
})

setMethod("ddCt", "ddCtExpression", function(object) {
  return(assayData(object)$ddCt)
})

setMethod("ddCtErr", "ddCtExpression", function(object) {
  return(assayData(object)$ddCt.error)
})

setMethod("numberCt", "ddCtExpression", function(object) {
  return(assayData(object)$number)
})

setMethod("numberNA", "ddCtExpression", function(object) {
  return(assayData(object)$numberNA)
})

################################################################################
## Class InputFrame
################################################################################

##----------------------------------------##
## constructor 
##----------------------------------------##

setMethod("initialize", "InputFrame",
          function(.Object, cd, file.name) {
            if(!missing(cd)) {
              if(!missing(cd))
                .Object@files  <- file.name
              if (! all (PRIMARY.INPUT.COLS %in% colnames(cd)))
                stop(gettextf("Your file does not contain the columns %s and '%s'.",
                              paste(sQuote(PRIMARY.INPUT.COLS[1:length(PRIMARY.INPUT.COLS)-1]), collapse=" "),
                              PRIMARY.INPUT.COLS[length(PRIMARY.INPUT.COLS)]))
              cd <- cd[,PRIMARY.INPUT.COLS]
              Platename <- rep(basename(file.name), nrow(cd))
              coreData(.Object) <- rbind(coreData(.Object),cbind(cd,Platename))
              ow <- getOption("warn")
              options(warn=-1)
              coreData(.Object)[,"Ct"] <- as.numeric(coreData(.Object)[,"Ct"])
              options(warn=ow)
              coreData(.Object)[,"Platename"] <- as.character(coreData(.Object)[,"Platename"])
              
              isNotNilStr <- coreData(.Object)$Sample!= "" & coreData(.Object)$Detector != ""
              coreData(.Object) <- coreData(.Object)[isNotNilStr,]
            }
            return(.Object)
          })

##----------------------------------------##
## getter methods
##----------------------------------------##
setMethod("coreData", "InputFrame", function(object) {
  return(object@coreData)
})

setMethod("fileNames", "InputFrame", function(object) {
  return(object@files)
})

setMethod("names", "InputFrame", function(x) {
  return(names(x@coreData))
})

setMethod("detectorNames", "InputFrame", function(object) {
  return(object@coreData$Detector)
})

setMethod("sampleNames", "InputFrame", function(object) {
  return(object@coreData$Sample)
})

setMethod("uniqueDetectorNames", "InputFrame", function(object) {
  return(unique(object@coreData$Detector))
})

setMethod("uniqueSampleNames", "InputFrame", function(object) {
  return(unique(object@coreData$Sample))
})

setMethod("$", "InputFrame", function(x, name) {
  return(x@coreData[,name])
})

setMethod("[", "InputFrame", function(x, i, j, ..., drop) {
  subCoreData <- x@coreData[i, j, ..., drop=FALSE]
  coreData(x) <- subCoreData
  return(x)
})

##currently disabled
##setMethod("[[", "InputFrame", function(x, i,j, ...,exact) {
##  return(x@coreData[[i=i, j=j, ..., exact=exact]])
##})

setMethod("Ct", "InputFrame", function(object) {
  coreData(object)$Ct
})


## assign method
setReplaceMethod("Ct", c("InputFrame", "numeric"), function(object, value) {
  coreData(object)$Ct <- value
  return(object)
})

##----------------------------------------##
## casting methods
##----------------------------------------##

setAs(from="InputFrame", to="data.frame", def=function(from) {
  return(coreData(from))
})

### S3 methods

as.data.frame.InputFrame <- function(x, row.names=NULL, optional=FALSE,...) {
  cd <- coreData(x)
  return(as.data.frame(cd, row.names=row.names, optional=optional,...))
}

##----------------------------------------##
## setter methods
##----------------------------------------##

setReplaceMethod("detectorNames", signature(object="InputFrame",value="character"),
                 function(object, value) {
                   object@coreData$Detector <- value
                   return(object)
                 })

setReplaceMethod("coreData", signature(object="InputFrame", value="data.frame"), function(object, value) {
  object@coreData <- value
  return(object)
})

setReplaceMethod("sampleNames", signature(object="InputFrame", value="character"),
                 function(object, value) {
                   object@coreData$Sample <- value
                   return(object)
                 })

setReplaceMethod("uniqueDetectorNames", signature(object="InputFrame",target="character", value="character"),
                 function(object, target, value) {
                   object@coreData$Detector <- replaceNames(object@coreData$Detector, target, value)
				   return(object)
                 })

setReplaceMethod("uniqueDetectorNames", signature(object="InputFrame",target="missing", value="character"),
                 function(object, target, value) {
                   target <- names(value)
                   value <- as.character(value)
                   object@coreData$Detector <- replaceNames(object@coreData$Detector, target, value)
				   return(object)
                 })

setReplaceMethod("uniqueSampleNames", signature(object="InputFrame", target="character", value="character"),
                 function(object, target, value) {
		   object@coreData$Sample <- replaceNames(object@coreData$Sample, target, value)
                   return(object)
                 })

setReplaceMethod("uniqueSampleNames", signature(object="InputFrame",target="missing", value="character"),
                 function(object, target, value) {
                   target <- names(value)
                   value <- as.character(value)
		   object@coreData$Sample <- replaceNames(object@coreData$Sample, target, value)
                   return(object)
                 })

setMethod("replaceDetector", signature(object="InputFrame", target="character", value="character"),
          function(object, target, value) {
            newDetectorNames <- replaceVectorByEquality(detectorNames(object), target, value)
            detectorNames(object) <- newDetectorNames
            return(object)
          })

setMethod("replaceSample", signature(object="InputFrame", target="character", value="character"),
          function(object, target, value) {
            newSampleNames <- replaceVectorByEquality(sampleNames(object), target, value)
            sampleNames(object) <- newSampleNames
            return(object)
          })
setMethod("removeSample", signature(object="InputFrame", sample="character"),
          function(object, sample) {
            object <- object[!sampleNames(object) %in% sample,]
            return(object)
          })
setMethod("removeDetector", signature(object="InputFrame", detector="character"),
          function(object, detector) {
            object <- object[!detectorNames(object) %in% detector,]
            return(object)
          })
setMethod("removeNTC", "InputFrame", function(object) {
  is.ntc <- grepl("^ntc",
                  as.character(sampleNames(object)),
                  ignore.case=TRUE)
  return(object[!is.ntc,])
})
setMethod("removeNTC", "ddCtExpression", function(object) {
  is.ntc <- grepl("^ntc",
                  as.character(sampleNames(object)),
                  ignore.case=TRUE)
  return(object[,!is.ntc])
})

##----------------------------------------##
## censoring
##----------------------------------------##
setMethod("rightCensoring", signature(object="InputFrame", threshold="numeric"),
          function(object, threshold, value) {
            if(missing(value))
              value <- as.numeric(NA)
            cData <- coreData(object)
            cData$Ct[cData$Ct > threshold] <- value
            coreData(object) <- cData
            return(object)
          })
                                      

##----------------------------------------##
## ddCtExpression public user method
##----------------------------------------##

setMethod("ddCtExpression", "InputFrame",
          function(object,
                   warningStream,
                   algorithm,
                   calibrationSample,
                   housekeepingGenes,
                   type,
                   sampleInformation,
                   toZero,
                   efficiencies,
                   efficiencies.error) {

            if(missing(algorithm))
              algorithm="ddCt"
            if(missing(warningStream))
              warningStream = "warnings.txt"
            if(missing(type))
              type="mean"
            if(missing(toZero))
              toZero = FALSE

            warningHandler <- function(x) {
              ww <- file(warningStream, open="a+")
              writeLines(as.character(x),con=ww)
              close(ww)
            }
            
            if(algorithm == "ddCt")
              return(withCallingHandlers(ddCtExec(object,
                                                  calibrationSample=calibrationSample,
                                                  housekeepingGenes=housekeepingGenes,
                                                  type=type,
                                                  sampleInformation=sampleInformation,
                                                  toZero=toZero),
                                         warning = warningHandler))
            else if(algorithm == "ddCtWithE")
              return(withCallingHandlers(ddCtWithEExec(object,
                                                       calibrationSample=calibrationSample,
                                                       housekeepingGenes=housekeepingGenes,
                                                       type=type,
                                                       sampleInformation=sampleInformation,
                                                       efficiencies=efficiencies,
                                                       efficiencies.error=efficiencies.error),
                                         warning = warningHandler))
            else
              stop(gettextf("'%s' is an invalid algorithm"))
          })

##----------------------------------------##
## ddCtwithEExec ddCt algorithm
##               with efficiencies
##----------------------------------------##
setMethod("ddCtWithEExec", "InputFrame",
          function(object,
                   calibrationSample,
                   housekeepingGenes,
                   type,
                   sampleInformation,
                   efficiencies,
                   efficiencies.error) {
  
            aaa <- object$Ct
            bbb <- object$Platename
            reduced.set <- object[,c("Sample","Detector")]

            if (!all(housekeepingGenes %in% reduced.set[,2]))
              stop("Not all of your housekeeping genes are in your table", call.=FALSE)
            if (! all(calibrationSample %in% reduced.set[,1]))
              stop("At least one of your reference samples is not in your table.", call.=FALSE)
            if (! type %in% c("median","mean"))
              stop("Type must be median or mean!", call.=FALSE)
            
            the.difference <- function(x) {if (any(is.na(x))|length(x) < 2 ) y <- NA else y <- max(diff(sort(x))); return (y)}
            sum.na         <- function(x) {sum(is.na(x))}
            unique.plate   <- function(x) {
              if (length(unique(x))!=1) warning(paste("g-s comb. on more than one plate:",paste(unique(x),collapse=",")))
              return (unique(x)[1])}
            
            number.of.na           <- tapply(aaa,reduced.set,sum.na)                  # number of points with NA
            number.of.all          <- tapply(aaa,reduced.set,length)                  # number of replication
            
            if (type=="median"){
              the.Ct.values         <- tapply(aaa,reduced.set,na.median,na.rm=TRUE)       # Median
              error.Ct.mad          <- tapply(aaa,reduced.set,na.mad,na.rm=TRUE,con=1)    # MAD 
              error.Ct              <- error.Ct.mad/sqrt(number.of.all - number.of.na )
            } else {
              the.Ct.values         <- tapply(aaa,reduced.set,na.mean,na.rm=TRUE)         # Mean
              error.Ct.sd           <- tapply(aaa,reduced.set,na.sd,na.rm=TRUE)           # SD 
              error.Ct              <- error.Ct.sd/sqrt(number.of.all - number.of.na ) # SEM 
            }
            
            the.difference.values  <- tapply(aaa,reduced.set,the.difference)          # ratio long distance short distance
            the.plate              <- tapply(bbb,reduced.set,unique.plate)
  
  
            ## warning messages if a reference sample or a housekeeping gene has no values
            
            for (sample in calibrationSample)
              for( Detector in unique(reduced.set[,2])) {
                if (is.na(the.Ct.values[sample,Detector]))
                  warning(paste("No value for gene",Detector,"in ref. sample",sample), call.=FALSE)
              }
  
            for (Detector in housekeepingGenes)
              for( sample in unique(reduced.set[,1])){
                if (is.na(the.Ct.values[sample,Detector]))
                  warning(paste("No value for housekeeping gene",Detector,"in sample",sample), call.=FALSE)
              }
            
            
            GeneNames   <- colnames(the.Ct.values)
            SampleNames <- rownames(the.Ct.values)
            NCOL <- length(GeneNames)
            NROW <- length(SampleNames)
            HKG <- which(GeneNames %in% housekeepingGenes)
            RS  <- which(SampleNames %in% calibrationSample)
            
            if (missing(efficiencies)){
              efficiencies <- rep(2,NCOL)
            } else {
              efficiencies <- efficiencies[GeneNames]
            }
            if (missing(efficiencies.error)) {
              efficiencies.error <- rep(0,NCOL)
            } else {
              efficiencies.error <- efficiencies.error[GeneNames]
            }
            
### functions for symbolic representation
            
            symb.mean <- function(x){
              x <- x[!is.na(x)]
              y <- paste(x,collapse="+")
              g <- length(x)
              return(paste("1/",g,"*(",y,")",sep=""))
            }
            
            
            symb.exp <- function(x,a){
              return(paste(a,"^",x,sep=""))
            }
            
            make.start.matrix<- function(xd,yd){
              f <- xd*yd
              The.Anfang <- matrix(paste("Y",1:f,sep=""),ncol=yd)
              res <-sapply(The.Anfang[,1],symb.exp,a=paste("X",1,sep=""))
              for(i in 2:yd){
                res <- cbind(res,sapply(The.Anfang[,i],symb.exp,a=paste("X",i,sep="")))
              }
              return(res)
            }
            
            A <- make.start.matrix(NROW,NCOL)

            ##---------------------------##
            ## The assignmen
            ##---------------------------##
            
            hhh <- as.vector(the.Ct.values)
            
            for ( i in 1:length(hhh))
              assign(paste("Y",i,sep=""),hhh[i])
            
            for(i in 1:length(efficiencies))
              assign(paste("X",i,sep=""),efficiencies[i])
            
            THEERRORS <- c(efficiencies.error,as.vector(error.Ct))
            
            ##---------------------------##
            ## Calculation 
            ##---------------------------##
  
            result           <- matrix(NA, ncol=NCOL,nrow=NROW)
            rownames(result) <- SampleNames
            colnames(result) <- GeneNames
            result.error     <- result
            
            for(i in 1:NROW)
              for(j in 1:NCOL){
                P1 <- paste(symb.mean(A[i,HKG]),"/",A[i,j],sep="")
                dummy <- c()
                for(k in RS)
                  dummy <- c(dummy,paste(symb.mean(A[k,HKG]),"/",A[k,j],sep=""))
                
                P2 <- paste("1/(",symb.mean(dummy),")",sep="")
                B <- eval(deriv(as.formula(paste("~",P1,"*",P2,sep="")),c(paste("X",1:NCOL,sep=""),paste("Y",1:(NCOL*NROW),sep=""))))
                result      [i,j] <- B[1]
                result.error[i,j] <- sqrt(sum(attr(B,"gradient")^2 * THEERRORS^2,na.rm=TRUE))
              }
            
            ##---------------------------##
            ## putting the stuff together 
            ##---------------------------##

            cali        <- rownames(the.Ct.values) %in% calibrationSample
            Samplenames <- rownames(the.Ct.values) 
            names(cali) <- Samplenames 
            
            ## create the annotated data frame
  
            ## create the annotated data frame
            a  <- new("AnnotatedDataFrame", data=data.frame(Sample=Samplenames, Calibrator=cali),
                      varMetadata=data.frame(labelDescription=c("given by user",
                                               "Is this sample used as reference?"),
                        row.names=c("Sample", "Calibrator")))
            
            
 #DF <-data.frame("labelDescription"=c("given by user","Is this sample used as reference?"))
 #rownames(DF) <-c("Sample","Calibrator")
 #a <- new("AnnotatedDataFrame", data=data.frame(Sample=Samplenames,Calibrator=cali), varMetadata=DF)

  
            theList <- list(exprs = t(result),
                            level.err  = t(result.error),
                            Ct         = t(the.Ct.values),
                            Ct.error   = t(error.Ct),
                            efficiencies =  matrix(rep(efficiencies,NROW), ncol=NROW),
                            efficiencies.error =  matrix(rep(efficiencies.error,NROW),ncol=NROW),
                            Difference = t(the.difference.values),
                            numberNA   = t(number.of.na),
                            number     = t(number.of.all),
                            Plate      = t(the.plate))
    
# result <- new("MultiSet",assayData=theList,
#                      phenoData  = a)
#                      sampleNames= rownames(pData(a)),
#                      reporterNames =colnames(the.level))
            
            result <- new("ddCtExpression",assayData=theList, phenoData=a,sampleNames= rownames(pData(a)))
  
  
            if (!missing(sampleInformation)) {
              if( !("Sample" %in% colnames(pData(sampleInformation))))
                stop("Your phenoData must contain a column named 'Sample'.")
              the.match <- match(rownames(pData(result)),as.character(pData(sampleInformation)$Sample))
              thePheno <- phenoData(result)
              pData(thePheno) <- cbind(pData(thePheno),pData(sampleInformation)[the.match,colnames(pData(sampleInformation))!="Sample",drop=FALSE])
              newVar <- as.data.frame(as.matrix(varLabels(sampleInformation)[names(varLabels(sampleInformation))!="Sample"]))
              colnames(newVar) <- "labelDescription"
              varMetadata(thePheno) <- rbind(varMetadata(thePheno),newVar)
              phenoData(result) <- thePheno
   #pData(result) <- cbind(pData(result),pData(sampleInformation)[the.match,colnames(pData(sampleInformation))!="Sample"])
   #phenoData(result)@varLabels<- c(varLabels(phenoData(result)),varLabels(sampleInformation)[names(varLabels(sampleInformation))!="Sample"])
   #colnames(pData(result)) <- names(phenoData(result)@varLabels)
            }
            return(result)
          })

##----------------------------------------##
## ddCtExec ddCt method calculates relative
##          expression with ddCt method
##----------------------------------------##

setMethod("ddCtExec", "InputFrame",
          function(object,
                   calibrationSample,
                   housekeepingGenes,
                   type,
                   sampleInformation,
                   toZero) {
 
            cts <- object$Ct
            platenames <- object$Platename

            ##--------------------##
            ## checking
            ##--------------------##
            type <- match.arg(type, c("median","mean"), several.ok=FALSE)
            
            ## unique detector and sample names
            uDetNames <- uniqueDetectorNames(object)
            uSamNames <- uniqueSampleNames(object)
            dsPair <- coreData(object[,c("Sample","Detector")])

            isHKGincluded <- all(housekeepingGenes %in% uDetNames)
            isCALincluded <- all(calibrationSample %in% uSamNames)
            if (!isHKGincluded) stop("Not all your housekeeping genes are in your table", call.=FALSE)
            if (!isCALincluded) stop("Not all reference samples are not in your table.", call.=FALSE)
            
            ## NA data points
            numberNA  <- tapply(cts,dsPair,function(x) sum(is.na(x)))
            ## All data points
            numberAll <- tapply(cts,dsPair,length)
            ## Effective data points
            numberEff <- numberAll - numberNA
            
            if (type=="median"){
              Ct         <- tapply(cts,dsPair,na.median,na.rm=TRUE)       # Median
              ## ?? question: should the constant here set to 1?
              ctMad          <- tapply(cts,dsPair,na.mad,na.rm=TRUE,con=1)    # MAD 
              CtError              <- ctMad/sqrt(numberEff)
            } else{
              Ct         <- tapply(cts,dsPair,na.mean,na.rm=TRUE)         # Mean
              CtSd           <- tapply(cts,dsPair,na.sd,na.rm=TRUE)
              if(toZero) CtSd[numberEff==1] <- 0             # SD 
              CtError              <- CtSd/sqrt(numberEff) # SEM 
            }
            difference  <- tapply(cts,dsPair,getDiff)          # ratio long distance short distance
            plate <- tapply(platenames,dsPair,uniquePlate)
            
            
            ## warning messages if a reference sample or a housekeeping gene has no values
            for (sample in calibrationSample)
              for( Detector in uDetNames){
                if (is.na(Ct[sample,Detector]))
                  warning(paste("No value for gene",Detector,"in ref. sample",sample), call.=FALSE)
              }
            
            for (Detector in housekeepingGenes)
              for( sample in uSamNames){
                if (is.na(Ct[sample,Detector]))
                  warning(paste("No value for housekeeping gene",Detector,"in sample",sample), call.=FALSE)
              }
            
# if (any (is.na( Ct.of.reference.gene))){
#   b <- names(Ct.of.reference.gene)[is.na(Ct.of.reference.gene)]
#   warning(paste("There is/are no Ct values of the reference gene for the following sample/s:",paste(b,collapse=",")))}
    
            ##------------------------------------------------------------##
            ## Ct and error calculation for the housekeeping gene(s) 
            ##------------------------------------------------------------##
            
            refGeneCt <- rowMeans(Ct[,housekeepingGenes,drop=FALSE])
            
            refGeneError  <- CtError[,housekeepingGenes,drop=FALSE]  
            refGeneLen <- length(unique(housekeepingGenes))
            refGeneCtError <- 1/refGeneLen * sqrt(rowSums((refGeneError)^2))
            
            
            ##------------------------------------------------------------##
            ## dCt and errors
            ##------------------------------------------------------------##
            
            dCt         <- Ct - refGeneCt
            hkgMat   <- matrix(1, ncol=ncol(Ct),nrow=nrow(Ct))
            hkgMat[,colnames(Ct) %in% housekeepingGenes] <- 1 - 2/refGeneLen
            red.error   <- CtError^2 *hkgMat
            dummy       <- red.error + (refGeneCtError)^2
            dummy[hkgMat==-1] <- 0
            dCtError   <- sqrt(dummy)
            
            ##------------------------------------------------------------##
            ## dCt and error calculation for the ref samples
            ##------------------------------------------------------------##
            
            calSampDCt  <- colMeans(dCt[calibrationSample,,drop=FALSE],na.rm=TRUE)
            isNA <- which((is.na(dCt[calibrationSample,,drop=FALSE])))
            
            ## the delta Ct error for the reference samples, NA values are not taken into account (set them zero)
            calSampError <- dCtError[calibrationSample,,drop=FALSE]
            calSampError[isNA] <- 0
            
            ## for every gene the number of reference samples without NA values
            CS <- length(calibrationSample) -apply(dCt[calibrationSample,,drop=FALSE],2, function(x) sum(is.na(x)))
            
            ## for every gene the error of the mean of the reference samples
            calSampDCtError <- 1/CS * sqrt(colSums((calSampError)^2,na.rm=TRUE))
            
            ##------------------------------------------------------------##
            ## ddCt value and errors
            ##------------------------------------------------------------##
            
            ddCt <- t (t(dCt)-calSampDCt)
            ctMat <- matrix(1, ncol=ncol(Ct),nrow=nrow(Ct))
            dummy1 <- ctMat[rownames(Ct) %in% calibrationSample,,drop=FALSE] 
            ctMat[rownames(Ct) %in% calibrationSample,] <- matrix((1-2/CS),nrow=nrow(dummy1),ncol=ncol(dummy1), byrow=TRUE)
            red.error2  <- dCtError^2 * ctMat
            dummy2 <- t(red.error2) + calSampDCtError^2
            dummy2[t(ctMat==-1)] <- 0
            ddCtError <- t(sqrt(dummy2))
            

            ##------------------------------------------------------------##
            ## level values and error 
            ##------------------------------------------------------------##
            
            level       <- apply(ddCt,c(1,2),levelfkt)
            levelError <- log(2) * level *  ddCtError
            
            ##------------------------------------------------------------##
            ## build ddCtExpression object
            ##------------------------------------------------------------##
            
            cali <- rownames(ddCt) %in% calibrationSample
            sampleNames <- rownames(ddCt) 
            names(cali) <-sampleNames 
            
            ## create the annotated data frame
            a  <- new("AnnotatedDataFrame", data=data.frame(Sample=sampleNames, Calibrator=cali),
                      varMetadata=data.frame(labelDescription=c("given by user",
                                               "Is this sample used as reference?"),
                        row.names=c("Sample", "Calibrator")))
            
            
            theList <- list(exprs = t(level),
                            level.err  = t(levelError),
                            Ct         = t(Ct),
                            Ct.error   = t(CtError),
                            dCt        = t(dCt),
                            dCt.error  = t(dCtError),
                            ddCt       = t(ddCt),
                            ddCt.error = t(ddCtError),
                            Difference = t(difference),
                            numberNA   = t(numberNA),
                            number     = t(numberAll),
                            Plate      = t(plate))
            
            result <- new("ddCtExpression",assayData=theList, phenoData=a,
                          sampleNames= rownames(pData(a)))
            
# mit Objekten die aus read.phenoData kommen scheint das nicht zu klappen
# also mach ich es auf die alte Art und Weise
#
# if (! is.null(sampleInformation)) {
#      the.match <- match(sampleNames(phenoData(result)), sampleNames(sampleInformation))
#      if(length(unique(the.match)) != ncol(result))
#        stop("Your sample annotation does not fit the samples", call.=FALSE)
#      phenoData(result) <- combine(phenoData(result), sampleInformation)
#    }

            if (!missing(sampleInformation)){
              if( !("Sample" %in% colnames(pData(sampleInformation)))) stop("Your phenoData must contain a column named 'Sample'.")
              the.match             <- match(rownames(pData(result)),as.character(pData(sampleInformation)$Sample))
              
              thePheno              <- phenoData(result)
              pData(thePheno)       <- cbind(pData(thePheno),pData(sampleInformation)[the.match,colnames(pData(sampleInformation))!="Sample",drop=FALSE])
              newVar                <- as.data.frame(as.matrix(varLabels(sampleInformation)[names(varLabels(sampleInformation))!="Sample"]))
              colnames(newVar)      <- "labelDescription"
              varMetadata(thePheno) <- rbind(varMetadata(thePheno),newVar)
              phenoData(result)     <- thePheno
   #phenoData(result)@varLabels<-c(varLabels(phenoData(result)),varLabels(sampleInformation)[names(varLabels(sampleInformation))!="Sample"])
   #colnames(pData(result)) <- names(phenoData(result)@varLabels)
            }
            
            return(result)
          })

##----------------------------------------##
## utility methods
##----------------------------------------##

setMethod("headtailPrint", "data.frame", function(object, head=2L, tail=2L, digits=NULL, quote=FALSE, right=TRUE, row.names=TRUE) {
  subhead <- head(object, head)
  subtail <- tail(object, tail)
  subomit <- object[1,]; subomit <- rep("...", ncol(object))

  subd <- rbind(subhead,subomit, subtail)
  print(subd, digits=digits, quote=quote, right=right, row.names=row.names)
})

setMethod("show", "InputFrame", function(object) {
  cat(gettextf("An instance of %s:\n", class(object)))
  cored <- coreData(object)
  wid <- options()$width; decwid <- round(wid*0.85) ## decorator width
  decs <- paste(rep("=", decwid),collapse="")

  if(nrow(coreData(object)) > 0) {
    cat(decs, "\n")
    headtailPrint(cored, row.names=FALSE)
    cat(decs, "\n")
  }
  
  cat(nrow(cored), "Sample-Detector pairs\n")

  ## files
  nfiles <- length(object@files)
  cat(sprintf(ngettext(nfiles,
                       "File: %s", "Files: %s"),
              paste(sQuote(object@files), collapse=", ")), "\n")

})

setMethod("rbind2", signature(x="InputFrame", y="InputFrame"), function(x,y) {
  coreData(x) <- rbind(coreData(x),coreData(y))
  x@files     <- c(x@files, y@files) 
  return(x)
})

setMethod("rbind2", signature(x="InputFrame", y="data.frame"), function(x,y) {
  coreData(x) <- rbind(coreData(x),new("InputFrame",y))
  x@files     <- c(x@files, class(y)) 
  return(x)
})

################################################################################
## Class InputReader
################################################################################

##----------------------------------------##
## factory methods 
##----------------------------------------##

setMethod("InputFrame", "InputReader",
          function(object) {
            inputFrame <- new("InputFrame");
            for (file.name in object@files){
              rawData <- readRawData(object, file.name)
              names(rawData)[match(object@colmap@sample,colnames(rawData))] <- DEFAULT.SAMPLE.COLNAME
              names(rawData)[match(object@colmap@feature,colnames(rawData))] <- DEFAULT.FEATURE.COLNAME
              names(rawData)[match(object@colmap@ct,colnames(rawData))] <- DEFAULT.CT.COLNAME
              extra <- new("InputFrame", rawData, file.name)
              inputFrame <- rbind2(inputFrame, extra)
            }
            return(inputFrame)
          })

setAs(from="data.frame", to="InputFrame", def=function(from) {
  ## TODO: add check columns
  return(new("InputFrame",from, file.name=as.character(NA)))
})

setMethod("InputFrame", "data.frame", function(object) {
  as(object, "InputFrame")
})
################################################################################
## Class SDMReader
################################################################################

##----------------------------------------##
## constructor 
##----------------------------------------##

setMethod("initialize", "SDMReader",
          function(.Object, files) {
            .Object@files <- files
            .Object@colmap <- ColMap()
            return(.Object)
          })

setMethod("readRawData", signature(object="SDMReader",file.name="character"),
          function(object, file.name) {
            x <- scan(file=file.name,what="character",sep="\n",blank.lines.skip=TRUE,quiet=TRUE)
            CTstart <- grep("^Well",x)
            CTend   <- grep("^Summary",x)
            if(length(CTstart)==0 | length(CTend)==0) stop("Your file does not seem to be a .sdm file!")
            number.of.skips <- CTstart - 1
            number.of.rows  <- (CTend - 1 ) - CTstart ## the first one becomes the Header and will not affect nrow
            rawdata <- read.table(file.name,sep="\t",
                                  nrows=number.of.rows,
                                  skip=number.of.skips,
                                  header=TRUE,
                                  as.is=TRUE,
                                  blank.lines.skip = TRUE,
                                  comment.char="")
            return(rawdata)
          })

################################################################################
## Class TSVReader
################################################################################

##----------------------------------------##
## constructor 
##----------------------------------------##

setMethod("initialize", "TSVReader",
          function(.Object, files, colmap=ColMap()) {
            .Object@files <- files
            .Object@colmap <- colmap
            return(.Object)
          })

setMethod("readRawData", signature(object="TSVReader",file.name="character"),
          function(object, file.name) {
            rawdata <- read.table(file.name,sep="\t",
                                  header=TRUE,
                                  as.is=TRUE,
                                  blank.lines.skip = TRUE,
                                  comment.char="")
            return(rawdata)
          })

################################################################################
## Class QuantStudioReader
################################################################################

QUANTSTUDIO_COLMAP <- ColMap(sample="Sample Name",
                             feature="Target Name",
                             ct="Ct")

setMethod("initialize", "QuantStudioReader",
          function(.Object, files) {
            .Object@files <- files
            .Object@colmap <- QUANTSTUDIO_COLMAP
            return(.Object)
          })

setMethod("readRawData", signature(object="QuantStudioReader",file.name="character"),
          function(object, file.name) {
            x <- scan(file = file.name, what = "character", sep = "\n", 
                      blank.lines.skip = FALSE, quiet = TRUE)
            CTstart <- grep("^Well", x)
            isCT <- grepl("^[0-9]", x)
            if (length(CTstart) == 0 | sum(isCT) == 0)  {
              if(length(CTstart)==0)
                message("No 'well' line found")
              if(sum(isCT)==0)
                message("No valid wells found")
              stop("Your file", file, "does not seem to be a valid QuantStudio TaqMan file!")
            }
            number.of.skips <- CTstart - 1
            number.of.rows <- sum(isCT)
            rawdata <- read.table(file.name, sep = "\t", nrows = number.of.rows, quote="", 
                                  skip = number.of.skips, header = TRUE, as.is = TRUE, 
                                  blank.lines.skip = TRUE, comment.char = "", check.names=FALSE)
            return(rawdata)
          })

################################################################################
## Class ddCtParam
################################################################################

##----------------------------------------##
## constructor for parameter initialisation
##----------------------------------------##
setMethod("initialize", "ddCtParam", function(.Object, type, default) {
            .Object@type <- type
            .Object@default <- default
            
            return(.Object)
          })

##----------------------------------------------------------------------------##
## Class errBarchartParameter
##----------------------------------------------------------------------------##
setMethod("exprsUndeterminedLabel", "errBarchartParameter", function(object) {
  object@exprsUndeterminedLabel
})

setMethod("show", "errBarchartParameter", function(object) {
  cat("Label for Undetermined:", exprsUndeterminedLabel(object), "\n")
})


