################################################################################
##
## This software is created by Molecular Genom Analysis Group
## Department of German Cancer Research Center in Heidelberg
##
##
##  absoluteQ.R
##  Created on: Oct 23, 2008
##      Author: Rudolf Biczok <r.biczok@dkfz-heidelberg.de>
##      Description: TODO
##
################################################################################

ddCtAbsolute <- function(raw.table,
                 addData,
                 type="mean",
                 ADD = -30.234,
                 DIV = - 1.6268,
                 sampleInformation=NULL,
                 toZero=FALSE,
                 filename="warning.output.txt"){
  
 withCallingHandlers({
 #require(Biobase) 
 if (! all(c("Ct","Sample","Detector","Platename")%in% colnames(raw.table))) stop ("Your table must include columns with the following names : 'Ct','Sample','Detector','Platename'.")

 if (! all(c("Sample","Concentration","Number.reactions","per.reaction","overall")%in% colnames(addData))) stop ("Your table with additional informations must include columns with the following names : 'Sample','Concentration','Number.reactions','per.reaction','overall'.")

 
 aaa <- raw.table$Ct
 bbb <- raw.table$Platename
 reduced.set <- raw.table[,c("Sample","Detector")]
 
 if (! type %in% c("median","mean"))                 stop("Type must be median or mean!")

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
 } else{
  the.Ct.values         <- tapply(aaa,reduced.set,na.mean,na.rm=TRUE)         # Mean
  error.Ct.sd           <- tapply(aaa,reduced.set,na.sd,na.rm=TRUE)
  if(toZero) error.Ct.sd[number.of.all - number.of.na==1] <- 0             # SD 
  error.Ct              <- error.Ct.sd/sqrt(number.of.all - number.of.na ) # SEM 
 }
 the.difference.values  <- tapply(aaa,reduced.set,the.difference)          # ratio long distance short distance
 the.plate              <- tapply(bbb,reduced.set,unique.plate)

 



# if (any (is.na( Ct.of.reference.gene))){
#   b <- names(Ct.of.reference.gene)[is.na(Ct.of.reference.gene)]
#   warning(paste("There is/are no Ct values of the reference gene for the following sample/s:",paste(b,collapse=",")))}

#########################################################
# Ct and error calculation for the housekeeping gene(s) #
#########################################################
theTransformation <- function(x) exp((x + ADD)/DIV)

ddd <- apply(the.Ct.values, c(1,2), theTransformation)
k <- rownames(ddd)
 
AD  <- addData$Concentration * addData$Number.reactions * 5 * addData$per.reaction / addData$overall
AD <- AD[match(k,addData$Sample)]


 
 ## FIXME add warnings here
 
 dCt         <- ddd / AD
 error.dCt   <- abs(error.Ct * dCt * 1/DIV)


 
# addData
# Ct.of.reference.gene       <- rowMeans(the.Ct.values[,housekeepingGenes,drop=FALSE])
# all.housekeeping.error     <- error.Ct[,housekeepingGenes,drop=FALSE]  
# HKG                        <- length(housekeepingGenes)
# error.Ct.of.reference.gene <- 1/HKG * sqrt(rowSums((all.housekeeping.error)^2))


#################################
# the delta CT value and errors #
#################################
 
 #  dCt         <- the.Ct.values - Ct.of.reference.gene
 #  the.hkg.c   <- matrix(1, ncol=ncol(the.Ct.values),nrow=nrow(the.Ct.values))
 #  the.hkg.c[,colnames(the.Ct.values) %in% housekeepingGenes] <- 1 - 2/HKG
 #  red.error   <- error.Ct^2 *the.hkg.c
 #  dummy       <- red.error + (error.Ct.of.reference.gene)^2
 #  dummy[the.hkg.c==-1] <- 0
 #  error.dCt   <- sqrt(dummy)
 
#########################################################
# dCt and error calculation for the ref samples         #
#########################################################
 
# dCt.calibration.sample               <- colMeans(dCt[calibrationSample,,drop=FALSE],na.rm=TRUE)
## makingNAgoaway                       <- which((is.na(dCt[calibrationSample,,drop=FALSE])))
# all.ref.sample.error                 <- error.dCt[calibrationSample,,drop=FALSE]
# all.ref.sample.error[makingNAgoaway] <- 0
# CS                                   <- length(calibrationSample) -apply(dCt[calibrationSample,,drop=FALSE],2, function(x) sum(is.na(x)))
# error.dCt.calibration.sample         <- 1/CS * sqrt(colSums((all.ref.sample.error)^2,na.rm=TRUE))

#######################################
# the delta delta CT value and errors #
#######################################
 
# ddCt                 <- t (t(dCt)-dCt.calibration.sample)
# the.cs.c             <- matrix(1, ncol=ncol(the.Ct.values),nrow=nrow(the.Ct.values))
# the.cs.c[rownames(the.Ct.values) %in% calibrationSample,] <- 1 - 2/CS
# red.error2           <- error.dCt^2 * the.cs.c
# dummy2               <- t(red.error2) + error.dCt.calibration.sample^2
# dummy2[t(the.cs.c==-1)] <- 0
# error.ddCt           <- t(sqrt(dummy2))
                  

##########################
# level values and error #
##########################

# levelfkt <- function(x) 2^(-x)
#
# the.level       <- apply(ddCt,c(1,2),levelfkt)
# the.level.error <- log(2) * the.level *  error.ddCt

##############################
# putting the stuff together #
##############################

 #cali <- rownames(ddCt) %in% calibrationSample
 Samplenames <- rownames(dCt) 
 #names(cali) <-Samplenames 

 DF <-data.frame("labelDescription"=c("given by user"))
 rownames(DF) <-c("Sample")
 a <- new("AnnotatedDataFrame", data=data.frame(Sample=Samplenames), varMetadata=DF)
  
     
     
 result <- new("MultiSet",      assayData=list(exprs        = t(dCt),
                                           error  = t(error.dCt),
                                           Ct         = t(the.Ct.values),
                                           Ct.error   = t(error.Ct),
                                           Difference = t(the.difference.values),
                                           numberNA   = t(number.of.na),
                                           number     = t(number.of.all),
                                           Plate      = t(the.plate)),
                             phenoData  = a)
#                             reporterNames=colnames(dCt),
#                             sampleNames=rownames(dCt))
 
 if (! is.null(sampleInformation)) {
   if( !("Sample" %in% colnames(pData(sampleInformation)))) stop("Your phenoData must contain a column named 'Sample'.")
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


 return(result)},
  warning = function(x){ww <- file(filename,open="a+");writeLines(as.character(x),con=ww);close(ww)})
}
