################################################################################
##
## This software is created by Molecular Genom Analysis Group
## Department of German Cancer Research Center in Heidelberg
##
##
##  visualization.R
##  Created on: Oct 23, 2008
##      Author: Rudolf Biczok <r.biczok@dkfz-heidelberg.de>, Jitao David Zhang
##              <jitao_david.zhang@roche.com>
##      Description: visualization utility
##
################################################################################

##----------------------------------------##
## Visualization functions with lattice
##----------------------------------------##

## panel function for barchart with error
panel.barchart.errbar <- function(x,y, subscripts, groups, 
                                  width.errbar=0.08, lwd.errbar=2, col.errbar="black",thr, ...) {
  isLevel <- 1:(length(x)/2)
  ## draw bar chart plot, na is replaced by 0
  y.na.tol <- y[isLevel]; y.na.tol[is.na(y.na.tol)] <- 0; y.na.tol[y.na.tol>thr] <- thr+1
  panel.barchart(x[isLevel],y.na.tol[isLevel],...)

  ## draw error bars
  xerr <- match(x[isLevel + length(isLevel)], levels(x))
  levelerr <- y[isLevel + length(isLevel)]
  ci.l <- y[isLevel] - levelerr
  ci.u <- y[isLevel] + levelerr
  width.errbar=0.08; lwd.errbar=2
  for(i in 1:length(isLevel)) {
    x0 <- xerr[i]
    panel.segments(x0=x0, y0=ci.l[i], x1=x0, y1=ci.u[i], col=col.errbar, lwd=lwd.errbar)
    panel.segments(x0=x0-width.errbar, y0=ci.l[i], x1=x0+width.errbar, y1=ci.l[i], col=col.errbar, lwd=lwd.errbar)
    panel.segments(x0=x0-width.errbar, y0=ci.u[i], x1=x0+width.errbar, y1=ci.u[i], col=col.errbar, lwd=lwd.errbar)
  }
}

panel.ddCtErrBarchart <- function(x,y, thr,round, outText, parameter,
                                  exprs.one.line=TRUE,
                                  ...) {
  extremes <- which(y>thr);
  extremes <- extremes[extremes <= length(y)/2]
  nas <- which(is.na(y)); nas <- nas[nas <= length(y)/2]  
  panel.grid(v=0,h=-1)
  if(exprs.one.line) {
    panel.abline(h=1, lty="dashed", lwd=2, col="red")
  }
  panel.barchart.errbar(x,y,col.errbar="black",thr=thr,...)
  if(outText & length(extremes)>0) {
    for ( i in seq(along=extremes)) {
      extreme <-extremes[i]
      panel.text(x[extreme], thr,
                 substitute(x %+-% y,
                            list(x=round(y[extreme], round),
                                 y=round(y[extreme + length(y)/2], round))),
                 col="black",cex=.75, pos=1)
    }
  }
  panel.text(x[nas], 0, exprsUndeterminedLabel(parameter), pos=3,font=2, col="darkgrey")
}

##assignInNamespace( "ddCtErrBarchart", ddCtErrBarchart, "ddCt")

ddCtErrBarchart <- function(x, by=c("Sample", "Detector"),
                            thr=3,
                            ylab="Expression fold change",
                            cols=brewer.pal(12, "Set3"),round=0, outText=TRUE, rot=45,
                            parameter=new("errBarchartParameter"),
                            detector.levels=levels(factor(as.character(x$Detector))),
                            sample.levels=levels(factor(as.character(x$Sample))),
                            ...) {

  ## if all exprs is NA, the plot will not be interpretable
  if(all(is.na(x$exprs)))
    stop("All expressions are NA!\n")
  
  x$Sample <- factor(as.character(x$Sample), sample.levels)
  x$Detector <- factor(as.character(x$Detector), detector.levels)
  
  by <- match.arg(by, choices=c("Sample", "Detector"))
  if(by=="Sample") {
    formula <- as.formula("exprs + level.err ~ Sample | Detector")
  } else  {
    formula <- as.formula("exprs + level.err ~ Detector | Sample")
  }
  xlab <- by
  barchart(formula, data=x, 
           scales=list(x=list(rot=rot), y=list(alternating=1, at=seq(0, thr, 0.5))),
           ylim=c(0, thr*1.1),
           panel=function(x,y,...) panel.ddCtErrBarchart(x=x,y=y,thr=thr,round=round, outText=outText,parameter=parameter,...),
           xlab=xlab, ylab=ylab, col=cols,...  
           )
}


##----------------------------------------##
## old barploterrbar: to be defuncted
##----------------------------------------##
barploterrbar<- function (y, yl, yh, barcol = "orange", errcol = "black", horiz = FALSE, w = 0.2, theCut=NULL,columnForDiffBars=TRUE,cex.axis = par("cex.axis"),zeroForNA=TRUE,legend=FALSE,groups=NULL,order=FALSE, ...){
   ##.Deprecated("errBarchart")
   if (!(is.null(groups) | is.factor(groups))) stop("'groups' must be a factor.")
   if(order){
     grps <- split(1:length(y), groups)
     for(g in seq(along=grps))
       grps[[g]] <- grps[[g]][order(y[grps[[g]]])]
     y <- y[unlist(grps)]
     yl <- yl[unlist(grps)]
     yh <- yh[unlist(grps)]
   }
   if(is.matrix(y))
      newsp <- c(0.2,1) else{
      newsp <- 0.2
      y <- t(as.matrix(y))
      yl <- t(as.matrix(yl))
      yh <- t(as.matrix(yh))
    }

    if(zeroForNA){
      gg <- is.na(y)
      y[gg] <- 0
      yl[gg] <- 0
      yh[gg] <- 0
     }

    if(columnForDiffBars){
      y <- t(y)
      yl <- t(yl)
      yh <- t(yh)
     }

     co <- ncol(y)
     ro <- nrow(y)
     lwd <- 3 * 1/(ro*co)
     oldP <- par("cex.axis")
     par("cex.axis"=cex.axis)


    if (is.null(theCut)){
       if(all (!is.finite(yh)))
         {newR <- range(0,1)
      }else
         newR <-range(0,max(yh,na.rm=TRUE))
       
    }else
       newR  <- range(0,theCut)
       #newR <- range(0,min(max(EE[Gred,Sred],na.rm=TRUE),CUTOFF))
    

    
    if(legend) {layout(matrix(c(1,1,2,1,1,2), 2, 3, byrow=TRUE), respect=TRUE)
                 oldM <-par("mar")
                 newp <- oldM
                 newp[4] <- 0
                 par(mar=newp)}
    # here we assume that the user is quite clever and chooses a groups vector with the correct dimension
    if(!is.null(groups)){
        oldcol <- barcol[1:nlevels(groups)]
        barcol <- barcol[as.numeric(groups)]
      }
    if (!horiz) {
        x<- barplot(y, space = newsp, width = 1, col = barcol, horiz = horiz,beside=TRUE,ylim=newR,
                      cex.axis =cex.axis,...)
        for (i in 1:ro){
        segments(x[i,], yl[i,], x[i,], yh[i,], lwd = lwd, col = errcol)
        segments(x[i,] - w, yl[i,], x[i,] + w, yl[i,], lwd = lwd, col = errcol)
        segments(x[i,] - w, yh[i,], x[i,] + w, yh[i,], lwd = lwd, col = errcol)
        }
      } else {
       barplot(y, space = newsp, width = 1, col = barcol, horiz = horiz,beside=TRUE,xlim=newR,cex.axis =cex.axis,...)
       for (i in 1:ro){
       x <- (2+ (ro-1)*1.2) * seq(1, co)  - 0.5 - 1.2 * (ro - i)
       segments(yl[i,], x, yh[i,], x, lwd = lwd, col = errcol)
       segments(yl[i,], x - w, yl[i,], x + w, lwd = lwd, col = errcol)
       segments(yh[i,], x - w, yh[i,], x + w, lwd = lwd, col = errcol)
    }
     }
      if(legend){
        newp <- oldM
        newp[2] <- 0
        par(mar=newp)
        plot(0,0,main="", axes=FALSE, type="n",xlab="",ylab="")
          if (is.null(groups))
            legend(-1,1,rownames(y),col=barcol,fill=barcol)
          else
            legend(-1,1,levels(groups), col=oldcol,fill=oldcol)
        par("mar"=oldM)}                
   par("cex.axis"=oldP)
}
