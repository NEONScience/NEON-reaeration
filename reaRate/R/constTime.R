##############################################################################################
#' @title Constant rate injection half plateau specific conductance location

#' @author 
#' Kaelin M. Cawley \email{kcawley@battelleecology.org} \cr

#' @description This function determines the location of the half plateau specific conductance for use 
#' calculating travel time between an upstream and downstream sensor set.
#' @importFrom grDevices dev.new
#' @importFrom grDevices dev.off
#' @importFrom graphics identify
#' @importFrom graphics abline 
#' @importFrom graphics axis
#' @importFrom graphics lines
#' @importFrom graphics mtext
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics title
#' @importFrom stats lm
#' @importFrom stats loess.smooth
#' @importFrom stats lsfit

#' @param reaCond Specific conductance values of the stream during the experiment [uS/cm @ 25 C]
#' @param plotName Plot title to be displayed [string]
#' @param bPlot User input to display a visual check, default is FALSE [boolean]

#' @return This function returns the peak tracer location [integer]

#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

#' @keywords surface water, streams, rivers, velocity, travel time, reaeration, metabolism

#' @examples
#' #Using an example file
#' s1peakLoc <- constTime(reaCond = condDataS1, plotName = "Example_S1", bPlot = TRUE)

#' @seealso def.calc.travelTime for calculating travel times

#' @export

# changelog and author contributions / copyrights
#   Kaelin M. Cawley (2017-08-03)
#     original creation
##############################################################################################
constTime <- function(
  reaCond,
  plotName,
  bPlot
){
  ##### Constants #####
  cSpan <- 1/20 #Range of data to smooth for a point, higher = smoother
  cDegree <- 2 #1 linear fit, 2 for 2nd order polynomial
  cEvaluation <- length(reaCond) #Total number of point after loess smoothing
  cCondTH <- 0.05 #Threshold to count as ~0 for derivative
  near0Len <- 3
  near0Fac <- 10 #Relative change for start check and end check
  absChange <- 2 #Absolute change for start check and end check
  
  reaMeasCount <- seq(along = reaCond)
  #smooth the tracer data
  #Warnings are created sometimes when it runs against the edge of data, so supressing warnings from the smoothing
  suppressWarnings(cond.loess <- loess.smooth(reaMeasCount, reaCond, span = cSpan, degree = cDegree, evaluation = cEvaluation, family = "symmetric"))
  #plot(reaMeasCount, reaCond, xlab = "Measurement Number", ylab = "Specific Conductance")
  #lines(cond.loess$x,cond.loess$y, col = "green", lwd = 2)
  
  #do a quick and dirty derivative
  for(i in 1:length(cond.loess$x)-1){
    cond.loess$z[i] <- (cond.loess$y[i+1]-cond.loess$y[i])/(cond.loess$x[i+1]-cond.loess$x[i])
  }
  
  #Find a flat spot, peak, flat spot
  pkRange = 12 #Number of points to test for peak slope
  pkNum = 8 #Nummber of point that must be increasing/decreasing for peak slope
  fltRange = 60 #Number of points that must be flat
  fltFrac = 0.7 #Fraction of flat points that can be above cCondTH
  peakLoc <- min(which(cond.loess$z == max(cond.loess$z, na.rm = T)))
  testedPeaks <- peakLoc
  dataSeq <- seq(along = cond.loess$z)
  truePeak = F
  
  for(j in dataSeq){
    if(peakLoc <= pkRange | peakLoc <= pkNum){
      peakLoc <- min(which(cond.loess$z == max(cond.loess$z[!(dataSeq %in% testedPeaks)], na.rm = T)), na.rm = T)
      testedPeaks[length(testedPeaks) +1] <- peakLoc
      next
    }
    
    #check to see if the peak is too high (>2, would be more than 20 uS/sec change in conductance)
    if(cond.loess$z[peakLoc] > absChange){
      peakLoc <- min(which(cond.loess$z == max(cond.loess$z[!(dataSeq %in% testedPeaks)], na.rm = T)), na.rm = T)
      testedPeaks[length(testedPeaks) +1] <- peakLoc
      next
    }
    
    peakRiseA <- cond.loess$z[(peakLoc-pkRange):(peakLoc-1)]
    peakRiseB <- cond.loess$z[(peakLoc-pkRange+1):(peakLoc)]
    peakFallA <- cond.loess$z[(peakLoc+1):(peakLoc+pkRange)]
    peakFallB <- cond.loess$z[peakLoc:(peakLoc+pkRange-1)]
    
    riseTest <- mean(peakRiseB-peakRiseA > 0, na.rm = T) >= (pkNum/pkRange)
    fallTest <- mean(peakFallB-peakFallA > 0, na.rm = T) >= (pkNum/pkRange)
    posTest <- sum(peakRiseA, na.rm = T) > 0 & sum(peakRiseB, na.rm = T) > 0 & 
      sum(peakFallA, na.rm = T) > 0 & sum(peakFallB, na.rm = T) > 0
    stepTest <- mean(peakRiseA, na.rm = T) < 1 & mean(peakRiseB, na.rm = T) < 1 & 
      mean(peakFallA, na.rm = T) < 1 & mean(peakFallB, na.rm = T) < 1
    
    if(riseTest & fallTest & posTest & stepTest){
      #Ok, now look for flat spots for at least 10 minutes/60 measurements on either side
      fltStart <- max(which(abs(cond.loess$z[1:peakLoc]) < cCondTH))
      
      #If the peak left flat is close to the beginning of the data reduce fltRange
      if((fltStart - fltRange) <=0){
        fltRange <- fltStart
      }
      
      leftFlat <- sum(abs(cond.loess$z[(fltStart - fltRange):fltStart]) > cCondTH)/fltRange < fltFrac
      fltEnd <- min(which(abs(cond.loess$z[peakLoc:length(cond.loess$z)]) < cCondTH))
      
      if((fltEnd + fltRange + peakLoc) >= length(cond.loess$z)){
        peakLoc <- min(which(cond.loess$z == max(cond.loess$z[!(dataSeq %in% testedPeaks)], na.rm = T)), na.rm = T)
        testedPeaks[length(testedPeaks) +1] <- peakLoc
        next
      }
      
      rightFlat <- sum(abs(cond.loess$z[(fltEnd + peakLoc):(fltEnd + fltRange + peakLoc)]) > cCondTH)/fltRange < fltFrac
      
      #If it's flat break out, you found a real peak
      if((peakLoc - fltStart) <= fltRange & (fltEnd - peakLoc) <= fltRange & leftFlat & rightFlat){
        break
      }
    }
    #If we're at the end we didn't find a peak at all, warn
    if (length(testedPeaks) == length(dataSeq)){
      warning("Peak or rising limb of salt tracer was not detected")
      return(-9999)
    }
    #Move on to the next highest value that we haven't checked already
    peakLoc <- min(which(cond.loess$z == max(cond.loess$z[!(dataSeq %in% testedPeaks)], na.rm = T)), na.rm = T)
    testedPeaks[length(testedPeaks) +1] <- peakLoc
  }
  
  risingLimbIdx <- peakLoc

  #Break it into before and after the inflection point
  beforeInf <- cond.loess$z[1:risingLimbIdx]
  afterInf <- cond.loess$z[risingLimbIdx:length(cond.loess$z)]
  
  #This cuts off the data if there is a large step after the inflection point
  #The step size to care about is the inflection point of the derivative divided by the near0fac, default 10x
  endCheck <- -999
  for(i in (near0Len+1):length(afterInf)){
    if(i == length(afterInf)){
      endCheck <- i
    }else if(abs(mean(afterInf[(i-near0Len):i], na.rm = T)) > (cond.loess$z[risingLimbIdx]*near0Fac) |
             abs(mean(afterInf[(i-near0Len):i], na.rm = T)) > (cond.loess$z[risingLimbIdx] + absChange)){
      endCheck <- i-near0Len-2
      break
    }
  }
  
  #This cuts off the data if there is a large step before the inflection point
  #The step size to care about is the inflection point of thederivative divided by the near0fac, default 10x
  startCheck <- -999
  for(i in 1:(length(beforeInf)-near0Len)){
    if(i == (length(beforeInf)-near0Len)){
      startCheck <- 1
    }else if(abs(mean(beforeInf[(length(beforeInf)-i-near0Len):(length(beforeInf)-i)], na.em = T)) > (cond.loess$z[risingLimbIdx]*near0Fac) |
             abs(mean(beforeInf[(length(beforeInf)-i-near0Len):(length(beforeInf)-i)], na.em = T)) > (cond.loess$z[risingLimbIdx] + absChange)){
      startCheck <- length(beforeInf)-i+1
      break
    }
  }
  
  afterInf <- afterInf[1:endCheck]
  plateauIdx <- which(afterInf == min(afterInf[afterInf > 0 & !is.na(afterInf)], na.rm = T))
  plateauStart <- ifelse((plateauIdx - 5) < 1, (1 + risingLimbIdx), (plateauIdx - 5 + risingLimbIdx))
  plateauEnd <- ifelse((plateauIdx + 5) > length(cond.loess$y), length(cond.loess$y), (plateauIdx + 5 + risingLimbIdx))
  meanPlateauCond <- mean(cond.loess$y[plateauStart:plateauEnd], na.rm = T)
  
  #Don't consider the data before the start check
  beforeInf <- beforeInf[startCheck:length(beforeInf)]
  #Find the index of the baseline
  baselineIdx <- which(beforeInf == min(beforeInf[beforeInf > 0 & !is.na(beforeInf)], na.rm = T))
  baselineStart <- ifelse((baselineIdx - 5) < 1, (1 + startCheck), (baselineIdx - 5 + startCheck))
  baselineEnd <- ifelse((baselineIdx + 5) > length(cond.loess$y), length(cond.loess$y), (baselineIdx + 5 + startCheck))
  meanBaselineCond <- mean(cond.loess$y[baselineStart:baselineEnd], na.rm = T)
  
  condChange <- meanPlateauCond - meanBaselineCond
  
  halfCondDiff <- abs(cond.loess$y - (meanPlateauCond - 0.5 *condChange))
  minCount <- which(halfCondDiff == min(halfCondDiff[1:plateauStart], na.rm = T))
  
  #pkTime <- reaTime[which(reaMeasCount == ceiling(cond.loess$x[cond.loess$x == minCount]))]
  
  if(bPlot == TRUE){
    # #Plotting the whole data range for troubleshooting and debugging
    # invisible(dev.new(noRStudioGD = TRUE))
    # plot(reaMeasCount, reaCond, xlab = "Measurement Number", ylab = "Specific Conductance")
    # abline(h = meanBaselineCond, col = "red", lwd = 2)
    # abline(h = meanPlateauCond, col = "blue", lwd = 2)
    # title(main = plotName)
    # lines(cond.loess$x,cond.loess$y, col = "green", lwd = 2) #plot of loess smoothed line
    # par(new = TRUE) #For adding second axis
    # plot(cond.loess$x[1:length(cond.loess$z)],cond.loess$z, type = "o", xlab = '', ylab = '', col = "orange", yaxt = 'n', xaxt = 'n') #Derivative dots
    # axis(side = 4, at = pretty(range(cond.loess$z))) #Derivative axis on the right
    # abline(v = cond.loess$x[minCount], col = "blue", lwd = 2)
    # abline(v = cond.loess$x[startCheck], col = "purple", lwd = 2)
    # abline(v = cond.loess$x[endCheck], col = "red", lwd = 2)
    # abline(v = cond.loess$x[peakLoc], col = "green", lwd = 2)

    
    #Plotting just the range of data used for the travel time calculations
    invisible(dev.new(noRStudioGD = TRUE))
    plot(round(cond.loess$x[startCheck]):round(cond.loess$x[endCheck + risingLimbIdx]),
         reaCond[round(cond.loess$x[startCheck]):round(cond.loess$x[endCheck + risingLimbIdx])],
         xlab = "Measurement Number",
         ylab = "Specific Conductance", 
         ylim = c((meanBaselineCond-meanBaselineCond*0.1),(meanPlateauCond+meanPlateauCond*0.1)))
    abline(h = meanBaselineCond, col = "red", lwd = 2)
    abline(h = meanPlateauCond, col = "red", lwd = 2)
    title(main = plotName)
    mtext(paste0("blue line = half plateau point, red lines = background and plateau concentrations", "\n Click anywhere to close and continue"), cex = 0.8)
    #lines(cond.loess$x,cond.loess$y, col = "green", lwd = 2) #plot of loess smoothed line
    par(new = TRUE) #For adding second axis
    # plot(cond.loess$x[startCheck:(endCheck + risingLimbIdx)],
    #      cond.loess$z[startCheck:(endCheck + risingLimbIdx)],
    #      type = "o",
    #      xlab = '',
    #      ylab = '',
    #      col = "orange",
    #      yaxt = 'n',
    #      xaxt = 'n') #Derivative dots
    #      axis(side = 4, at = pretty(range(cond.loess$z))) #Derivative axis on the right
    abline(v = cond.loess$x[minCount], col = "blue", lwd = 2) #travel time line
    #print("Click anywhere on the plot to close and continue")
    ans <- identify(round(cond.loess$x[startCheck]):round(cond.loess$x[endCheck + risingLimbIdx]),
                    reaCond[round(cond.loess$x[startCheck]):round(cond.loess$x[endCheck + risingLimbIdx])], 
                    n = 1, 
                    tolerance = 1000, 
                    plot = F)
    invisible(dev.off())
    
  }
  pkIdx <- round(cond.loess$x[minCount])
  return(pkIdx)
}
