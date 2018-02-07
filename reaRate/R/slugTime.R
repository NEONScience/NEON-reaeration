##############################################################################################
#' @title Slug peak concentration location

#' @author 
#' Kaelin M. Cawley \email{kcawley@battelleecology.org} \cr

#' @description This function determines the location of the peak specific conductance for use 
#' calculating travel time between an upstream and downstream sensor set.

#' @param condData data to be fit for peak
#' @param pkRange Number of points on either side of peak to test for slope  [integer]
#' @param pkNum Number of points that must be increasing/decreasing for peak slope [integer]

#' @return This function returns the peak tracer index [integer]

#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

#' @keywords surface water, streams, rivers, velocity, travel time, reaeration, metabolism

#' @examples
#' #Using an example file
#' s4peakLoc <- slugTime(condDataS4)

#' @seealso def.calc.reaeration.R for calculating reaeration rates and 
#' def.calc.travelTime.R for determining travel times

#' @export

# changelog and author contributions / copyrights
#   Kaelin M. Cawley (2017-11-01)
#     original creation
##############################################################################################
slugTime <- function(
  condData,
  pkRange = 5, #Number of points to test for peak slope
  pkNum = 4 #Nummber of point that must be increasing/decreasing for peak slope
){
  
  #Find slug peak, start with the max and look for 5 on either side that decrease
  peakLoc <- min(which(condData == max(condData, na.rm = T)))
  testedPeaks <- peakLoc
  dataSeq <- seq(along = condData)
  
  for(j in dataSeq){
    if(peakLoc <= pkRange | peakLoc <= pkNum){
      peakLoc <- min(which(condData == max(condData[!(dataSeq %in% testedPeaks)], na.rm = T)), na.rm = T)
      testedPeaks[length(testedPeaks) +1] <- peakLoc
      next
    }
    peakRiseA <- condData[(peakLoc-pkRange):(peakLoc-1)]
    peakRiseB <- condData[(peakLoc-pkRange+1):(peakLoc)]
    peakFallA <- condData[(peakLoc+1):(peakLoc+pkRange)]
    peakFallB <- condData[peakLoc:(peakLoc+pkRange-1)]
    
    riseTest <- mean(peakRiseB-peakRiseA > 0, na.rm = T) >= (pkNum/pkRange)
    fallTest <- mean(peakFallB-peakFallA > 0, na.rm = T) >= (pkNum/pkRange)
    posTest <- sum(peakRiseA, na.rm = T) > 0 & sum(peakRiseB, na.rm = T) > 0 & 
      sum(peakFallA, na.rm = T) > 0 & sum(peakFallB, na.rm = T) > 0
    stepTest <- mean(peakRiseA, na.rm = T) < 1 & mean(peakRiseB, na.rm = T) < 1 & 
      mean(peakFallA, na.rm = T) < 1 & mean(peakFallB, na.rm = T) < 1
    
    if(riseTest & fallTest & posTest & stepTest){
      break
    }else{
      if (length(testedPeaks) == length(dataSeq)){
        warning("Peak or rising limb of dalt tracer was not detected")
        return(-9999)
      }
      peakLoc <- min(which(condData == max(condData[!(dataSeq %in% testedPeaks)], na.rm = T)), na.rm = T)
      testedPeaks[length(testedPeaks) +1] <- peakLoc
    }
    
  }
  return(peakLoc)
}

