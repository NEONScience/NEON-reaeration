##############################################################################################
#' @title Constant rate or slug travel time between two stations

#' @author 
#' Kaelin M. Cawley \email{kcawley@battelleecology.org} \cr

#' @description This function determines the timestamp of the peak tracer conductance for use 
#' calculating travel time between an upstream and downstream sensor set.

#' @importFrom neonUtilities stackByTable
#' @importFrom graphics points
#' @importFrom pracma trapz
#' @importFrom stats loess.smooth

#' @param loggerData User input of the R data object holding the conductivity time series 
#' for a site, date, and station [string]
#' @param currEventID User input of the eventID of the tracer experiment [string]
#' @param injectionType User input of the injection type either "constant" or "slug" [string]

#' @return This function returns the peak tracer timestamp [dateTime]

#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

#' @keywords surface water, streams, rivers, velocity, travel time, reaeration, metabolism

#' @examples
#' #Using an example file
#' #travelTime <- def.calc.travelTime(
#' #dataDir = paste(path.package("reaRate"),"inst\\extdata", sep = "\\"), 
#' #currEventID = "GUIL.20150129", injectionType = "constant", bPlot = T)

#' @seealso def.calc.reaeration.R for calculating reaeration rates

#' @export

# changelog and author contributions / copyrights
#   Kaelin M. Cawley (2017-08-03)
#     original creation
##############################################################################################
def.calc.peakTime <- function(
  loggerData,
  currEventID,
  injectionType
  ){
  
  #Create a plot where users select the range to pick the peak
  invisible(dev.new(noRStudioGD = TRUE))
  plot(seq(along=loggerData), loggerData, xlab = "Measurement Number", ylab = "Specific Conductance")
  
  #Have users choose if the plot has a defined peak
  points(x = c(length(loggerData)*.1,length(loggerData)*.9),
         y = c(max(loggerData,na.rm = T)*.8,max(loggerData,na.rm = T)*.8), 
         col = c("green","red"), 
         lwd = 2, 
         pch = 19,
         cex = 2)
  title(main = paste0("Click green dot (upper lefthand) if the peak/plateau is identifiable. \nClick red dot (upper righthand) if not identifiable.\n",currEventID))
  badPlotBox <- identify(x = c(length(loggerData)*.1,length(loggerData)*.9),
                         y = c(max(loggerData,na.rm = T)*.8,max(loggerData,na.rm = T)*.8), 
                         n = 1, 
                         tolerance = 0.25, 
                         labels = c("Good", "Bad"))
  Sys.sleep(1)
  invisible(dev.off())
  
  if(badPlotBox==1){
    #If things look good, move on
    invisible(dev.new(noRStudioGD = TRUE))
    plot(seq(along=loggerData), loggerData, xlab = "Measurement Number", ylab = "Specific Conductance")
    title(main = paste0("Click left and right of of peak/plateau. \n Keep at least the width of the '+' cursor on either side.\n"))
    ans <- identify(seq(along=loggerData), n = 2, loggerData, tolerance = 0.25)
    Sys.sleep(1)
    invisible(dev.off())
    beginHere <- min(ans)
    endHere <- max(ans)
  }else{
    return(NULL)
  }
  
  #Trim the loggerData to just the area specified
  loggerData <- loggerData[beginHere:endHere]
  
  #Slug injections, find the index/measurement number of the tracer peak
  if(injectionType == "NaBr"){
    peakLoc <- which(loggerData == max(loggerData,na.rm = T))
    #Handle when the peakTime is more than one value
    if(length(peakLoc)>1){
      peakLoc <- round(mean(peakLoc,na.rm = T), digits = 0)
    }
    peakLocOut <- peakLoc + beginHere
    
    #Background correct the logger data
    loggerData <- loggerData - mean(loggerData[1:5])
    #Calculate the area under the conductivity time series
    areaResponse <- trapz(1:length(loggerData), loggerData)
    #Convert from integration over measurement number to time
    areaResponse <- areaResponse*10 #Each measurement is 10 seconds apart
    
    peakInfoOut <- list("peakLocOut"=peakLocOut,"peakStart"=beginHere,"peakEnd"=endHere, "peakArea"=areaResponse)
  }else if(injectionType == "NaCl"){
    ##### Constants #####
    cSpan <- 1/10 #Range of data to smooth for a point, higher = smoother
    if(length(loggerData)<40){
      cSpan <- 0.5
    }
    cDegree <- 2 #1 linear fit, 2 for 2nd order polynomial
    cEvaluation <- length(loggerData) #Total number of point after loess smoothing
    #cCondTH <- 0.05 #Threshold to count as ~0 for derivative
    #near0Len <- 3
    #near0Fac <- 10 #Relative change for start check and end check
    #absChange <- 2 #Absolute change for start check and end check
    
    reaMeasCount <- seq(along = loggerData)
    #smooth the tracer data
    #Warnings are created sometimes when it runs against the edge of data, so supressing warnings from the smoothing
    suppressWarnings(cond.loess <- loess.smooth(reaMeasCount, loggerData, span = cSpan, degree = cDegree, evaluation = cEvaluation, family = "symmetric"))
    if(!exists("cond.loess")){
      print("Error finding peak time, smoothing could not be applied")
      return(NA)
    }
    #plot(reaMeasCount, loggerData, xlab = "Measurement Number", ylab = "Specific Conductance")
    #lines(cond.loess$x,cond.loess$y, col = "green", lwd = 2)
    
    #do a quick and dirty derivative
    for(i in 1:length(cond.loess$x)-1){
      cond.loess$z[i] <- (cond.loess$y[i+1]-cond.loess$y[i])/(cond.loess$x[i+1]-cond.loess$x[i])
    }
    
    peakLoc <- which(cond.loess$z == max(cond.loess$z,na.rm = T))
    
    peakCondVal <- cond.loess$y[peakLoc]
    peakLocOut <- which(abs(loggerData-peakCondVal) == min(abs(loggerData-peakCondVal))) + beginHere
    #Handle when the peakTime is more than one value
    if(length(peakLocOut)>1){
      peakLocOut <- round(mean(peakLocOut,na.rm = T), digits = 0)
    }
    
    peakInfoOut <- list("peakLocOut"=peakLocOut,"peakStart"=beginHere,"peakEnd"=endHere)
  }
  
  return(peakInfoOut)
}

  



                                           