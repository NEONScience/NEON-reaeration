##############################################################################################
#' @title Constant rate or slug travel time between two stations

#' @author
#' Kaelin M. Cawley \email{kcawley@battelleecology.org} \cr

#' @description This function determines the timestamp of the peak tracer conductance for use
#' calculating travel time between an upstream and downstream sensor set.

#' @importFrom graphics points
#' @importFrom pracma trapz
#' @importFrom stats loess.smooth

#' @param loggerData User input of the R data object holding the conductivity time series
#' for a site, date, and station [dataframe]
#' @param currEventID User input of the eventID of the tracer experiment [string]
#' @param injectionType User input of the injection type either "constant" or "slug" [string]
#' @param expStartTime User input of the experiment start time, in UTC timezone [posixct]

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
#   Kaelin M. Cawley (2021-05-25)
#     updated to handle model injection types
##############################################################################################
def.calc.peakTime <- function(
  loggerData,
  currEventID,
  injectionType,
  expStartTime
){

  #Trim the data for only after the experiment started
  loggerData <- loggerData[loggerData$dateTimeLogger > expStartTime,]
  #plot(loggerData$dateTimeLogger, loggerData$spCond)
  #lines(trimmedLoggerData$dateTimeLogger, trimmedLoggerData$spCond, col = "blue")

  #Create a plot where users select the range to pick the peak
  medCond <- median(loggerData$spCond, na.rm =TRUE)
  stdCond <- mad(loggerData$spCond, na.rm = TRUE)
  lowPlot <- ifelse(medCond-30 < 0, 0, medCond-30)
  highPlot <- ifelse(medCond+30 > max(loggerData$spCond, na.rm = TRUE),
                     max(loggerData$spCond, na.rm = TRUE),
                     medCond+30)
  invisible(dev.new(noRStudioGD = TRUE))
  plot(loggerData$dateTimeLogger,
       loggerData$spCond,
       xlab = "Measurement Number",
       ylab = "Specific Conductance",
       ylim = c(lowPlot, highPlot))

  #Have users choose if the plot has a defined peak
  points(x = c(min(loggerData$dateTimeLogger),max(loggerData$dateTimeLogger)),
         y = c(highPlot*.8,highPlot*.8),
         col = c("green", "red"),
         lwd = 2,
         pch = 19,
         cex = 2)
  title(main = paste0("Click green dot (upper lefthand) if the peak/plateau is identifiable. \nClick red dot (upper righthand) if not identifiable.\n",currEventID))
  badPlotBox <- identify(x = c(min(loggerData$dateTimeLogger),max(loggerData$dateTimeLogger)),
                         y = c(highPlot*.8,highPlot*.8),
                         n = 1,
                         tolerance = 0.25,
                         labels = c("Good", "Bad"))
  Sys.sleep(1)
  invisible(dev.off())

  if(length(badPlotBox) && badPlotBox==1){
    #If things look good, move on
    invisible(dev.new(noRStudioGD = TRUE))
    plot(loggerData$dateTimeLogger,
         loggerData$spCond,
         xlab = "Measurement Number",
         ylab = "Specific Conductance",
         ylim = c(lowPlot, highPlot))
    title(main = paste0("The plot stars at the time of the exeriment.\nClick right of of peak/plateau to where the useful timeseries ends.\n"))
    ans <- identify(x = loggerData$dateTimeLogger,
                    y = loggerData$spCond,
                    n = 1,
                    tolerance = 0.25)
    Sys.sleep(1)
    invisible(dev.off())
    beginHere <- ans
  }else{
    return(NULL)
  }

  #Trim the loggerData to just the area specified
  loggerData <- loggerData[1:endHere,]

  slugInjTypes <- c("NaBr","model","model - slug")
  criInjTypes <- c("NaCl","model - CRI")

  loggerData$originalSpCond <- loggerData$spCond
  # If it's a CRI convert to a peak by taking the derivative
  if(injectionType %in% criInjTypes){
    ##### Constants #####
    cSpan <- 1/10 #Range of data to smooth for a point, higher = smoother
    if(length(loggerData$spCond)<60){
      cSpan <- 0.3 # Used to be 0.5 for things under 40
    }
    if(length(loggerData$spCond)<40){
      cSpan <- 0.5 # Used to be 0.5 for things under 40
    }
    cDegree <- 2 #1 linear fit, 2 for 2nd order polynomial
    cEvaluation <- length(loggerData$spCond) #Total number of point after loess smoothing
    #cCondTH <- 0.05 #Threshold to count as ~0 for derivative
    #near0Len <- 3
    #near0Fac <- 10 #Relative change for start check and end check
    #absChange <- 2 #Absolute change for start check and end check

    reaMeasCount <- seq(along = loggerData$spCond)
    #smooth the tracer data
    #Warnings are created sometimes when it runs against the edge of data, so supressing warnings from the smoothing
    suppressWarnings(cond.loess <- loess.smooth(reaMeasCount, loggerData$spCond, span = cSpan, degree = cDegree, evaluation = cEvaluation, family = "symmetric"))
    if(!exists("cond.loess")){
      print("Error finding peak time, smoothing could not be applied")
      return(NA)
    }
    #plot(reaMeasCount, loggerData$spCond, xlab = "Measurement Number", ylab = "Specific Conductance")
    #lines(cond.loess$x,cond.loess$y, col = "green", lwd = 2)

    #do a quick and dirty derivative
    for(i in 1:length(cond.loess$x)-1){
      cond.loess$z[i] <- (cond.loess$y[i+1]-cond.loess$y[i])/(cond.loess$x[i+1]-cond.loess$x[i])
    }
    # plot(cond.loess$x, cond.loess$y)
    # par(new = TRUE)
    # plot(cond.loess$x[2:length(cond.loess$x)],cond.loess$z, col = "green", ylim = c(-1,1))
    loggerData$spCond[2:length(loggerData$spCond)] <- cond.loess$z
    loggerData$spCond[1] <- mean(cond.loess$z[1:4])

  }

  # Now that you have a peak, calculate the max location, temporal centroid, and harmonic mean
  peakTime <- loggerData$dateTimeLogger[which(loggerData$spCond == max(loggerData$spCond,na.rm = T))]
  #Handle when the peakTime is more than one value
  if(length(peakTime)>1){
    peakTime <- mean(peakTime,na.rm = TRUE)
  }

  #Now work on the ones the require integrating the area under the peak

  #Assume the first 5 points are all background concentrations
  backgroundConc <- mean(loggerData$spCond[1:5], na.rm = TRUE)

  #Background correct the logger data
  loggerData$corrSpCond <- loggerData$spCond - backgroundConc

  #Add the time difference column
  loggerData$timeDiff <- as.numeric(difftime(loggerData$dateTimeLogger, loggerData$dateTimeLogger[1], units = "sec"))

  #Now loop through and calculate the centriod time (tc) and harmonic mean time (thm)
  areaUnderCurve <- pracma::trapz(loggerData$timeDiff, loggerData$spCond)

  loggerData$tc <- NA
  loggerData$hm <- NA
  ti <- loggerData$timeDiff[1]
  tc <- 0
  hm <- 0
  for(i in 2:length(loggerData$spCond)){
    t <- loggerData$timeDiff[i]
    dt <- t - ti
    px <- (1/areaUnderCurve) * loggerData$spCond[i]

    tcToAdd <- t * px * dt
    loggerData$tc[i] <- tcToAdd
    tc <- tc + tcToAdd

    hmToAdd <- (1/t) * px * dt
    loggerData$hm[i] <- hmToAdd
    hm <- hm + hmToAdd

    ti <- t
  }
  thm <- 1/hm
  centroidTime <- loggerData$dateTimeLogger[1] + tc
  harmonicMeanTime <- loggerData$dateTimeLogger[1] + thm

  # #Take a look at the outputs
  # plot(loggerData$dateTimeLogger, loggerData$originalSpCond)
  # par(new = TRUE)
  # plot(loggerData$dateTimeLogger, loggerData$spCond, xlab = "", ylab = "", col = "blue", type = "l", axes = FALSE)
  # abline(v = peakTime, col = "cyan")
  # abline(v = tc+loggerData$dateTimeLogger[1], col = "green")
  # abline(v = thm+loggerData$dateTimeLogger[1], col = "purple")

  peakInfoOut <- list("peakTime"=peakTime,"peakArea"=areaUnderCurve,"centroidTime"=endHere, "harmonicMeanTime"=areaResponse)

  return(peakInfoOut)
}






