##############################################################################################
#' @title Constant rate or slug travel time between two stations

#' @author
#' Kaelin M. Cawley \email{kcawley@battelleecology.org} \cr

#' @description This function determines the timestamp of the peak tracer conductance for use
#' calculating travel time between an upstream and downstream sensor set.

#' @importFrom graphics points
#' @importFrom pracma trapz
#' @importFrom stats loess.smooth

#' @param loggerDataIn User input of the R data object holding the conductivity time series
#' for a site, date, and station [dataframe]
#' @param currEventID User input of the eventID of the tracer experiment [character]
#' @param injectionType User input of the injection type either "constant" or "slug" [character]
#' @param expStartTime User input of the experiment start time, in UTC timezone [posixct]
#' @param expBufferTime User input of the amount of time (in seconds) before the startTime
#' to show in plots, defaults to 600 seconds (10 minutes) [numeric]
#' @param backgroundCond User inpute of the background conductivity for making plots
#' look nicer [numeric]

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
  loggerDataIn,
  currEventID,
  injectionType,
  expStartTime,
  expBufferTime = 3900, #seconds
  backgroundCond = NA
){

  #Trim the data for only after the experiment started
  trimTime <- ifelse(min(loggerDataIn$dateTimeLogger) < (expStartTime - expBufferTime),
                                (expStartTime - expBufferTime),
                                min(loggerDataIn$dateTimeLogger))
  #trimTime <- min(loggerDataIn$dateTimeLogger)
  loggerDataTrimIn <- loggerDataIn[loggerDataIn$dateTimeLogger > trimTime,]
  #plot(loggerDataIn$dateTimeLogger, loggerDataIn$spCond)
  #lines(loggerDataTrimIn$dateTimeLogger, loggerDataTrimIn$spCond, col = "blue")

  # Need to pull in conductivity of the background water here, really
  #GUIL had an issue where the median reflected the concductivity in the air :(
  # #Create a plot where users select the range to pick the peak
  # medCond <- stats::median(loggerDataTrimIn$spCond, na.rm =TRUE)
  # stdCond <- stats::mad(loggerDataTrimIn$spCond, na.rm = TRUE)
  # lowPlot <- ifelse(medCond-30 < 0, 0, medCond-30)
  # highPlot <- ifelse(medCond+30 > max(loggerDataTrimIn$spCond, na.rm = TRUE),
  #                    max(loggerDataTrimIn$spCond, na.rm = TRUE),
  #                    medCond+30)

  if(nrow(loggerDataTrimIn) < 1){
    print(paste0("Trimmed logger data has no data for: ", currEventID))
    return(NULL)
  }
  
  if(is.na(backgroundCond) | is.nan(backgroundCond) | is.null(backgroundCond)){
    backgroundCond <- 0
  }
  
  
  lowPlot <- ifelse(min(loggerDataTrimIn$spCond, na.rm = TRUE) > backgroundCond - 50,
                    min(loggerDataTrimIn$spCond, na.rm = TRUE),
                    backgroundCond - 50)
  highPlot <- ifelse(max(loggerDataTrimIn$spCond, na.rm = TRUE) < backgroundCond + 50,
                     backgroundCond + 50,
                     max(loggerDataTrimIn$spCond, na.rm = TRUE))
  #Not trimming at all based on background conductivity
  # lowPlot <- min(loggerDataTrimIn$spCond, na.rm = TRUE)
  # highPlot <- max(loggerDataTrimIn$spCond, na.rm = TRUE)

  invisible(dev.new(noRStudioGD = TRUE))
  plot(loggerDataTrimIn$dateTimeLogger,
       loggerDataTrimIn$spCond,
       xlab = "Logger Date",
       ylab = "Specific Conductance",
       ylim = c(lowPlot, highPlot))

  #Have users choose if the plot has a defined peak
  points(x = c(min(loggerDataTrimIn$dateTimeLogger),max(loggerDataTrimIn$dateTimeLogger)),
         y = c(highPlot-(highPlot-lowPlot)*.2,highPlot-(highPlot-lowPlot)*.2),
         col = c("green", "red"),
         lwd = 2,
         pch = 19,
         cex = 2)
  
  if(injectionType == "NaBr"){
    titleText <- paste0("Click green dot (upper lefthand) if the peak is identifiable. \nClick red dot (upper righthand) if not identifiable.\n",
                        currEventID)
  }else{
    titleText <- paste0("Click green dot (upper lefthand) if the plateau rising limb is identifiable. \nClick red dot (upper righthand) if not identifiable.\n",
                        currEventID)
  }
  
  title(main = titleText)
  badPlotBox <- identify(x = c(min(loggerDataTrimIn$dateTimeLogger),max(loggerDataTrimIn$dateTimeLogger)),
                         y = c(highPlot-(highPlot-lowPlot)*.2,highPlot-(highPlot-lowPlot)*.2),
                         n = 1,
                         tolerance = 0.25,
                         labels = c("Good", "Bad"))
  Sys.sleep(1)
  invisible(dev.off())

  if(injectionType == "NaBr"){
    pickText <- paste0("Click left and right of peak. \n Keep at least the width of the '+' cursor on either side.\n", 
                       currEventID)
  }else{
    pickText <- paste0("Click left and right of plateau rising limb. \n Keep at least the width of the '+' cursor on either side.\n", 
                       currEventID)
  }
  
  if(length(badPlotBox) && badPlotBox==1){
    #If things look good, move on
    invisible(dev.new(noRStudioGD = TRUE))
    plot(loggerDataTrimIn$dateTimeLogger,
         loggerDataTrimIn$spCond,
         xlab = "Measurement Number",
         ylab = "Specific Conductance",
         ylim = c(lowPlot, highPlot))
    title(main = pickText)
    ans <- identify(x = loggerDataTrimIn$dateTimeLogger,
                    y = loggerDataTrimIn$spCond,
                    n = 2,
                    tolerance = 0.25)
    Sys.sleep(1)
    invisible(dev.off())
    beginHere <- min(ans)
    endHere <- max(ans)
  }else{
    return(NULL)
  }

  #Trim the loggerData to just the area specified
  loggerDataTrim <- loggerDataTrimIn[beginHere:endHere,]

  slugInjTypes <- c("NaBr","model","model - slug")
  criInjTypes <- c("NaCl","model - CRI")

  loggerDataTrim$originalSpCond <- loggerDataTrim$spCond
  # If it's a CRI convert to a peak by taking the derivative
  if(injectionType %in% criInjTypes){
    loggerDataTrim$derivative <- NA
    ##### Constants #####
    cSpan <- 1/10 #Range of data to smooth for a point, higher = smoother
    if(length(loggerDataTrim$spCond)<110){
      cSpan <- 0.17 # Need something in between for the shorter timeseries
    }
    if(length(loggerDataTrim$spCond)<60){
      cSpan <- 0.3 # Used to be 0.5 for things under 40
    }
    if(length(loggerDataTrim$spCond)<40){
      cSpan <- 0.5 # Used to be 0.5 for things under 40
    }
    cDegree <- 2 #1 linear fit, 2 for 2nd order polynomial
    cEvaluation <- length(loggerDataTrim$spCond) #Total number of point after loess smoothing
    #cCondTH <- 0.05 #Threshold to count as ~0 for derivative
    #near0Len <- 3
    #near0Fac <- 10 #Relative change for start check and end check
    #absChange <- 2 #Absolute change for start check and end check

    #Probably delete this reaMeasCount <- seq(along = loggerDataTrim$spCond)
    #smooth the tracer data
    #Warnings are created sometimes when it runs against the edge of data, so suppressing warnings from the smoothing
    suppressWarnings(cond.loess <- loess.smooth(loggerDataTrim$dateTimeLogger, loggerDataTrim$spCond, span = cSpan, degree = cDegree, evaluation = cEvaluation, family = "symmetric"))
    if(!exists("cond.loess")){
      print("Error finding peak time, smoothing could not be applied")
      return(NA)
    }
    plot(cond.loess$x, loggerDataTrim$originalSpCond, xlab = "Measurement Number", ylab = "Specific Conductance")
    lines(cond.loess$x, cond.loess$y, col = "green", lwd = 2)

    # #Manually smoothing the data
    # loggerDataTrim$test <- NA
    # for(i in 4:(length(loggerDataTrim$originalSpCond)-3)){
    #   loggerDataTrim$test[i] <- mean(loggerDataTrim$originalSpCond[(i-3):(i+3)])
    # }
    # loggerDataTrim$test[1:3] <- mean(loggerDataTrim$test[4:7], na.rm = TRUE)
    # loggerDataTrim$test[(length(loggerDataTrim$originalSpCond)-3):length(loggerDataTrim$originalSpCond)] <- mean(loggerDataTrim$test[(length(loggerDataTrim$originalSpCond)-7):(length(loggerDataTrim$originalSpCond)-3)], na.rm = TRUE)
    # plot(loggerDataTrim$dateTimeLogger, loggerDataTrim$originalSpCond, ylim = c(min(loggerDataTrim$originalSpCond, na.rm = TRUE), max(loggerDataTrim$originalSpCond, na.rm = TRUE)))
    # par(new = TRUE)
    # plot(loggerDataTrim$dateTimeLogger, loggerDataTrim$test, col = "blue", ylim = c(min(loggerDataTrim$originalSpCond, na.rm = TRUE), max(loggerDataTrim$originalSpCond, na.rm = TRUE)))

    #do a quick and dirty derivative
    for(i in 1:length(cond.loess$x)-1){
      cond.loess$z[i] <- (cond.loess$y[i+1]-cond.loess$y[i])/as.numeric(cond.loess$x[i+1]-cond.loess$x[i])
    }

    # peakLoc <- which(cond.loess$z == max(cond.loess$z,na.rm = T))
    # peakCondVal <- cond.loess$y[peakLoc]
    # peakLocOut <- which(abs(loggerDataTrim$spCond-peakCondVal) == min(abs(loggerDataTrim$spCond-peakCondVal))) + beginHere
    # #Handle when the peakTime is more than one value
    # if(length(peakLocOut)>1){
    #   peakLocOut <- round(mean(peakLocOut,na.rm = T), digits = 0)
    #   peakInfoOut <- list("peakLocOut"=peakLocOut,"peakStart"=beginHere,"peakEnd"=endHere)
    # }else{
    #   stop("Invalid injection type, stopping")
    # }
    #This is for the manual smoothing that isn't great
    # loggerDataTrim$derivative <- NA
    # for(i in 1:(length(loggerDataTrim$test)-1)){
    #   loggerDataTrim$derivative[i] <- (loggerDataTrim$test[i+1]-loggerDataTrim$test[i])/as.numeric(difftime(loggerDataTrim$dateTimeLogger[i+1],loggerDataTrim$dateTimeLogger[i], units = "min"))
    # }
#
#     plot(loggerDataTrim$dateTimeLogger, loggerDataTrim$derivative)
#     loggerDataTrim$spCond <- loggerDataTrim$derivative
#     loggerDataTrim$spCond[length(loggerDataTrim$spCond)] <- 0

    loggerDataTrim$derivative[1:length(cond.loess$z)] <- cond.loess$z
    loggerDataTrim$derivative[length(loggerDataTrim$derivative)] <- 0
    plot(loggerDataTrim$dateTimeLogger, loggerDataTrim$derivative)

    # loggerDataTrim$spCond[2:length(loggerDataTrim$spCond)] <- cond.loess$z
    # loggerDataTrim$spCond[1] <- mean(cond.loess$z[1:4])

  }else{
    loggerDataTrim$derivative <- loggerDataTrim$spCond
  }

  # Now that you have a peak, calculate the max location, temporal centroid, and harmonic mean
  peakTime <- loggerDataTrim$dateTimeLogger[which(loggerDataTrim$derivative == max(loggerDataTrim$derivative,na.rm = T))]
  #Handle when the peakTime is more than one value
  if(length(peakTime)>1){
    peakTime <- mean(peakTime,na.rm = TRUE)
  }

  #Now work on the ones the require integrating the area under the peak

  #Assume the first 5 points are all background concentrations
  backgroundConc <- mean(loggerDataTrim$derivative[1:5], na.rm = TRUE)

  #Background correct the logger data
  loggerDataTrim$corrSpCond <- loggerDataTrim$derivative - backgroundConc

  #Add the time difference column
  loggerDataTrim$timeDiff <- as.numeric(difftime(loggerDataTrim$dateTimeLogger, loggerDataTrim$dateTimeLogger[1], units = "sec"))

  #Now loop through and calculate the centriod time (tc) and harmonic mean time (thm)
  areaUnderCurve <- pracma::trapz(loggerDataTrim$timeDiff, loggerDataTrim$derivative)

  loggerDataTrim$tc <- NA
  loggerDataTrim$hm <- NA
  ti <- loggerDataTrim$timeDiff[1]
  tc <- 0
  hm <- 0
  for(i in 2:length(loggerDataTrim$derivative)){
    t <- loggerDataTrim$timeDiff[i]
    dt <- t - ti
    px <- (1/areaUnderCurve) * loggerDataTrim$derivative[i]

    tcToAdd <- t * px * dt
    loggerDataTrim$tc[i] <- tcToAdd
    tc <- tc + tcToAdd

    hmToAdd <- (1/t) * px * dt
    loggerDataTrim$hm[i] <- hmToAdd
    hm <- hm + hmToAdd

    ti <- t
  }
  thm <- 1/hm
  centroidTime <- loggerDataTrim$dateTimeLogger[1] + tc
  harmonicMeanTime <- loggerDataTrim$dateTimeLogger[1] + thm

  #Take a look at the outputs
  plot(loggerDataTrim$dateTimeLogger, loggerDataTrim$originalSpCond)
  lines(loggerDataTrim$dateTimeLogger, loggerDataTrim$test)
  par(new = TRUE)
  plot(loggerDataTrim$dateTimeLogger, loggerDataTrim$derivative, xlab = "", ylab = "", col = "blue", type = "l", axes = FALSE)
  abline(v = peakTime, col = "cyan")
  abline(v = tc+loggerDataTrim$dateTimeLogger[1], col = "green")
  abline(v = thm+loggerDataTrim$dateTimeLogger[1], col = "purple")
  graphics::legend(x = "right", legend = c("smoothed","derivative","peakTime","centroidTime","harmonicMeanTime"), lty = c(1,1,1,1), col = c("black","blue","cyan","green","purple"))

  peakInfoOut <- list("backgroundConc"=backgroundConc,
                      "peakArea"=areaUnderCurve,
                      "peakTime"=peakTime,
                      "centroidTime"=centroidTime,
                      "harmonicMeanTime"=harmonicMeanTime,
                      "startPlotTime"=loggerDataTrimIn$dateTimeLogger[beginHere],
                      "endPlotTime"=loggerDataTrimIn$dateTimeLogger[endHere])

  return(peakInfoOut)
}






