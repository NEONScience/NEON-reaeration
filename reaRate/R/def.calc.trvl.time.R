##############################################################################################
#' @title Travel time and mean depth calculations

#' @author
#' Kaelin M. Cawley \email{kcawley@battelleecology.org} \cr

#' @description This function calculates travel time and mean depth from reaeration experiments.

#' @importFrom grDevices dev.new
#' @importFrom grDevices dev.off
#' @importFrom grDevices dev.copy
#' @importFrom grDevices png
#' @importFrom grDevices jpeg
#' @importFrom graphics identify
#' @importFrom graphics abline
#' @importFrom graphics axis
#' @importFrom graphics lines
#' @importFrom graphics mtext
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics title
#' @importFrom stats lm
#' @importFrom stats lsfit
#' @importFrom methods is

#' @param inputFile Name of the data frame containing the information needed to calculate the
#' reaeration parameters. If the headers are named: "injectionType", "eventID",
#' "stationToInjectionDistance", "plateauGasConc", "corrPlatSaltConc", "hoboSampleID",
#' "wettedWidth", respectively, no other inputs are required. Otherwise, the names of the
#' columns need to be input for the function to work. [string]
#' @param loggerData User identified filename of logger data [string]
#' @param namedLocation A string identifier for the station where data was collected [string]
#' @param injectionTypeName Either constant rate or slug [string]
#' @param eventID A string identifier to link records collected as part of the same experiment,
#' SITE.YYYYMMDD for NEON [string]
#' @param slopeRaw Dataframe column name for SF6 slope from raw data [string]
#' @param slopeClean Dataframe column name for SF6 slope from data with outliers removed [string]
#' @param slopeSaltCorr Dataframe column name for SF6 slope from plateau salt corrected data 
#' with outliers removed [string]
#' @param slopeBackCorr Dataframe column name for SF6 slope from background corrected plateau 
#' salt data with outliers removed [string]
#' @param stationToInjectionDistance Dataframe column name for distance from station to
#' injection [string]
#' @param hoboSampleID Dataframe column name for ID to link to conductivity timeseries data [string]
#' @param slugPourTime Dataframe column name for dateTime when slug was poured [string]
#' @param dripStartTime Dataframe column name for dateTiem when drip was started [string]
#' @param meanBackgroundCond Dataframe column name for mean background specific conductance [string]
#' @param discharge Dataframe column name for stream discharge in literPerSecond [string]
#' @param waterTemp Dataframe column name for mean water temperature data [string]
#' @param wettedWidth Dataframe column name for mean wetted width for the stream reach [string]
#' @param plot User input to plot the SF6/corrected salt concentration versus distance downstream,
#' defaults to TRUE [boolean]
#' @param savePlotPath If a user specifies a path the plots will be saved to this location [string]

#' @return This function returns a list of two dataframes, the input dataframe of data for up to
#' 4 stations per site per date and an output dataframe appended with travel time and mean depth
#' for a given site and date

#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

#' @keywords surface water, streams, rivers, reaeration, deaeration, SF6, metabolism, tracer

#' @examples
#' #TBD

#' @seealso def.calc.peakTime for calculating travel times and def.format.reaeration for
#' formatting reaeration data

#' @export

# changelog and author contributions / copyrights
#   Kaelin M. Cawley (2022-03-02)
#     original creation
##############################################################################################
#This code is for calculating reaeration rates and Schmidt numbers
def.calc.trvl.time <- function(
  inputFile = NULL,
  loggerData = NULL,
  namedLocation = "namedLocation",
  injectionTypeName = "injectionType",
  eventID = "eventID",
  slopeRaw = 'slopeRaw',
  slopeClean = 'slopeClean',
  slopeSaltCorr = 'slopeSaltCorr',
  slopeBackCorr = 'slopeBackCorr',
  stationToInjectionDistance = "stationToInjectionDistance",
  hoboSampleID = "hoboSampleID",
  slugPourTime = "slugPourTime",
  dripStartTime = "dripStartTime",
  meanBackgroundCond = "meanBackgroundCond",
  discharge = "fieldDischarge_lps",
  waterTemp = "waterTemp",
  wettedWidth = "wettedWidth",
  plot = TRUE,
  savePlotPath = NULL
){
  
  if(!plot && !is.null(savePlotPath)){
    stop("Please turn plotting on (plot = T) in order to save plots.")
  }
  
  namLocIdx <- which(names(inputFile) == namedLocation)
  injTypeIdx <- which(names(inputFile) == injectionTypeName)
  eventIDIdx <- which(names(inputFile) == eventID)
  slopeRawIdx <- which(names(inputFile) == slopeRaw)
  slopeCleanIdx <- which(names(inputFile) == slopeClean)
  slopeSaltCorrIdx <- which(names(inputFile) == slopeSaltCorr)
  slopeBackCorrIdx <- which(names(inputFile) == slopeBackCorr)
  staDistIdx <- which(names(inputFile) == stationToInjectionDistance)
  loggerIdx <- which(names(inputFile) == hoboSampleID)
  slugTimeIdx <- which(names(inputFile) == slugPourTime)
  dripStartIdx <- which(names(inputFile) == dripStartTime)
  backCondIdx <- which(names(inputFile) == meanBackgroundCond)
  QIdx <- which(names(inputFile) == discharge)
  watTempIdx <- which(names(inputFile) == waterTemp)
  wwIdx <- which(names(inputFile) == wettedWidth)
  
  ##### Constants #####
  convLpsCms = 1/1000 #Conversion from litersPerSecond to cubicMetersPerSecond
  
  #Create output file
  outputDFNames <- c(
    'siteID',
    'eventID',
    'slopeRaw',
    'slopeClean',
    'slopeSaltCorr',
    'slopeBackCorr',
    # 'noSaltCorrLossRateSF6',
    # 'platSaltCorrLossRateSF6',
    # 'backSaltCorrLossRateSF6',
    'S1PeakTime',
    'S4PeakTime',
    'peakMaxTravelTime',
    # 'S1CentroidTime', # These aren't really worth doing with the peak detection algorithm we have
    # 'S4CentroidTime',
    # 'centroidTravelTime',
    # 'S1HarmonicMeanTime',
    # 'S4HarmonicMeanTime',
    # 'harmonicMeanTravelTime',
    'btwStaDist',
    'peakMaxVelocity',
    # 'centroidVelocity',
    # 'harmonicMeanVelocity',
    'meanDepth',
    'meanQ_lps',
    'meanQ_cms',
    'meanTemp'
  )
  
  #Only use the unique eventIDs
  allEventID <- unique(inputFile[[eventIDIdx]])
  outputDF <- data.frame(matrix(data=NA, ncol=length(outputDFNames), nrow=length(allEventID)))
  names(outputDF) <- outputDFNames
  
  outputDF$eventID <- unique(inputFile[[eventIDIdx]])
  
  #Check for correct date format
  if(all(grepl("20[0-9]{2}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z",loggerData$dateTimeLogger))){
    dateFormat <- "%Y-%m-%dT%H:%M:%SZ"
    loggerData$dateTimeLogger <- as.POSIXct(loggerData$dateTimeLogger, format = dateFormat, tz = "UTC")
  }else if(all(grepl("20[0-9]-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}\\.000\\+0000",loggerData$dateTimeLogger))){
    dateFormat <- "%Y-%m-%dT%H:%M:%S.000+0000"
    loggerData$dateTimeLogger <- as.POSIXct(loggerData$dateTimeLogger, format = dateFormat, tz = "UTC")
  }else if(!"POSIXct" %in% is(loggerData$dateTimeLogger)){
    stop("Inconsistent or unidentified date formats in conductivity logger data.")
  }
  
  for(i in seq(along = outputDF$eventID)){
    #for(i in 21:22){
    
    modelInjType <- FALSE
    currEventID <- outputDF$eventID[i]
    #Uncomment this if you'd like to see a list of all the eventIDs for troubleshooting or debugging
    #print(paste0(i, " - ", currEventID))
    injectionType <- unique(inputFile[inputFile[[eventIDIdx]] == currEventID & !is.na(inputFile[[injTypeIdx]]), injTypeIdx])
    if(length(injectionType)<1){
      cat("Warning - Injection type unknown for",currEventID,"\n")
      next
    }
    #Calculations for the "model" slug injections only TBD
    #For the moment just skip those
    if(injectionType %in% c("model","model - slug","model - CRI")){
      print(paste0("Model injection type, cannot calculate loss rate for ", currEventID))
      modelInjType <- TRUE
    }
    
    #Use drip of slug time for the experiment start time
    slugTime <- unique(inputFile[inputFile[[eventIDIdx]] == currEventID, slugTimeIdx])
    injTime <- unique(inputFile[inputFile[[eventIDIdx]] == currEventID, dripStartIdx])
    if(is.na(slugTime)){
      currExpStartTime <- unique(inputFile[inputFile[[eventIDIdx]] == currEventID, dripStartIdx])
    }else{
      currExpStartTime <- unique(inputFile[inputFile[[eventIDIdx]] == currEventID, slugTimeIdx])
    }
    
    if(is.na(currExpStartTime)){
      cat('Warning, experiment startTime could not be determined for',currEventID)
    }
    
    outputDF$siteID[i] <- unique(substr(inputFile[[namLocIdx]][inputFile[[eventIDIdx]] == currEventID], 1, 4))
    S1 <- paste(outputDF$siteID[i], "AOS.reaeration.station.01", sep = ".")
    S2 <- paste(outputDF$siteID[i], "AOS.reaeration.station.02", sep = ".")
    S3 <- paste(outputDF$siteID[i], "AOS.reaeration.station.03", sep = ".")
    S4 <- paste(outputDF$siteID[i], "AOS.reaeration.station.04", sep = ".")
    
    backCond <- mean(inputFile[inputFile$eventID == currEventID, backCondIdx], na.rm = TRUE)

    #New section that requires the user to pick the range of data for the peak or plateau rising limb
    #currEventID <- currEventID
    s1LoggerData <- loggerData[loggerData$hoboSampleID == paste0(substr(currEventID, 1, 4), "_S1_", substr(currEventID, 6, 13)),]
    s1LoggerData <- s1LoggerData[order(s1LoggerData$measurementNumber),]
    s2LoggerDeployed <- FALSE
    
    s4LoggerData <- loggerData[loggerData$hoboSampleID == paste0(substr(currEventID, 1, 4), "_S4_", substr(currEventID, 6, 13)),]
    s4LoggerData <- s4LoggerData[order(s4LoggerData$measurementNumber),]
    
    if(length(s1LoggerData[[1]]) <= 0 & length(s4LoggerData[[1]]) <= 0){
      print(paste0("Conductivity logger data not available for ", currEventID, ", stations S1 & S4"))
      next
    }else if(length(s1LoggerData[[1]]) <= 0){
      #This is added in for times when the first logger is at station 2 instead of station 1
      print(paste0("Conductivity logger data not available for ", currEventID, ", station S1, looking for S2 logger data."))
      
      s2LoggerData <- loggerData[loggerData$hoboSampleID == paste0(substr(currEventID, 1, 4), "_S2_", substr(currEventID, 6, 13)),]
      s2LoggerData <- s2LoggerData[order(s2LoggerData$measurementNumber),]
      
      if(length(s2LoggerData[[1]]) <= 0|length(s2LoggerData[[1]]) < 10){
        print(paste0("Conductivity logger data not available/sufficient for ", currEventID, ", station S2"))
        next
      }else{
        print(paste0("Conductivity logger data found available for ", currEventID, ", station S2"))
        s1LoggerData <- s2LoggerData
        s2LoggerDeployed <- TRUE
      }
      
    }else if(length(s4LoggerData[[1]]) <= 0){
      print(paste0("Conductivity logger data not available for ", currEventID, ", station S4"))
      next
    }
    
    if(length(s1LoggerData[[1]]) < 10){
      print(paste0("Conductivity logger data has less than ten points for ", currEventID, ", station S1"))
      next
    }else if(length(s4LoggerData[[1]]) < 10){
      print(paste0("Conductivity logger data has less than ten points for ", currEventID, ", station S4"))
      next
    }
    
    #If low range isn't collected use the full range
    if(!all(is.na(s1LoggerData$lowRangeSpCondNonlinear))){
      condDataS1 <- s1LoggerData[,c("dateTimeLogger","lowRangeSpCondNonlinear")]
      s1RangeFull <- FALSE
    }else if(!all(is.na(s1LoggerData$fullRangeSpCondNonlinear))){
      condDataS1 <- s1LoggerData[,c("dateTimeLogger","fullRangeSpCondNonlinear")]
      s1RangeFull <- TRUE
    }else{
      print(paste0("Conductivity logger data not available for ", currEventID, ", station S1"))
      next
    }
    
    if(!all(is.na(s4LoggerData$lowRangeSpCondNonlinear))){
      condDataS4 <- s4LoggerData[,c("dateTimeLogger","lowRangeSpCondNonlinear")]
      s4RangeFull <- FALSE
    }else if(!all(is.na(s4LoggerData$fullRangeSpCondNonlinear))){
      condDataS4 <- s4LoggerData[,c("dateTimeLogger","fullRangeSpCondNonlinear")]
      s4RangeFull <- TRUE
    }else{
      print(paste0("Conductivity logger data not available for ", currEventID, ", station S4"))
      next
    }
    names(condDataS1) <- c("dateTimeLogger","spCond")
    names(condDataS4) <- c("dateTimeLogger","spCond")
    
    if(is.na(currExpStartTime)){
      currExpStartTime <- max(condDataS1$dateTimeLogger[1], condDataS4$dateTimeLogger[1])
    }
    
    #Find the peak locations
    s1peakLoc <- reaRate::def.calc.peakTime(loggerDataIn = condDataS1,
                                            currEventID = currEventID,
                                            injectionType = injectionType,
                                            expStartTime = currExpStartTime,
                                            backgroundCond = backCond) # index to get date and time of peak/plateau half max
    s4peakLoc <- reaRate::def.calc.peakTime(loggerDataIn = condDataS4,
                                            currEventID = currEventID,
                                            injectionType = injectionType,
                                            expStartTime = currExpStartTime,
                                            backgroundCond = backCond) # index to get date and time of peak/plateau half max
    
    #If either of the peakTimes are NULL move on to the next eventID
    if(is.null(s1peakLoc)){
      print(paste0("Conductivity logger data peak/plateau cannot be identified for ", currEventID, ", station S1"))
      next
    }
    if(is.null(s4peakLoc)){
      print(paste0("Conductivity logger data peak/plateau cannot be identified for ", currEventID, ", station S4"))
      next
    }
    
    #Get the dates from the indices and subtract to get travel time
    outputDF$S1PeakTime[i] <- s1peakLoc$peakTime
    outputDF$S4PeakTime[i] <- s4peakLoc$peakTime
    outputDF$peakMaxTravelTime[i] <- difftime(s4peakLoc$peakTime,
                                              s1peakLoc$peakTime,
                                              units = "secs")
    
    # #Get the dates from the indices and subtract to get travel time
    # outputDF$S1CentroidTime[i] <- s1peakLoc$centroidTime
    # outputDF$S4CentroidTime[i] <- s4peakLoc$centroidTime
    # outputDF$centroidTravelTime[i] <- difftime(s4peakLoc$centroidTime,
    #                                            s1peakLoc$centroidTime,
    #                                            units = "secs")
    # 
    # #Get the dates from the indices and subtract to get travel time
    # outputDF$S1HarmonicMeanTime[i] <- s1peakLoc$harmonicMeanTime
    # outputDF$S4HarmonicMeanTime[i] <- s4peakLoc$harmonicMeanTime
    # outputDF$harmonicMeanTravelTime[i] <- difftime(s4peakLoc$harmonicMeanTime,
    #                                                s1peakLoc$harmonicMeanTime,
    #                                                units = "secs")
    
    #Plot the travel times to check
    if(plot==T){
      # if(s1RangeFull){
      #   s1YData <- s1LoggerData$fullRangeSpCondNonlinear[s1peakLoc$peakStart:s1peakLoc$peakEnd]
      # }else{
      #   s1YData <- s1LoggerData$lowRangeHobo[s1peakLoc$peakStart:s1peakLoc$peakEnd]
      # }
      #
      # if(s4RangeFull){
      #   s4YData <- s4LoggerData$fullRangeSpCondNonlinear[s4peakLoc$peakStart:s4peakLoc$peakEnd]
      # }else{
      #   s4YData <- s4LoggerData$lowRangeHobo[s4peakLoc$peakStart:s4peakLoc$peakEnd]
      # }
      
      s1YData <- condDataS1$spCond[condDataS1$dateTimeLogger > s1peakLoc$startPlotTime & condDataS1$dateTimeLogger < s1peakLoc$endPlotTime]
      s4YData <- condDataS4$spCond[condDataS4$dateTimeLogger > s4peakLoc$startPlotTime & condDataS4$dateTimeLogger < s4peakLoc$endPlotTime]
      
      invisible(dev.new(noRStudioGD = TRUE))
      x <- condDataS1$dateTimeLogger[condDataS1$dateTimeLogger > s1peakLoc$startPlotTime & condDataS1$dateTimeLogger < s1peakLoc$endPlotTime]
      #y <- s1LoggerData$fullRangeSpCondNonlinear[s1peakLoc$peakStart:s1peakLoc$peakEnd]
      minTime <- min(s1peakLoc$startPlotTime,s4peakLoc$startPlotTime)
      maxTime <- max(s1peakLoc$endPlotTime,s4peakLoc$endPlotTime)
      minY <- min(s1YData,s4YData,na.rm = TRUE)
      maxY <- max(s1YData,s4YData,na.rm = TRUE)
      
      #Save out plot of loss rate to specified directory
      # if(!is.null(savePlotPath)){
      #   png(paste0(savePlotPath,"/travelTime_",currEventID,".png"))
      #   plot(x,
      #        s1YData,
      #        xlim = c(minTime,maxTime),
      #        ylim = c(minY,maxY),
      #        ylab = "Conductivity, uS",
      #        xlab = "Time (UTC)")
      #   mtext(paste0("Travel Time = ",outputDF$peakMaxTravelTime[i]," seconds, (",round(as.numeric(outputDF$peakMaxTravelTime[i])/60,digits=1) ," min)\n Click anywhere to close and continue"), cex = 1.2)
      #   points(condDataS4$dateTimeLogger[condDataS4$dateTimeLogger > s4peakLoc$startPlotTime & condDataS4$dateTimeLogger < s4peakLoc$endPlotTime],
      #          s4YData,
      #          col = "blue")
      #   abline(v = s1peakLoc$peakTime)
      #   abline(v = s4peakLoc$peakTime, col = "blue")
      #   graphics::legend(x = "bottomright", legend = c("upstream","downstream"), lty = c(1,1), col = c("black","blue"))
      #   dev.off()
      # }
      plot(x,
           s1YData,
           xlim = c(minTime,maxTime),
           ylim = c(minY,maxY),
           ylab = "Conductivity, uS",
           xlab = "Time (UTC)")
      mtext(paste0("Travel Time = ",outputDF$peakMaxTravelTime[i]," seconds, (",round(as.numeric(outputDF$peakMaxTravelTime[i])/60,digits=1) ," min)\n Click anywhere to close and continue\n", currEventID), cex = 1.2)
      points(condDataS4$dateTimeLogger[condDataS4$dateTimeLogger > s4peakLoc$startPlotTime & condDataS4$dateTimeLogger < s4peakLoc$endPlotTime],
             s4YData,
             col = "blue")
      abline(v = s1peakLoc$peakTime)
      abline(v = s4peakLoc$peakTime, col = "blue")
      graphics::legend(x = "bottomright", legend = c("upstream","downstream"), lty = c(1,1), col = c("black","blue"))
      ans <- identify(x, s1YData, n = 1, tolerance = 100, plot = F)
      if(!is.null(savePlotPath)){
        dev.copy(jpeg,paste0(savePlotPath,"/travelTime_",currEventID,".jpg"))
      }
      invisible(dev.off())
    }
    
    #More calculations related to hydrology
    if(s2LoggerDeployed){
      outputDF$btwStaDist[i] <- inputFile[inputFile[[namLocIdx]] == S4 & inputFile[[eventIDIdx]] == currEventID, staDistIdx] -
        inputFile[inputFile[[namLocIdx]] == S2 & inputFile[[eventIDIdx]] == currEventID, staDistIdx] # meters
    }else{
      outputDF$btwStaDist[i] <- inputFile[inputFile[[namLocIdx]] == S4 & inputFile[[eventIDIdx]] == currEventID, staDistIdx] -
        inputFile[inputFile[[namLocIdx]] == S1 & inputFile[[eventIDIdx]] == currEventID, staDistIdx] # meters
    }
    outputDF$peakMaxVelocity[i] <- outputDF$btwStaDist[i]/as.numeric(outputDF$peakMaxTravelTime[i]) # m/s
    outputDF$meanQ_lps[i] <- mean(inputFile[inputFile[[eventIDIdx]] == currEventID, QIdx], na.rm = T) #lps
    outputDF$meanQ_cms[i] <- mean(inputFile[inputFile[[eventIDIdx]] == currEventID, QIdx], na.rm = T)*convLpsCms # m^3 s^-1
    outputDF$meanDepth[i] <- outputDF$meanQ_cms[i]/(inputFile[inputFile[[namLocIdx]] == S4 & inputFile[[eventIDIdx]] == currEventID, wwIdx]*outputDF$peakMaxVelocity[i]) # meters
    outputDF$meanTemp[i] <- inputFile[inputFile[[namLocIdx]] == S4 & inputFile[[eventIDIdx]] == currEventID, watTempIdx]
    
    # Carry over SF6 loss rate slopes from the inputFile
    try(outputDF$slopeRaw[i] <- unique(inputFile[inputFile[[eventIDIdx]] == currEventID, slopeRawIdx]))
    try(outputDF$slopeClean[i] <- unique(inputFile[inputFile[[eventIDIdx]] == currEventID, slopeCleanIdx]))
    try(outputDF$slopeSaltCorr[i] <- unique(inputFile[inputFile[[eventIDIdx]] == currEventID, slopeSaltCorrIdx]))
    try(outputDF$slopeBackCorr[i] <- unique(inputFile[inputFile[[eventIDIdx]] == currEventID, slopeBackCorrIdx]))
    
  }
  
  outputList <- list("outputDF"=outputDF,"inputFile"=inputFile)
  return(outputList)
}

