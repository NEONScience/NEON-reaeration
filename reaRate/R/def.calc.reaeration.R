##############################################################################################
#' @title Reaeration rate and Schmidt number calculations

#' @author
#' Kaelin M. Cawley \email{kcawley@battelleecology.org} \cr

#' @description This function calculates loss rate, travel time, SF6 reaeration rate, O2
#' gas transfer velocity, and Schmidt number 600.
#' @importFrom grDevices dev.new
#' @importFrom grDevices dev.off
#' @importFrom grDevices dev.copy
#' @importFrom grDevices png
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

#' @param inputFile Name of the data fram containing the information needed to calculate the
#' reaeration parameters. If the headers are named: "injectionType", "eventID",
#' "stationToInjectionDistance", "plateauGasConc", "corrPlatSaltConc", "hoboSampleID",
#' "wettedWidth", respectively, no other inputs are required. Otherwise, the names of the
#' columns need to be input for the function to work. [string]
#' @param loggerData User identified filename of logger data [string]
#' @param namedLocation A string identifier for the station where data was collected [string]
#' @param injectionTypeName Either constant rate or slug [string]
#' @param eventID A string identifier to link records collected as part of the same experiment,
#' SITE.YYYYMMDD for NEON [string]
#' @param stationToInjectionDistance Dataframe column name for distance from station to
#' injection [string]
#' @param plateauGasConc Dataframe column name for natural log of gas concentration normalized to
#' background corrected salt concentration [string]
#' @param corrPlatSaltConc Dataframe column name for natural log of gas concentration normalized to
#' background corrected salt concentration [string]
#' @param hoboSampleID Dataframe column name for ID to link to conductivity timeseries data [string]
#' @param discharge Dataframe column name for stream discharge in literPerSecond [string]
#' @param waterTemp Dataframe column name for mean water temperature data [string]
#' @param wettedWidth Dataframe column name for mean wetted width for the stream reach [string]
#' @param plot User input to plot the SF6/corrected salt concentration versus distance downstream,
#' defaults to TRUE [boolean]
#' @param savePlotPath If a user specifies a path the plots will be saved to this location [string]
#' @param processingInfo Metadata about choices made while processing data [dataframe]

#' @return This function returns a list of two dataframes, the input dataframe of data for up to
#' 4 stations per site per date and an output dataframe appended with loss rate, travel time,
#' SF6 reaeration rate, O2 gas transfer velocity, and Schmidt number 600 for a given site and date

#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

#' @keywords surface water, streams, rivers, reaeration, deaeration, SF6, metabolism, tracer

#' @examples
#' #where the data frame "reaFormatted" is already read in
#' #reaRatesCalc <- def.calc.reaeration(inputFile = reaFormatted,
#' #dataDir = paste(path.package("reaRate"),"inst\\extdata", sep = "\\"), plot = TRUE)
#' #where the data is read in from a file in the working directory (also works with a full path)
#' #reaRatesCalc <- def.calc.reaeration(inputFile =
#' #system.file("extdata", "reaTestData.csv", package = "reaRate"))

#' @seealso def.calc.peakTime for calculating travel times and def.format.reaeration for
#' formatting reaeration data

#' @export

# changelog and author contributions / copyrights
#   Kaelin M. Cawley (2017-08-03)
#     original creation
#   Kaelin M. Cawley (2018-05-03)
#     added functionality for saving plots to a specified directory
#   Kaelin M. Cawley (2020-12-21)
#     updated to only plot after slugPourTime or dripStartTime
#   Kaelin M. Cawley (2021-05-25)
#     updated to calculate supset of metrics for model injection types
##############################################################################################
#This code is for calculating reaeration rates and Schmidt numbers
def.calc.reaeration <- function(
  inputFile = NULL,
  loggerData = NULL,
  namedLocation = "namedLocation",
  injectionTypeName = "injectionType",
  eventID = "eventID",
  stationToInjectionDistance = "stationToInjectionDistance",
  plateauGasConc = "plateauGasConc",
  corrPlatSaltConc = "corrPlatSaltConc",
  hoboSampleID = "hoboSampleID",
  discharge = "fieldDischarge",
  waterTemp = "waterTemp",
  wettedWidth = "wettedWidth",
  plot = TRUE,
  savePlotPath = NULL,
  processingInfo = NULL
){

  if(!plot && !is.null(savePlotPath)){
    stop("Please turn plotting on (plot = T) in order to save plots.")
  }

  namLocIdx <- which(names(inputFile) == namedLocation)
  injTypeIdx <- which(names(inputFile) == injectionTypeName)
  eventIDIdx <- which(names(inputFile) == eventID)
  staDistIdx <- which(names(inputFile) == stationToInjectionDistance)
  plGasIdx <- which(names(inputFile) == plateauGasConc)
  plSaltIdx <- which(names(inputFile) == corrPlatSaltConc)
  loggerIdx <- which(names(inputFile) == hoboSampleID)
  QIdx <- which(names(inputFile) == discharge)
  watTempIdx <- which(names(inputFile) == waterTemp)
  wwIdx <- which(names(inputFile) == wettedWidth)

  ##### Constants #####
  #Coefficients for Least Squares Third-Order Polynomial Fits of Schmidt Number Versus Temperature
  #Valid for 0 - 30 Celsius temperature range
  #Table A1, Fresh water, Wanninkhof (1992), DOI: 10.1029/92JC00188
  #See also Jahne et al. (1987), DOI: 10.4236/jep.2014.511103
  A_O2 = 1800.6
  B_O2 = 120.10
  C_O2 = 3.7818
  D_O2 = 0.047608

  A_CO2 = 1911.1
  B_CO2 = 118.11
  C_CO2 = 3.4527
  D_CO2 = 0.041320

  A_SF6 = 3255.3
  B_SF6 = 217.13
  C_SF6 = 6.8370
  D_SF6 = 0.086070

  Sc_CO2 = 600 #Schmidt number of O2 at 20 C in fresh water

  convLpsCms = 1/1000 #Conversion from litersPerSecond to cubicMetersPerSecond

  #Reaeration Rate Conversion
  #Equation 7, Wanninkhof (1990), DOI: 10.1029/WR026i007p01621
  Sc_O2_25 <- A_O2 - B_O2 * 25 + C_O2 * 25^2 - D_O2 * 25^3
  Sc_SF6_25 <- A_SF6 - B_SF6 * 25 + C_SF6 * 25^2 - D_SF6 * 25^3
  reaRateConv <- (Sc_O2_25/Sc_SF6_25) ^ (-0.5)

  #Create output file
  outputDFNames <- c(
    'siteID',
    'startDate',
    'eventID',
    'lossRateSF6',
    'S1PeakTime',
    'S4PeakTime',
    'travelTime',
    'btwStaDist',
    'velocity',
    'reaRateSF6',
    'meanDepth',
    'reaRateO2',
    'gasTransVelO2',
    'meanQ',
    'meanTemp',
    'k600'
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

    outputDF$siteID[i] <- unique(substr(inputFile[[namLocIdx]][inputFile[[eventIDIdx]] == currEventID], 1, 4))
    S1 <- paste(outputDF$siteID[i], "AOS.reaeration.station.01", sep = ".")
    S2 <- paste(outputDF$siteID[i], "AOS.reaeration.station.02", sep = ".")
    S3 <- paste(outputDF$siteID[i], "AOS.reaeration.station.03", sep = ".")
    S4 <- paste(outputDF$siteID[i], "AOS.reaeration.station.04", sep = ".")

    if(!modelInjType){

      #Background correct salt samples, normalize gas concentration, and natural log transform the plateau gas concentrations
      backSalt <- inputFile$backgroundSaltConc[inputFile$eventID == currEventID]
      platSalt <- as.character(inputFile$plateauSaltConc[inputFile$eventID == currEventID])
      platGas <- as.character(inputFile$plateauGasConc[inputFile$eventID == currEventID])
      statDist <- inputFile$stationToInjectionDistance[inputFile$eventID == currEventID]

      #If the background values are below detection, just use 0
      if(any(is.na(backSalt))){
        backSalt[is.na(backSalt)] <- 0
      }

      x <- NA
      y <- NA
      meanY <- NA
      for(j in 1:length(statDist)){
        currStart <- (j-1)*5

        currBack <- backSalt[j]
        currPlatSalt <- as.numeric(strsplit(platSalt[j],"\\|")[[1]])
        currPlatGas <- as.numeric(strsplit(platGas[j],"\\|")[[1]])

        #Background correct plateau salt concentrations
        corrPlatSalt <- NA
        if(length(currPlatSalt)>0 && length(currBack)>0){
          corrPlatSalt <- currPlatSalt-currBack
        }

        #Normalize plateau gas concentration to corrected plateau salt concentration
        normPlatGas <- NA
        if(length(currPlatGas)>0 && length(corrPlatSalt)>0 && any(!is.na(currPlatGas)) && any(!is.na(corrPlatSalt)) && length(currPlatGas)==length(corrPlatSalt)){
          normPlatGas <- currPlatGas/corrPlatSalt
        }

        if(length(currPlatSalt)<1 || length(currBack)<1 || length(currPlatGas)<1 || all(is.na(currPlatGas)) || all(is.na(currPlatGas))){
          print(paste0("Tracer data for station ",j,", eventID ",currEventID," not available."))
          next
        }

        if(min(normPlatGas, na.rm = T) <= 0 | min(corrPlatSalt, na.rm = T) <= 0){
          print("A gas concentration or background corrected salt concentration is zero or negative producing NaNs for LNgasNormalizedToSalt")
        }

        normPlatGas[normPlatGas <= 0] <- NA
        corrPlatSalt[corrPlatSalt <= 0] <- NA

        logNormPlatGas <- try(log(normPlatGas))

        numVals <- min(length(corrPlatSalt),length(normPlatGas))

        x[(1+currStart):(numVals+currStart)] <- statDist[j]
        y[(1+currStart):(numVals+currStart)] <- logNormPlatGas
        meanY[j] <- log(mean(currPlatGas, na.rm = T)/mean(corrPlatSalt, na.rm = T))
      }

      #Calculate the Loss Rate, slope of the salt corrected SF6 over the reach
      lineFit <- NA
      #Warnings when there isn't data suppressed
      suppressWarnings(try(lineFit <- lsfit(statDist,meanY), silent = T))

      if(sum(is.na(lineFit))){
        print(paste0("Warning, loss rate could not be determined for ", currEventID))
        next
      }

      #Clean up y for plotting if there are Inf values
      x <- x[!is.infinite(y)]
      y <- y[!is.infinite(y)]

      try(outputDF$lossRateSF6[i] <- lineFit$coefficients[[2]], silent = T)

      if(plot == T & !all(is.na(x)) & !all(is.na(y))){
        #Save out plot of loss rate to specified directory
        if(!is.null(savePlotPath)){
          png(paste0(savePlotPath,"/lossRate_",currEventID,".png"))
          plot(x,y,main = currEventID, xlab = "meters downstream of injection", ylab = "LN(Tracer Gas/Background Corrected Tracer Salt)", col = "blue")
          points(statDist,meanY, pch = 19)
          abline(a = lineFit$coefficients[["Intercept"]], b = lineFit$coefficients[["X"]])
          mtext(paste("y = ", lineFit$coefficients[[2]], "x +", lineFit$coefficients[[1]], "\n Click anywhere to close and continue"), cex = 0.8)
          dev.off()
        }

        invisible(dev.new(noRStudioGD = TRUE))
        plot(x,y,main = currEventID, xlab = "meters downstream of injection", ylab = "LN(Tracer Gas/Background Corrected Tracer Salt)", col = "blue")
        points(statDist,meanY, pch=19)
        abline(a = lineFit$coefficients[["Intercept"]], b = lineFit$coefficients[["X"]])
        mtext(paste("y = ", lineFit$coefficients[[2]], "x +", lineFit$coefficients[[1]], "\n Click anywhere to close and continue"), cex = 0.8)
        #print("Click anywhere on the plot to close and continue")
        ans <- identify(x, y, n = 1, tolerance = 100, plot = F)

        invisible(dev.off())
      }
    }

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
      condDataS1 <- s1LoggerData$lowRangeSpCondNonlinear
      s1RangeFull <- FALSE
    }else if(!all(is.na(s1LoggerData$fullRangeSpCondNonlinear))){
      condDataS1 <- s1LoggerData$fullRangeSpCondNonlinear
      s1RangeFull <- TRUE
    }else{
      print(paste0("Conductivity logger data not available for ", currEventID, ", station S1"))
      next
    }

    if(!all(is.na(s4LoggerData$lowRangeSpCondNonlinear))){
      condDataS4 <- s4LoggerData$lowRangeSpCondNonlinear
      s4RangeFull <- FALSE
    }else if(!all(is.na(s4LoggerData$fullRangeSpCondNonlinear))){
      condDataS4 <- s4LoggerData$fullRangeSpCondNonlinear
      s4RangeFull <- TRUE
    }else{
      print(paste0("Conductivity logger data not available for ", currEventID, ", station S4"))
      next
    }

    #Find the peak locations
    s1peakLoc <- reaRate::def.calc.peakTime(loggerData = condDataS1,
                                   currEventID = currEventID,
                                   injectionType = injectionType) # index to get date and time of peak/plateau half max
    s4peakLoc <- reaRate::def.calc.peakTime(loggerData = condDataS4,
                                   currEventID = currEventID,
                                   injectionType = injectionType) # index to get date and time of peak/plateau half max

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
    outputDF$S1PeakTime[i] <- s1LoggerData$dateTimeLogger[s1peakLoc$peakLocOut]
    outputDF$S4PeakTime[i] <- s4LoggerData$dateTimeLogger[s4peakLoc$peakLocOut]
    outputDF$travelTime[i] <- difftime(s4LoggerData$dateTimeLogger[s4peakLoc$peakLocOut],
                                       s1LoggerData$dateTimeLogger[s1peakLoc$peakLocOut],
                                       units = "secs")

    #Plot the travel times to check
    if(plot==T){
      if(s1RangeFull){
        s1YData <- s1LoggerData$fullRangeSpCondNonlinear[s1peakLoc$peakStart:s1peakLoc$peakEnd]
      }else{
        s1YData <- s1LoggerData$lowRangeHobo[s1peakLoc$peakStart:s1peakLoc$peakEnd]
      }

      if(s4RangeFull){
        s4YData <- s4LoggerData$fullRangeSpCondNonlinear[s4peakLoc$peakStart:s4peakLoc$peakEnd]
      }else{
        s4YData <- s4LoggerData$lowRangeHobo[s4peakLoc$peakStart:s4peakLoc$peakEnd]
      }

      invisible(dev.new(noRStudioGD = TRUE))
      x <- as.POSIXct(s1LoggerData$dateTimeLogger[s1peakLoc$peakStart:s1peakLoc$peakEnd],format = dateFormat)
      #y <- s1LoggerData$fullRangeSpCondNonlinear[s1peakLoc$peakStart:s1peakLoc$peakEnd]
      minTime <- min(as.POSIXct(s1LoggerData$dateTimeLogger[s1peakLoc$peakStart:s1peakLoc$peakEnd],format = dateFormat),
                     as.POSIXct(s4LoggerData$dateTimeLogger[s4peakLoc$peakStart:s4peakLoc$peakEnd],format = dateFormat))
      maxTime <- max(as.POSIXct(s1LoggerData$dateTimeLogger[s1peakLoc$peakStart:s1peakLoc$peakEnd],format = dateFormat),
                     as.POSIXct(s4LoggerData$dateTimeLogger[s4peakLoc$peakStart:s4peakLoc$peakEnd],format = dateFormat))
      minY <- min(s1YData,s4YData,na.rm = T)
      maxY <- max(s1YData,s4YData,na.rm = T)

      #Save out plot of loss rate to specified directory
      if(!is.null(savePlotPath)){
        png(paste0(savePlotPath,"/travelTime_",currEventID,".png"))
        plot(x,
             s1YData,
             xlim = c(minTime,maxTime),
             ylim = c(minY,maxY),
             ylab = "Conductivity, uS",
             xlab = "Time (UTC)")
        mtext(paste0("Travel Time = ",outputDF$travelTime[i]," seconds, (",round(as.numeric(outputDF$travelTime[i])/60,digits=1) ," min)\n Click anywhere to close and continue"), cex = 1.2)
        points(as.POSIXct(s4LoggerData$dateTimeLogger[s4peakLoc$peakStart:s4peakLoc$peakEnd],format = dateFormat),
               s4YData,
               col = "blue")
        abline(v = as.POSIXct(s1LoggerData$dateTimeLogger[s1peakLoc$peakLocOut],format = dateFormat))
        abline(v = as.POSIXct(s4LoggerData$dateTimeLogger[s4peakLoc$peakLocOut],format = dateFormat), col = "blue")
        dev.off()
      }
      plot(x,
           s1YData,
           xlim = c(minTime,maxTime),
           ylim = c(minY,maxY),
           ylab = "Conductivity, uS",
           xlab = "Time (UTC)")
      mtext(paste0("Travel Time = ",outputDF$travelTime[i]," seconds, (",round(as.numeric(outputDF$travelTime[i])/60,digits=1) ," min)\n Click anywhere to close and continue"), cex = 1.2)
      points(as.POSIXct(s4LoggerData$dateTimeLogger[s4peakLoc$peakStart:s4peakLoc$peakEnd],format = dateFormat),
             s4YData,
             col = "blue")
      abline(v = as.POSIXct(s1LoggerData$dateTimeLogger[s1peakLoc$peakLocOut],format = dateFormat))
      abline(v = as.POSIXct(s4LoggerData$dateTimeLogger[s4peakLoc$peakLocOut],format = dateFormat), col = "blue")
      ans <- identify(x, s1YData, n = 1, tolerance = 100, plot = F)
      invisible(dev.off())
    }

    #More calculations to get to the reaeration rate
    if(s2LoggerDeployed){
      outputDF$btwStaDist[i] <- inputFile[inputFile[[namLocIdx]] == S4 & inputFile[[eventIDIdx]] == currEventID, staDistIdx] -
        inputFile[inputFile[[namLocIdx]] == S2 & inputFile[[eventIDIdx]] == currEventID, staDistIdx] # meters
    }else{
      outputDF$btwStaDist[i] <- inputFile[inputFile[[namLocIdx]] == S4 & inputFile[[eventIDIdx]] == currEventID, staDistIdx] -
        inputFile[inputFile[[namLocIdx]] == S1 & inputFile[[eventIDIdx]] == currEventID, staDistIdx] # meters
    }
    outputDF$velocity[i] <- outputDF$btwStaDist[i]/as.numeric(outputDF$travelTime[i]) # m/s
    outputDF$reaRateSF6[i] <- outputDF$lossRateSF6[i] * outputDF$velocity[i] * -1 * 86400# m^-1 * m/s * -1 for negative slope and 86400 for number of seconds in a day

    #Calculate the gas transfer velocity for oxygen
    outputDF$reaRateO2[i] <- outputDF$reaRateSF6[i] * reaRateConv #convert from SF6 to O2 reaeration rate coefficient

    #Determine gas transfer velocity for O2
    outputDF$meanQ[i] <- mean(inputFile[inputFile[[eventIDIdx]] == currEventID, QIdx], na.rm = T)*convLpsCms # m^3 s^-1
    outputDF$meanDepth[i] <- outputDF$meanQ[i]/(inputFile[inputFile[[namLocIdx]] == S4 & inputFile[[eventIDIdx]] == currEventID, wwIdx]*outputDF$velocity[i]) # meters
    outputDF$gasTransVelO2[i] <- outputDF$reaRateO2[i] * outputDF$meanDepth[i] # d^-1 * m

    #Normalize to schmidt number of 600
    outputDF$meanTemp[i] <- inputFile[inputFile[[namLocIdx]] == S4 & inputFile[[eventIDIdx]] == currEventID, watTempIdx]
    scO2 <- A_O2 - B_O2 * outputDF$meanTemp[i] + C_O2 * outputDF$meanTemp[i]^2 - D_O2 * outputDF$meanTemp[i]^3
    outputDF$k600[i] <- (Sc_CO2/scO2)^(-0.5) * outputDF$gasTransVelO2[i] #Equation 1, Wanninkhof (1992)

  }
  outputList <- list("outputDF"=outputDF,"inputFile"=inputFile)
  return(outputList)
}

