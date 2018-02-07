##############################################################################################
#' @title Reaeration rate and Schmidt number calculations

#' @author 
#' Kaelin M. Cawley \email{kcawley@battelleecology.org} \cr

#' @description This function calculates loss rate, travel time, SF6 reaeration rate, O2 gas transfer velocity, and Schmidt number 600.
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

#' @param inputFile Name of the data fram containing the information needed to calculate the dissolved gas concentrations. If the headers are named: "injectionType", "eventID", "stationToInjectionDistance", "plateauGasConc", "corrPlatSaltConc", "hoboSampleID", "wettedWidth", respectively, no other inputs are required. Otherwise, the names of the columns need to be input for the function to work.
#' @param dataDir User identifies the directory that contains the unzipped data
#' @param namedLocation A string identifier for the station where data was collected [string]
#' @param injectionType Either constant rate or slug [string]
#' @param eventID A string identifier to link records collected as part of the same experiment, SITE.YYYYMMDD for NEON [string]
#' @param stationToInjectionDistance Dataframe column name for distance from station to injection [string]
#' @param plateauGasConc Dataframe column name for natural log of gas concentration normalized to background corrected salt concentration [string]
#' @param corrPlatSaltConc Dataframe column name for natural log of gas concentration normalized to background corrected salt concentration [string]
#' @param hoboSampleID Dataframe column name for ID to link to conductivity timeseries data [string]
#' @param discharge Dataframe column name for stream discharge in literPerSecond [string]
#' @param waterTemp Dataframe column name for mean water temperature data [string]
#' @param wettedWidth Dataframe column name for mean wetted width for the stream reach [string]
#' @param plot User input to plot the SF6/corrected salt concentration versus distance downstream, defaults to FALSE [boolean]

#' @return This function returns a dataframe appended with loss rate, travel time, 
#' SF6 reaeration rate, O2 gas transfer velocity, and Schmidt number 600

#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

#' @keywords surface water, streams, rivers, reaeration, deaeration, SF6, metabolism, tracer 

#' @examples
#' #where the data frame "reaFormatted" is already read in
#' reaRatesCalc <- def.calc.reaeration(inputFile = reaFormatted, 
#' dataDir = paste(path.package("reaRate"),"inst\\extdata", sep = "\\"), plot = TRUE)
#' #where the data is read in from a file in the working directory (also works with a full path)
#' reaRatesCalc <- def.calc.reaeration(inputFile = 
#' system.file("extdata", "reaTestData.csv", package = "reaRate"))

#' @seealso def.calc.travelTime for calculating travel times and def.format.reaeration for formatting reaeration data

#' @export

# changelog and author contributions / copyrights
#   Kaelin M. Cawley (2017-08-03)
#     original creation
##############################################################################################
#This code is for calculating reaeration rates and Schmidt numbers
def.calc.reaeration <- function(
  inputFile,
  dataDir,
  namedLocation = "namedLocation",
  injectionType = "injectionType",
  eventID = "eventID",
  stationToInjectionDistance = "stationToInjectionDistance",
  plateauGasConc = "plateauGasConc",
  corrPlatSaltConc = "corrPlatSaltConc",
  hoboSampleID = "hoboSampleID",
  discharge = "discharge",
  waterTemp = "waterTemp",
  wettedWidth = "wettedWidth",
  plot = F
) {
  
  namLocIdx <- which(names(inputFile) == namedLocation)
  injTypeIdx <- which(names(inputFile) == injectionType)
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
  
  if(min(inputFile[[plGasIdx]], na.rm = T) < 0 | min(inputFile[[plSaltIdx]], na.rm = T) < 0){
    print("A gas concentration or background corrected salt concentration is negative producing NaNs for LNgasNormalizedToSalt")
  }
  try(inputFile$LNgasNormalizedTOSalt <- log(inputFile[[plGasIdx]]/inputFile[[plSaltIdx]]), silent = T)
  
  #Only use the unique eventIDs
  allEventID <- unique(inputFile[[eventIDIdx]])
  outputDF <- data.frame(matrix(data=NA, ncol=length(outputDFNames), nrow=length(allEventID)))
  names(outputDF) <- outputDFNames
  
  outputDF$eventID <- unique(inputFile[[eventIDIdx]])
  
  for(i in seq(along = outputDF$eventID)){
    #Uncomment this if you'd like to see a list of all the eventIDs for troubleshooting or debugging
    print(paste0(i, " - ", outputDF$eventID[i]))
    
    outputDF$siteID[i] <- unique(substr(inputFile[[namLocIdx]][inputFile[[eventIDIdx]] == outputDF$eventID[i]], 1, 4))
    S1 <- paste(outputDF$siteID[i], "AOS.reaeration.station.01", sep = ".")
    S4 <- paste(outputDF$siteID[i], "AOS.reaeration.station.04", sep = ".")
    
    #Calculate the Loss Rate, slope of the salt corrected SF6 over the reach
    x <- inputFile[inputFile[[eventIDIdx]] == outputDF$eventID[i], staDistIdx]
    y <- inputFile$LNgasNormalizedTOSalt[inputFile[[eventIDIdx]] == outputDF$eventID[i]]
    lineFit <- NA
    #Warnings when there isn't data suppressed
    suppressWarnings(try(lineFit <- lsfit(x,y), silent = T))
    
    if(sum(is.na(lineFit))){
      print(paste0("Warning, loss rate could not be determined for ", outputDF$eventID[i]))
      next
    }
    
    try(outputDF$lossRateSF6[i] <- lineFit$coefficients[[2]], silent = T)
    
    if(plot == T){
      invisible(dev.new(noRStudioGD = TRUE))
      plot(x,y,main = outputDF$eventID[i], xlab = "meters downstream of injection", ylab = "LN(Tracer Gas/Background Corrected Tracer Salt)", col = "blue")
      abline(lm(y ~ x, data = structure(list(x = x, y = y))))
      mtext(paste("y = ", lineFit$coefficients[[2]], "x +", lineFit$coefficients[[1]], "\n Click anywhere to close and continue"), cex = 0.8)
      #print("Click anywhere on the plot to close and continue")
      ans <- identify(x, y, n = 1, tolerance = 100, plot = F)
      invisible(dev.off())
    }
    
    #Calculate velocity to determine K_SF6
    outputDF$travelTime[i] <- def.calc.travelTime(dataDir = dataDir,
                                  currEventID = outputDF$eventID[i], 
                                  injectionType = unique(inputFile[inputFile[[eventIDIdx]] == outputDF$eventID[i], injTypeIdx]), 
                                  bPlot = plot) # seconds
    if(outputDF$travelTime[i] == -9999){
      outputDF$travelTime[i] = "Conductivity logger data not available, peak/plateau not found"
    }
    
    if(grepl("Conductivity logger data not available", outputDF$travelTime[i])){
      print(outputDF$travelTime[i])
      next
    }
    
    outputDF$btwStaDist[i] <- inputFile[inputFile[[namLocIdx]] == S4 & inputFile[[eventIDIdx]] == outputDF$eventID[i], staDistIdx] - 
      inputFile[inputFile[[namLocIdx]] == S1 & inputFile[[eventIDIdx]] == outputDF$eventID[i], staDistIdx] # meters
    outputDF$velocity[i] <- outputDF$btwStaDist[i]/as.numeric(outputDF$travelTime[i]) # m/s
    outputDF$reaRateSF6[i] <- outputDF$lossRateSF6[i] * outputDF$velocity[i] # m s^-2
    
    #Calculate the gas transfer velocity for oxygen
    outputDF$reaRateO2[i] <- outputDF$reaRateSF6[i] * reaRateConv #convert from SF6 to O2 reaeration rate coefficient
    
    #Determine gas transfer velocity for O2
    outputDF$meanQ[i] <- mean(inputFile[inputFile[[eventIDIdx]] == outputDF$eventID[i], QIdx], na.rm = T)*convLpsCms # m^3 s^-1
    outputDF$meanDepth[i] <- outputDF$meanQ[i]/inputFile[inputFile[[namLocIdx]] == S4 & inputFile[[eventIDIdx]] == outputDF$eventID[i], wwIdx]/outputDF$velocity[i] # meters
    outputDF$gasTransVelO2[i] <- outputDF$reaRateO2[i] * outputDF$meanDepth[i] # m^2 s^-2
    
    #Normalize to schmidt number of 600
    outputDF$meanTemp[i] <- inputFile[inputFile[[namLocIdx]] == S4 & inputFile[[eventIDIdx]] == outputDF$eventID[i], watTempIdx]
    scO2 <- A_O2 - B_O2 * outputDF$meanTemp[i] + C_O2 * outputDF$meanTemp[i]^2 - D_O2 * outputDF$meanTemp[i]^3
    outputDF$k600[i] <- (Sc_CO2/scO2)^(-0.5) * outputDF$gasTransVelO2[i] #Equation 1, Wanninkhof (1992)
    
  }
  
  return(outputDF)
  
}

