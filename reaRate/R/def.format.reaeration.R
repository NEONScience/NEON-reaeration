##############################################################################################
#' @title Formats reaeration data for rate calculations

#' @author 
#' Kaelin M. Cawley \email{kcawley@battelleecology.org} \cr

#' @description This function reads in data from the NEON reaeration data product to calculate 
#' loss rate, travel time, SF6 reaeration rate, O2 gas transfer velocity, and Schmidt number 600. 
#' Either the basic or expanded package can be downloaded. No need to unzip the downloaded files, 
#' just place them all in the same directory.
#' @importFrom streamQ def.format.Q
#' @importFrom streamQ def.calc.Q.inj
#' @importFrom streamQ def.calc.Q.slug
#' @importFrom neonUtilities stackByTable
#' @importFrom utils read.csv

#' @param dataDir User identifies the directory that contains the zipped data

#' @return This function returns one data frame formatted for use with def.calc.reaeration.R

#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

#' @keywords surface water, streams, rivers, reaeration, gas transfer velocity, schmidt number

#' @examples
#' #where the data .zip file is in the working directory and has the default name, 
#' #reaFormatted <- def.format.reaeration()
#' #where the data.zip file is in the downloads folder and has default name, 
#' #reaFormatted <- 
#' #def.format.reaeration(dataDir = path.expand("~/Downloads/NEON_reaeration.zip"))
#' #where the data.zip file is in the downloads folder and has a specified name,
#' #reaFormatted <- def.format.reaeration(dataDir = path.expand("~/Downloads/non-standard-name.zip"))
#' #Using the example data in this package
#' dataDirectory <- paste(path.package("reaRate"),"inst\\extdata", sep = "\\")
#' reaFormatted <- def.format.reaeration(dataDir = dataDirectory)

#' @seealso def.calc.tracerTime.R for calculating the stream travel time, 
#' def.plot.reaQcurve.R for plotting reaeration rate versusu stream flow

#' @export

# changelog and author contributions / copyrights
#   Kaelin M. Cawley (2017-10-30)
#     original creation
##############################################################################################
#This code is for calculating salt-based discharge for a slug
def.format.reaeration <- function(
  dataDir = paste0(getwd(),"/NEON_reaeration.zip")
) {
  
  #Stack field and external lab data
  if(!dir.exists(substr(dataDir, 1, (nchar(dataDir)-4)))){
    stackByTable(dpID="DP1.20190.001",filepath=dataDir)
  }
  
  #Read in stacked data
  backgroundDataLogger <- read.csv(
    paste(gsub("\\.zip","",dataDir), "stackedFiles", "rea_backgroundFieldCondData.csv", sep = "/"), 
    stringsAsFactors = F)
  
  backgroundDataSalt <- read.csv(
    paste(gsub("\\.zip","",dataDir), "stackedFiles", "rea_backgroundFieldSaltData.csv", sep = "/"), 
    stringsAsFactors = F)
  
  fieldDataSite <- read.csv(
    paste(gsub("\\.zip","",dataDir), "stackedFiles", "rea_fieldData.csv", sep = "/"), 
    stringsAsFactors = F)
  fieldDataSite$namedLocation <- NULL #So that merge goes smoothly
  
  plateauDataCond <- read.csv(
    paste(gsub("\\.zip","",dataDir),"stackedFiles","rea_plateauMeasurementFieldData.csv", sep = "/"), 
    stringsAsFactors = F)
  
  plateauDataSalt <- read.csv(
    paste(gsub("\\.zip","",dataDir),"stackedFiles","rea_plateauSampleFieldData.csv", sep = "/"), 
    stringsAsFactors = F)
  
  externalLabDataSalt <- read.csv(
    paste(gsub("\\.zip","",dataDir),"stackedFiles","rea_externalLabDataSalt.csv", sep = "/"), 
    stringsAsFactors = F)
  
  externalLabDataGas <- read.csv(
    paste(gsub("\\.zip","",dataDir),"stackedFiles","rea_externalLabDataGas.csv", sep = "/"), 
    stringsAsFactors = F)
  
  wettedWidths <- read.csv(
    paste(gsub("\\.zip","",dataDir),"stackedFiles","rea_widthFieldData.csv", sep = "/"), 
    stringsAsFactors = F)
  
  loggerData <- read.csv(
    paste(gsub("\\.zip","",dataDir),"stackedFiles","rea_conductivityFieldData.csv", sep = "/"), 
    stringsAsFactors = F)
  
  #Merge the backgroundDataSalt and fieldDataSite tables
  loggerSiteData <- merge(backgroundDataSalt, 
                          fieldDataSite, 
                          by = c('siteID', 'startDate'), 
                          all = T)
  
  #Create input file for reaeration calculations
  outputDFNames <- c(
    'siteID',
    'namedLocation', #Station at this point
    'startDate',
    'stationToInjectionDistance',
    'injectionType',
    'slugPourTime',
    'dripStartTime',
    'backgroundSaltConc',
    'plateauSaltConc',
    'corrPlatSaltConc',
    'plateauGasConc',
    'wettedWidth',
    'waterTemp',
    'hoboSampleID',
    'discharge',
    'eventID'
  )
  outputDF <- data.frame(matrix(data=NA, ncol=length(outputDFNames), nrow=length(loggerSiteData$siteID)))
  names(outputDF) <- outputDFNames
  
  #Fill in the fields from the loggerSiteData table
  for(i in seq(along = names(outputDF))){
    if(names(outputDF)[i] %in% names(loggerSiteData)){
      outputDF[,i] <- loggerSiteData[,which(names(loggerSiteData) == names(outputDF)[i])]
    }
  }
  
  #Change to more generic injection types
  outputDF$injectionType[outputDF$injectionType == "NaCl"] <- "constant"
  outputDF$injectionType[outputDF$injectionType == "NaBr"|outputDF$injectionType == "model"] <- "slug"
  
  outputDF$eventID <- paste0(outputDF$siteID, ".", substr(outputDF$startDate,1,4), substr(outputDF$startDate,6,7), substr(outputDF$startDate,9,10))
  
  QFile <- def.format.Q(dataDir)
  QFile <- def.calc.Q.inj(QFile)
  QFile <- def.calc.Q.slug(inputFile = QFile, dataDir = dataDir)
  
  for(i in seq(along = outputDF$siteID)){
    
    siteID <- outputDF$siteID[i]
    startDate <- outputDF$startDate[i]
    station <- outputDF$namedLocation[i]
    stationType <- substr(station, 6, nchar(station))
    
    repRegex <- switch(stationType,
                       "AOS.reaeration.station.01" = "\\.0[12345]\\.",
                       "AOS.reaeration.station.02" = "\\.0[6789]\\.|\\.10\\.",
                       "AOS.reaeration.station.03" = "\\.1[12345]\\.",
                       "AOS.reaeration.station.04" = "\\.1[6789]\\.|\\.20\\."
    )
    
    if(!is.na(outputDF$injectionType[i]) && outputDF$injectionType[i] == "constant"){
      try(outputDF$discharge[i] <- QFile$Q.inj[QFile$namedLocation == station
                                           & QFile$startDate == startDate], silent = T)
    }else if(!is.na(outputDF$injectionType[i]) && 
             (outputDF$injectionType[i] == "slug")){
      try(outputDF$discharge[i] <- QFile$Q.slug[QFile$namedLocation == station
                                            & QFile$startDate == startDate], silent = T)
    }
    
    #Fill in hoboSampleID from background logger table
    try(outputDF$hoboSampleID[i] <- backgroundDataLogger$hoboSampleID[
      backgroundDataLogger$namedLocation == station &
        backgroundDataLogger$startDate == startDate], silent = T)
    
    #Fill in background concentration data
    try(outputDF$backgroundSaltConc[i] <- 
          unique(externalLabDataSalt$finalConcentration[
            externalLabDataSalt$namedLocation == station &
              externalLabDataSalt$startDate == startDate & 
              grepl(paste0(".B",substr(station,nchar(station),nchar(station)),"."),
                    externalLabDataSalt$saltSampleID)]), silent = T)
    
    #Fill in plateau concentration data for constant rate injection
    pSaltConc <- externalLabDataSalt$finalConcentration[
      externalLabDataSalt$namedLocation == station &
        externalLabDataSalt$startDate == startDate & 
        grepl(repRegex, externalLabDataSalt$saltSampleID)]
    
    #Remove outliers TBD
    #Calculate the mean plateau concentration
    outputDF$plateauSaltConc[i] <- ifelse(!is.nan(mean(pSaltConc, na.rm = T)),mean(pSaltConc, na.rm = T),NA)
    
    #Fill in plateau gas concentration
    pGasConc <- externalLabDataGas$gasTracerConcentration[
      externalLabDataGas$namedLocation == station &
        externalLabDataGas$startDate == startDate & 
        grepl(repRegex, externalLabDataGas$gasSampleID)]
    
    #Remove outliers TBD
    #Calculate the mean plateau concentration
    outputDF$plateauGasConc[i] <- ifelse(!is.nan(mean(pGasConc, na.rm = T)),mean(pGasConc, na.rm = T),NA)
    
    #Fill in mean wetted width
    wettedWidthVals <- wettedWidths$wettedWidth[
      wettedWidths$namedLocation == siteID &
        grepl(substr(startDate, 1, 10), wettedWidths$startDate)]
    
    #Remove outliers TBD
    #Calculate the mean wetted width
    outputDF$wettedWidth[i] <- ifelse(!is.nan(mean(wettedWidthVals, na.rm = T)),mean(wettedWidthVals, na.rm = T),NA)
    
    #Fill in mean water temp
    tempVals <- plateauDataCond$waterTemp[
      plateauDataCond$siteID == siteID &
        grepl(substr(startDate, 1, 10), plateauDataCond$startDate)]
    
    #Remove outliers TBD
    #Calculate the mean water temp
    outputDF$waterTemp[i] <- ifelse(!is.nan(mean(tempVals, na.rm = T)),mean(tempVals, na.rm = T),NA)
    
  }
  
  outputDF$corrPlatSaltConc[!is.na(outputDF$plateauSaltConc) & !is.na(outputDF$backgroundSaltConc)] <- 
    as.numeric(outputDF$plateauSaltConc[!is.na(outputDF$plateauSaltConc) & !is.na(outputDF$backgroundSaltConc)]) - 
    as.numeric(outputDF$backgroundSaltConc[!is.na(outputDF$plateauSaltConc) & !is.na(outputDF$backgroundSaltConc)])
  
  return(outputDF)

  
}

