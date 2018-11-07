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
#' @importFrom neonUtilities zipsByProduct
#' @importFrom utils read.csv

#' @param dataDir User identifies the directory that contains the zipped data or sets to 
#' "API" to pull data from the NEON API [string]
#' @param site User identifies the site(s), defaults to "all" [string]
#' @param fieldQ specifies whether or no field discharge data should be included [boolean]

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
#' #dataDirectory <- paste(path.package("reaRate"),"inst\\extdata", sep = "\\")
#' #reaFormatted <- def.format.reaeration(dataDir = dataDirectory)

#' @seealso def.calc.tracerTime.R for calculating the stream travel time, 
#' def.plot.reaQcurve.R for plotting reaeration rate versusu stream flow

#' @export

# changelog and author contributions / copyrights
#   Kaelin M. Cawley (2017-10-30)
#     original creation
#   Kaelin M. Cawley (2018-05-03)
#     added the option of getting data from the API rather than a file download
##############################################################################################
def.format.reaeration <- function(
  dataDir = paste0(getwd(),"/NEON_reaeration.zip"),
  site = "all",
  fieldQ = FALSE
) {
  
  reaDPID <- "DP1.20190.001"
  qDPID <- "DP1.20048.001"
  folder <- FALSE
  #Pull files from the API to stack
  if(dataDir == "API"&&!dir.exists(paste(getwd(), "/filesToStack", substr(reaDPID, 5, 9), sep=""))){
    dataFromAPI <- zipsByProduct(reaDPID,site,package="expanded",check.size=TRUE)
    if(fieldQ){
      fieldQAPI <- zipsByProduct(qDPID,site,package="basic",check.size=TRUE)
    }
  }
  
  if(dataDir == "API"){
    filepath <- paste(getwd(), "/filesToStack", substr(reaDPID, 5, 9), sep="")
    qFilepath <- paste(getwd(), "/filesToStack", substr(qDPID, 5, 9), sep="")
    folder <- TRUE
  } else{
    filepath = dataDir
    qFilepath = dataDir
  }
  
  #Stack field and external lab data
  if(!dir.exists(paste(gsub("\\.zip","",filepath), "/stackedFiles", sep = "/"))&&
     file.exists(filepath)){
    stackByTable(dpID=reaDPID,filepath=filepath,folder=folder)
    filepath <- paste(gsub("\\.zip","",filepath), "stackedFiles", sep = "/")
  }
  
  #Stack discharge files if needed
  if(!dir.exists(paste(gsub("\\.zip","",qFilepath), "/stackedFiles", sep = "/"))&&
     file.exists(qFilepath)&&
     fieldQ){
    stackByTable(dpID=qDPID,filepath=qFilepath,folder=TRUE)
    qFilepath <- paste(gsub("\\.zip","",qFilepath), "stackedFiles", sep = "/")
  }else if(dir.exists(paste(gsub("\\.zip","",filepath), "/stackedFiles", sep = "/"))){
    filepath <- paste(gsub("\\.zip","",filepath), "stackedFiles", sep = "/")
    if(fieldQ){
      qFilepath <- paste(gsub("\\.zip","",qFilepath), "stackedFiles", sep = "/")
    }
  }

  #Read in stacked files
  if(dir.exists(filepath)){
    #Read in stacked data
    rea_backgroundFieldCondData <- read.csv(
      paste(filepath,"rea_backgroundFieldCondData.csv", sep = "/"), 
      stringsAsFactors = F)
    
    try(rea_backgroundFieldSaltData <- read.csv(
      paste(filepath, "rea_backgroundFieldSaltData.csv", sep = "/"), 
      stringsAsFactors = F))
    
    rea_fieldData <- read.csv(
      paste(filepath,"rea_fieldData.csv", sep = "/"), 
      stringsAsFactors = F)
    
    try(rea_plateauMeasurementFieldData <- read.csv(
      paste(filepath,"rea_plateauMeasurementFieldData.csv", sep = "/"), 
      stringsAsFactors = F))
    
    # #This isn't used anywhere else
    # rea_plateauSampleFieldData <- read.csv(
    #   paste(filepath,"rea_plateauSampleFieldData.csv", sep = "/"), 
    #   stringsAsFactors = F)
    
    try(rea_externalLabDataSalt <- read.csv(
      paste(filepath,"rea_externalLabDataSalt.csv", sep = "/"), 
      stringsAsFactors = F))
    
    try(rea_externalLabDataGas <- read.csv(
      paste(filepath,"rea_externalLabDataGas.csv", sep = "/"), 
      stringsAsFactors = F))
    
    rea_widthFieldData <- read.csv(
      paste(filepath,"rea_widthFieldData.csv", sep = "/"), 
      stringsAsFactors = F)
    
    # #This isn't used anywhere else
    # rea_conductivityFieldData <- read.csv(
    #   paste(filepath,"rea_conductivityFieldData.csv", sep = "/"), 
    #   stringsAsFactors = F)
  } else{
    stop("Error, stacked files could not be read in reaeration data")
  }
  
  #Read in stacked field discharge data
  if(fieldQ&&dir.exists(qFilepath)){
    #Read in stacked data
    dsc_fieldData <- read.csv(
      paste(qFilepath,"dsc_fieldData.csv", sep = "/"), 
      stringsAsFactors = F)
    dsc_fieldData$eventID <- paste(dsc_fieldData$siteID,gsub("-","",substr(dsc_fieldData$startDate,1,10)),sep = ".")
  } else{
    stop("Error, stacked discharge files could not be read in reaeration data")
  }
  
  rea_fieldData$namedLocation <- NULL #So that merge goes smoothly
  
  #Merge the rea_backgroundFieldSaltData and rea_fieldData tables
  if(exists("rea_backgroundFieldSaltData")){
    loggerSiteData <- merge(rea_backgroundFieldSaltData, 
                            rea_fieldData, 
                            by = c('siteID', 'startDate'), 
                            all = T)
  }else{
    loggerSiteData <- merge(rea_backgroundFieldCondData, 
                            rea_fieldData, 
                            by = c('siteID', 'startDate'), 
                            all = T)
  }
  
  #Create input file for reaeration calculations
  outputDFNames <- c(
    'siteID',
    'namedLocation', #Station at this point
    'startDate',
    'stationToInjectionDistance',
    'injectionType',
    'slugTracerMass',
    'slugPourTime',
    'dripStartTime',
    'backgroundSaltConc',
    'plateauSaltConc',
    'corrPlatSaltConc',
    'plateauGasConc',
    'wettedWidth',
    'waterTemp',
    'hoboSampleID',
    'fieldDischarge',
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
  
  #Eliminated this so that the calc function can differentiate between NaBr and model
  #Change to more generic injection types
  #outputDF$injectionType[outputDF$injectionType == "NaCl"] <- "constant"
  #outputDF$injectionType[outputDF$injectionType == "NaBr"|outputDF$injectionType == "model"] <- "slug"
  
  outputDF$eventID <- paste0(outputDF$siteID, ".", substr(outputDF$startDate,1,4), substr(outputDF$startDate,6,7), substr(outputDF$startDate,9,10))
  
  #Remove data for model type injections
  outputDF <- outputDF[outputDF$injectionType!="model"&!is.na(outputDF$injectionType),]
  
  QFile <- def.format.Q(dataDir = dataDir, site = site)
  QFile <- def.calc.Q.inj(QFile)
  #Move the slug calculations for Q to peakTime so that users only have to click on the plots once
  #QFile <- def.calc.Q.slug(inputFile = QFile, dataDir = filepath)
  
  if(fieldQ){
    for(i in unique(outputDF$eventID)){
      #print(i)
      currQ <- dsc_fieldData$totalDischarge[dsc_fieldData$eventID == i]
      currUnits <- dsc_fieldData$totalDischargeUnits[dsc_fieldData$eventID == i]
      if(length(currUnits)>0 && currUnits == "cubicMetersPerSecond"){
        currQ <- currQ * 1000
      }
      try(outputDF$fieldDischarge[outputDF$eventID == i] <- currQ, silent = T)
    }
  }
  
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
      try(outputDF$constDischarge[i] <- QFile$Q.inj[QFile$namedLocation == station
                                           & QFile$startDate == startDate], silent = T)
    }
    #Moved this to the peakTime function step in the calc reaeration
    # else if(!is.na(outputDF$injectionType[i]) && 
    #          (outputDF$injectionType[i] == "slug")){
    #   try(outputDF$slugDischarge[i] <- QFile$Q.slug[QFile$namedLocation == station
    #                                         & QFile$startDate == startDate], silent = T)
    # }
    
    #Fill in hoboSampleID from background logger table
    try(outputDF$hoboSampleID[i] <- rea_backgroundFieldCondData$hoboSampleID[
      rea_backgroundFieldCondData$namedLocation == station &
        rea_backgroundFieldCondData$startDate == startDate], silent = T)
    
    #Fill in background concentration data
    try(outputDF$backgroundSaltConc[i] <- 
          unique(rea_externalLabDataSalt$finalConcentration[
            rea_externalLabDataSalt$namedLocation == station &
              rea_externalLabDataSalt$startDate == startDate & 
              grepl(paste0(".B",substr(station,nchar(station),nchar(station)),"."),
                    rea_externalLabDataSalt$saltSampleID)]), silent = T)
    
    #Fill in 0 for background salt concentration if it is below detection
    if(length(rea_externalLabDataSalt$saltBelowDetectionQF[
      rea_externalLabDataSalt$namedLocation == station &
      rea_externalLabDataSalt$startDate == startDate & 
      grepl(paste0(".B",substr(station,nchar(station),nchar(station)),"."),
            rea_externalLabDataSalt$saltSampleID)])>0&&(
       unique(rea_externalLabDataSalt$saltBelowDetectionQF[
      rea_externalLabDataSalt$namedLocation == station &
      rea_externalLabDataSalt$startDate == startDate & 
      grepl(paste0(".B",substr(station,nchar(station),nchar(station)),"."),
            rea_externalLabDataSalt$saltSampleID)])==1&
      is.na(unique(rea_externalLabDataSalt$finalConcentration[
        rea_externalLabDataSalt$namedLocation == station &
        rea_externalLabDataSalt$startDate == startDate & 
        grepl(paste0(".B",substr(station,nchar(station),nchar(station)),"."),
              rea_externalLabDataSalt$saltSampleID)])))){
      outputDF$backgroundSaltConc[i] <- 0
    }
    
    #Fill in plateau concentration data for constant rate injection
    pSaltConc <- rea_externalLabDataSalt$finalConcentration[
      rea_externalLabDataSalt$namedLocation == station &
        rea_externalLabDataSalt$startDate == startDate & 
        grepl(repRegex, rea_externalLabDataSalt$saltSampleID)]
    
    #Remove outliers TBD
    #Instead of calculating the mean, concatenate all values for plotting and assessment
    outputDF$plateauSaltConc[i] <- paste(pSaltConc, collapse = "|")
    #Calculate the mean plateau concentration
    #outputDF$plateauSaltConc[i] <- ifelse(!is.nan(mean(pSaltConc, na.rm = T)),mean(pSaltConc, na.rm = T),NA)
    
    #Fill in plateau gas concentration
    pGasConc <- rea_externalLabDataGas$gasTracerConcentration[
      rea_externalLabDataGas$namedLocation == station &
        rea_externalLabDataGas$startDate == startDate & 
        grepl(repRegex, rea_externalLabDataGas$gasSampleID)]
    
    #Remove outliers TBD
    #Instead of calculating the mean, concatenate all values for plotting and assessment
    outputDF$plateauGasConc[i] <- paste(pGasConc, collapse = "|")
    #Calculate the mean plateau concentration
    #outputDF$plateauGasConc[i] <- ifelse(!is.nan(mean(pGasConc, na.rm = T)),mean(pGasConc, na.rm = T),NA)
    
    #Fill in mean wetted width
    wettedWidthVals <- rea_widthFieldData$wettedWidth[
      rea_widthFieldData$namedLocation == siteID &
        grepl(substr(startDate, 1, 10), rea_widthFieldData$startDate)]
    
    #Remove outliers TBD
    #Calculate the mean wetted width
    outputDF$wettedWidth[i] <- ifelse(!is.nan(mean(wettedWidthVals, na.rm = T)),mean(wettedWidthVals, na.rm = T),NA)
    
    #Fill in mean water temp
    tempVals <- rea_plateauMeasurementFieldData$waterTemp[
      rea_plateauMeasurementFieldData$siteID == siteID &
        grepl(substr(startDate, 1, 10), rea_plateauMeasurementFieldData$startDate)]
    
    #Remove outliers TBD
    #Calculate the mean water temp
    outputDF$waterTemp[i] <- ifelse(!is.nan(mean(tempVals, na.rm = T)),mean(tempVals, na.rm = T),NA)
    
  }
  
  #Assume that if the background concentration is below detection that the value is 0 for corrections
  # outputDF$corrPlatSaltConc[!is.na(outputDF$plateauSaltConc) & !is.na(outputDF$backgroundSaltConc)] <- 
  #   as.numeric(outputDF$plateauSaltConc[!is.na(outputDF$plateauSaltConc) & !is.na(outputDF$backgroundSaltConc)]) - 
  #   as.numeric(outputDF$backgroundSaltConc[!is.na(outputDF$plateauSaltConc) & !is.na(outputDF$backgroundSaltConc)])
  
  return(outputDF)

  
}

