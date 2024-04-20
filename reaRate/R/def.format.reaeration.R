##############################################################################################
#' @title Formats reaeration data for rate calculations

#' @author
#' Kaelin M. Cawley \email{kcawley@battelleecology.org} \cr

#' @description This function formats data from the NEON reaeration data product to calculate
#' loss rate, travel time, SF6 reaeration rate, O2 gas transfer velocity, and Schmidt number 600.
#' Either the basic or expanded package can be downloaded. The data files need to be loaded
#' into the R environment

#' @importFrom stageQCurve conv.calc.Q
#' @importFrom geoNEON getLocBySite
#' @importFrom utils read.csv

#' @param rea_backgroundFieldCondData This dataframe contains the data for the NEON rea_backgroundFieldCondData table [dataframe]
#' @param rea_backgroundFieldSaltData This dataframe contains the data for the NEON rea_backgroundFieldSaltData table [dataframe]
#' @param rea_fieldData This dataframe contains the data for the NEON rea_fieldData table [dataframe]
#' @param rea_plateauMeasurementFieldData This dataframe contains the data for the NEON rea_plateauMeasurementFieldData table [dataframe]
#' @param rea_plateauSampleFieldData This dataframe contains the data for the NEON rea_plateauSampleFieldData table [dataframe]
#' @param rea_externalLabDataSalt This dataframe contains the data for the NEON rea_externalLabDataSalt table [dataframe]
#' @param rea_externalLabDataGas This dataframe contains the data for the NEON rea_externalLabDataGas table [dataframe]
#' @param rea_widthFieldData This dataframe contains the data for the NEON rea_widthFieldData table [dataframe]
#' @param dsc_fieldData This dataframe contains the data for the NEON dsc_fieldData table, optional if there is a dsc_fieldDataADCP table [dataframe]
#' @param dsc_individualFieldData This dataframe contains the data for the NEON dsc_individualFieldData table, optional [dataframe]
#' @param dsc_fieldDataADCP This dataframe contains the data for the NEON dsc_fieldDataADCP table, optional[dataframe]
#' @param waq_instantaneous This dataframe contains the data for the NEON Water Quality sensor data, optional[dataframe]

#' @return This function returns one data frame formatted for use with def.calc.reaeration.R

#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

#' @keywords surface water, streams, rivers, reaeration, gas transfer velocity, schmidt number

#' @examples
#' #TBD

#' @seealso def.calc.tracerTime.R for calculating the stream travel time,
#' def.plot.reaQcurve.R for plotting reaeration rate versus stream flow

#' @export

# changelog and author contributions / copyrights
#   Kaelin M. Cawley (2017-10-30)
#     original creation
#   Kaelin M. Cawley (2018-05-03)
#     added the option of getting data from the API rather than a file download
#   Kaelin M. Cawley (2020-12-03)
#     updated to allow users to use data already loaded to R since there are so many options
#     of how to get it there now
#   Kaelin M. Cawley
#     Updated to fix a few user bugs (may have been mac specific) and include model injection types
#   Kaelin M. Cawley (2021-06-21)
#     Update to fix a bug with merging model injections
#   Kaelin M. Cawley (2022-02-11)
#     Additional data cleaning and bug fix for https://github.com/NEONScience/NEON-reaeration/issues/9
##############################################################################################
def.format.reaeration <- function(
  rea_backgroundFieldCondData,
  rea_backgroundFieldSaltData = NULL,
  rea_fieldData,
  rea_plateauMeasurementFieldData,
  rea_plateauSampleFieldData,
  rea_externalLabDataSalt,
  rea_externalLabDataGas,
  rea_widthFieldData,
  dsc_fieldData = NULL,
  dsc_individualFieldData = NULL,
  dsc_fieldDataADCP = NULL,
  waq_instantaneous = NULL
) {
  
  if(is.null(dsc_fieldData) & is.null(dsc_fieldDataADCP)){
    resp <- readline("Discharge data not loaded or available. Reaeration rates cannot be determined. Do you want to continue to calculate travel time and SF6 loss rate only? y/n: ")
    if(resp %in% c("n","N")) {
      stop("Input data will not be used to make any calculations. Exiting.")
    }
    # if(!(resp %in% c("y","Y"))) {
    #   stop("Input data will not be used to make any calculations. Exiting.")
    # }
  }
  
  # #Hopefully I can comment this out in the future when we get a siteID column, but for now, I'll make one
  # dsc_fieldDataADCP$siteID <- dsc_fieldDataADCP$stationID
  
  # Remove all samplingImpractical records from rea_fieldData
  rea_fieldData <- rea_fieldData[is.na(rea_fieldData$samplingImpractical),]
  
  # Pull the timezone for the site(s) for making sure the eventIDs match depending on the time of day, need to convert to local time.
  allSites <- unique(rea_fieldData$siteID)
  
  rea_fieldData$localDate <- NA
  dsc_fieldData$localDate <- NA
  dsc_fieldDataADCP$localDate <- NA
  rea_plateauMeasurementFieldData$localDate <- NA
  rea_backgroundFieldCondData$localDate <- NA
  for(currSite in allSites){
    currLocInfo <- geoNEON::getLocBySite(site = currSite)
    currTimeZone <- as.character(currLocInfo$siteTimezone)
    
    rea_fieldData$localDate[rea_fieldData$siteID == currSite] <- format(rea_fieldData$collectDate, tz = currTimeZone, format = "%Y%m%d")
    dsc_fieldData$localDate[dsc_fieldData$siteID == currSite] <- format(dsc_fieldData$collectDate, tz = currTimeZone, format = "%Y%m%d")
    dsc_fieldDataADCP$localDate[dsc_fieldDataADCP$siteID == currSite] <- format(dsc_fieldDataADCP$endDate, tz = currTimeZone, format = "%Y%m%d")
    rea_plateauMeasurementFieldData$localDate[rea_plateauMeasurementFieldData$siteID == currSite] <- format(rea_plateauMeasurementFieldData$collectDate, tz = currTimeZone, format = "%Y%m%d")
    rea_backgroundFieldCondData$localDate[rea_backgroundFieldCondData$siteID == currSite] <- format(rea_backgroundFieldCondData$startDate, tz = currTimeZone, format = "%Y%m%d")
  }
  
  # Add an eventID for later
  rea_fieldData$eventID <- paste(rea_fieldData$siteID, rea_fieldData$localDate, sep = ".")
  dsc_fieldData$eventID <- paste(dsc_fieldData$siteID, dsc_fieldData$localDate, sep = ".")
  rea_plateauMeasurementFieldData$eventID <- paste(rea_plateauMeasurementFieldData$siteID, rea_plateauMeasurementFieldData$localDate, sep = ".")
  rea_backgroundFieldCondData$eventID <- paste(rea_backgroundFieldCondData$siteID, rea_backgroundFieldCondData$localDate, sep = ".")
  dsc_fieldDataADCP$eventID <- paste(dsc_fieldDataADCP$siteID, dsc_fieldDataADCP$localDate, sep = ".")
  
  rea_fieldData$namedLocation <- NULL #So that merge goes smoothly
  rea_backgroundFieldCondData$collectDate <- rea_backgroundFieldCondData$startDate #Also to smooth merging
  
  # Populate the saltBelowDetectionQF if it isn't there and remove any values with flags of 1
  rea_externalLabDataSalt$saltBelowDetectionQF[is.na(rea_externalLabDataSalt$saltBelowDetectionQF)] <- 0
  rea_externalLabDataSalt$finalConcentration[rea_externalLabDataSalt$saltBelowDetectionQF == 1] <- NA
  
  #Format the date for the sensor data so that it matches the REA data
  waq_instantaneous$startDateTimeTrim <- format(waq_instantaneous$startDateTime, format = "%Y-%m-%d %H:%M")
  
  #Merge the rea_backgroundFieldSaltData, rea_backgroundFieldCondData, and rea_fieldData tables to handle the model injections
  if(!is.null(rea_backgroundFieldSaltData)){
    loggerSiteData <- merge(rea_backgroundFieldSaltData,
                            rea_fieldData,
                            by = c('siteID', 'collectDate'),
                            all = TRUE)
  }else{
    loggerSiteData <- merge(rea_backgroundFieldCondData,
                            rea_fieldData,
                            by = c('siteID', 'collectDate', 'eventID'),
                            all = TRUE)
  }
  
  #Add in station if it's missing for a model injectionType
  missingStations <- loggerSiteData$eventID[which(is.na(loggerSiteData$namedLocation))]
  if(length(missingStations) > 0){
    loggerSiteData <- merge(loggerSiteData,
                            rea_backgroundFieldCondData[rea_backgroundFieldCondData$eventID %in% missingStations,],
                            by = c("siteID","collectDate"),
                            all = TRUE)
    loggerSiteData$namedLocation <- loggerSiteData$namedLocation.x
    loggerSiteData$namedLocation[is.na(loggerSiteData$namedLocation.x)] <- loggerSiteData$namedLocation.y[is.na(loggerSiteData$namedLocation.x)]
    
    #Add back in a few variables that got messed up with the bonus merge step
    loggerSiteData$stationToInjectionDistance <- loggerSiteData$stationToInjectionDistance.x
    loggerSiteData$eventID <- loggerSiteData$eventID.x
  }
  
  # Remove loggerSiteData rows with no injectionType
  loggerSiteData <- loggerSiteData[!is.na(loggerSiteData$injectionType),]
  
  
  #Add the injectate, background, or plateau type to external lab data
  #This is probably going to have to change with the switch to using only barcodes!
  rea_externalLabDataSalt$sampleType <- NA
  rea_externalLabDataSalt$sampleType[rea_externalLabDataSalt$saltSampleID %in% rea_fieldData$injectateSampleID] <- "injectate"
  rea_externalLabDataSalt$sampleType[rea_externalLabDataSalt$saltSampleID %in% rea_backgroundFieldSaltData$saltBackgroundSampleID] <- "background"
  rea_externalLabDataSalt$sampleType[rea_externalLabDataSalt$saltSampleID %in% rea_plateauSampleFieldData$saltTracerSampleID] <- "plateau"
  
  #Create input file for reaeration calculations
  outputDFNames <- c(
    'siteID',
    'namedLocation', #Station at this point
    'collectDate',
    'collectDateTrim',
    'stationToInjectionDistance',
    'injectionType',
    'slugTracerMass',
    'slugPourTime',
    'dripStartTime',
    'plateauCollectTime',
    'backgroundSaltConc', #One sample per station per experiment
    'meanBackgroundCond', # Average of field measurements
    'backgroundSensorCond', # Conductivity pulled from sensor data prior to drip/slug
    'plateauSaltConc', # Concatenated string of plateau salt from lab
    'meanPlatSaltConc', # Mean of plateauSaltConc entries
    'sdPlatSaltConc', # Standard deviation of plateauSaltConc entries
    'plateauSaltConcClean', # Outliers removed from plateauSaltConc entries, concatenated string
    'meanPlatSaltConcClean', # Mean of plateauSaltConcClean
    'sdPlatSaltConcClean', # Standard deviation of plateauSaltConcClean
    'plateauSaltConcCleanCorr', # plateauSaltConcClean - backgroundSaltConc, concatenated string
    'meanPlatSaltConcCleanCorr', # Mean plateauSaltConcCleanCorr
    'sdPlatSaltConcCleanCorr', # Standard deviation plateauSaltConcCleanCorr
    'platSensorCond', # Conductivity pulled from sensor data at plateau time
    'plateauGasConc', # Concatenated string of plateau gas from lab
    'meanPlatGasConc', # Mean of plateauGasConc
    'sdPlatGasConc', # Standard deviation of plateauGasConc
    'plateauGasConcClean', #Outliers removed from plateauGasConc entries, concatenated string
    'meanPlatGasConcClean', # Mean of plateauGasConcClean
    'sdPlatGasConcClean', # Standard deviation of plateauGasConcClean
    'corrPlatGasConcClean', # plateauGasConcClean/(plateauSaltConcClean - backgroundSaltConc), matched gas and salt syringes
    'meanCorrPlatGasConcClean', # plateauGasConcClean/meanPlatSaltConcCleanCorr, all gas divided by the same mean salt from station
    'unmixedStationFlag',
    'wettedWidth',
    'waterTemp',
    'hoboSampleID',
    'fieldDischarge_lps',
    'stage',
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
  
  # #Remove data for model type injections since we can't get k values from those anyway
  # modelInjectionTypes <- c("model","model - slug","model - CRI")
  # outputDF <- outputDF[!outputDF$injectionType%in%modelInjectionTypes & !is.na(outputDF$injectionType),]
  
  # This is no longer needed since the finalDischarge field now contains the recalculated discharge values
  # #Recalculate wading survey discharge using the stageQCurve package and then add to the output dataframe
  # dsc_fieldData_calc <- stageQCurve::conv.calc.Q(stageData = dsc_fieldData,
  #                                                dischargeData = dsc_individualFieldData)
  
  #Populate Q and stage from wading surveys
  for(i in unique(outputDF$eventID)){
    #print(i)
    currQ <- dsc_fieldData$finalDischarge[dsc_fieldData$eventID == i]
    try(outputDF$fieldDischarge_lps[outputDF$eventID == i] <- currQ, silent = T)
    
    currStage <- dsc_fieldData$streamStage[dsc_fieldData$eventID == i]
    try(outputDF$stage[outputDF$eventID == i] <- currStage, silent = T)
  }
  
  #Populate Q and stage from ADCP data, if applicable
  if(exists('dsc_fieldDataADCP')){
    for(i in unique(outputDF$eventID)){
      #print(i)
      currQ <- dsc_fieldDataADCP$totalDischarge[dsc_fieldDataADCP$eventID == i]
      currQUnits <- dsc_fieldDataADCP$totalDischargeUnits[dsc_fieldDataADCP$eventID == i]
      if(length(currQUnits) > 0 && currQUnits == "cubicMetersPerSecond"){
        currQ <- currQ * 1000 # Convert to lps
      }
      try(outputDF$fieldDischarge_lps[outputDF$eventID == i] <- currQ, silent = T)
      
      currStage <- dsc_fieldDataADCP$streamStage[dsc_fieldDataADCP$eventID == i]
      try(outputDF$stage[outputDF$eventID == i] <- currStage, silent = T)
    }
  }
  
  
  #Format the collect date for matching with sensor data
  outputDF$collectDateTrim <- format(outputDF$collectDate, format = "%Y-%m-%d %H:%M")
  
  #Loop through all the records to populate the other fields
  for(i in seq(along = outputDF$siteID)){
    siteID <- outputDF$siteID[i]
    startDate <- outputDF$collectDate[i]
    station <- outputDF$namedLocation[i]
    stationType <- substr(station, 6, nchar(station))
    
    #Fill in hoboSampleID from background logger table
    try(outputDF$hoboSampleID[i] <- rea_backgroundFieldCondData$hoboSampleID[
      rea_backgroundFieldCondData$namedLocation == station &
        rea_backgroundFieldCondData$startDate == startDate], silent = T)
    
    #Fill in background concentration data
    try(outputDF$backgroundSaltConc[i] <- rea_externalLabDataSalt$finalConcentration[
      rea_externalLabDataSalt$namedLocation == station &
        rea_externalLabDataSalt$startDate == startDate &
        rea_externalLabDataSalt$sampleType == "background"], silent = T)
    
    #Fill in background conductivity data
    try(outputDF$meanBackgroundCond[i] <- mean(rea_backgroundFieldSaltData$specificConductanceRep1[rea_backgroundFieldSaltData$namedLocation == station &
                                                                                                    rea_backgroundFieldSaltData$collectDate == startDate],
                                              rea_backgroundFieldSaltData$specificConductanceRep2[rea_backgroundFieldSaltData$namedLocation == station &
                                                                                                    rea_backgroundFieldSaltData$collectDate == startDate],
                                              rea_backgroundFieldSaltData$specificConductanceRep3[rea_backgroundFieldSaltData$namedLocation == station &
                                                                                                    rea_backgroundFieldSaltData$collectDate == startDate],
                                              rea_backgroundFieldSaltData$specificConductanceRep4[rea_backgroundFieldSaltData$namedLocation == station &
                                                                                                    rea_backgroundFieldSaltData$collectDate == startDate],na.rm = TRUE), silent = T)
    
    #Populate sensor data for the background time (20 minutes prior to the injection start time)
    backgroundDate <- ifelse(outputDF$injectionType[i] %in% c("NaCl","model - CRI"),format(outputDF$dripStartTime[i] - 20*60, format = "%Y-%m-%d %H:%M"),format(outputDF$slugPourTime[i] - 20*60, format = "%Y-%m-%d %H:%M"))
    
    #Determine the best horizontal and vertical indices to use
    
    HOROptions <- waq_instantaneous$horizontalPosition[waq_instantaneous$startDateTimeTrim == backgroundDate]
    
    if("101" %in% HOROptions){
      S1Hor <- "101"
    }else{
      S1Hor <- "111"
    }
    
    if("102" %in% HOROptions){
      S2Hor <- "102"
    }else{
      S2Hor <- "112"
    }
    
    # Pull out the background data
    if(station == paste0(currSite,".AOS.reaeration.station.01")){
      try(outputDF$backgroundSensorCond[i] <- waq_instantaneous$specificConductance[waq_instantaneous$horizontalPosition == S1Hor & waq_instantaneous$startDateTimeTrim == backgroundDate], silent = TRUE)
    }else if(station == paste0(currSite,".AOS.reaeration.station.04")){
      try(outputDF$backgroundSensorCond[i] <- waq_instantaneous$specificConductance[waq_instantaneous$horizontalPosition == S2Hor & waq_instantaneous$startDateTimeTrim == backgroundDate], silent = TRUE)
    }
    
    #Fill in plateau sample collection time
    try(outputDF$plateauCollectTime[i] <- unique(rea_plateauSampleFieldData$plateauCollectTime[
      rea_plateauSampleFieldData$namedLocation == station &
        rea_plateauSampleFieldData$startDate == startDate]), silent = T)
    
    #Populate sensor data for the sample collection time
    platCollectDateTrim <- format(as.POSIXct(outputDF$plateauCollectTime[i], origin = "1970-01-01", tz = "GMT"), format = "%Y-%m-%d %H:%M")
    if(station == paste0(currSite,".AOS.reaeration.station.01")){
      try(outputDF$platSensorCond[i] <- waq_instantaneous$specificConductance[waq_instantaneous$horizontalPosition == S1Hor & waq_instantaneous$startDateTimeTrim == platCollectDateTrim], silent = TRUE)
    }else if(station == paste0(currSite,".AOS.reaeration.station.04")){
      try(outputDF$platSensorCond[i] <- waq_instantaneous$specificConductance[waq_instantaneous$horizontalPosition == S2Hor & waq_instantaneous$startDateTimeTrim == platCollectDateTrim], silent = TRUE)
    }
    
    #Fill in plateau concentration data for constant rate injection
    # Need to join with the field data rather than use the sampleID since we're switching to barcodes only
    pSaltConc <- rea_externalLabDataSalt$finalConcentration[
      rea_externalLabDataSalt$namedLocation == station &
        rea_externalLabDataSalt$startDate == startDate &
        rea_externalLabDataSalt$sampleType == "plateau"]
    
    #Calculate a mean concentration for plateau salt
    outputDF$meanPlatSaltConc[i] <- mean(pSaltConc, na.rm = TRUE)
    outputDF$sdPlatSaltConc[i] <- stats::sd(pSaltConc, na.rm = TRUE)
    
    #Concatenate all values for plotting and assessment
    outputDF$plateauSaltConc[i] <- paste(pSaltConc, collapse = "|")
    
    #Fill in plateau gas concentration
    pGasConc <- rea_externalLabDataGas$gasTracerConcentration[
      rea_externalLabDataGas$namedLocation == station &
        rea_externalLabDataGas$startDate == startDate]
    
    #Calculate a mean concentration for plateau gas
    outputDF$meanPlatGasConc[i] <- mean(pGasConc, na.rm = TRUE)
    outputDF$sdPlatGasConc[i] <- stats::sd(pGasConc, na.rm = TRUE)
    
    #Concatenate all values for plotting and assessment
    outputDF$plateauGasConc[i] <- paste(pGasConc, collapse = "|")
    
    #Remove outliers from plateau salt and be sure to remove gas as well to maintain matches
    validSaltOut <- TRUE
    
    if(all(pSaltConc == 0)|all(is.na(pSaltConc))){
      cat("\t",outputDF$eventID[i],"contains all 0 or NA values for salt data and cannot run outlier detection.\n")
      validSaltOut <- FALSE
    }
    if(length(pSaltConc) < 3 & validSaltOut){
      cat("\t",outputDF$eventID[i],"contains fewer than 3 replicates for salt data and cannot run outlier detection.\n")
      validSaltOut <- FALSE
    }
    #if(validSaltOut){
      platSaltTest <- graphics::boxplot(pSaltConc, plot = FALSE)
    #}
    
    #Remove outliers from plateau gas
    validGasOut <- TRUE
    
    if(all(pGasConc == 0)|all(is.na(pGasConc))){
      cat("\t",outputDF$eventID[i],"contains all 0 or NA values for gas data and cannot run outlier detection.\n")
      validGasOut <- FALSE
    }
    if(length(pGasConc) < 3 & validGasOut){
      cat("\t",outputDF$eventID[i],"contains fewer than 3 replicates for gas data and cannot run outlier detection.\n")
      validGasOut <- FALSE
    }
    #if(validGasOut){
      platGasTest <- graphics::boxplot(pGasConc, plot = FALSE)
    #}
    
    # Remove both salt and gas if either is an outlier for a pair
    if(length(platSaltTest$out) == 0 && length(platGasTest$out) == 0){
      outputDF$plateauSaltConcClean[i] <- outputDF$plateauSaltConc[i]
      outputDF$meanPlatSaltConcClean[i] <- outputDF$meanPlatSaltConc[i]
      outputDF$sdPlatSaltConcClean[i] <- outputDF$sdPlatSaltConc[i]
      
      outputDF$plateauGasConcClean[i] <- outputDF$plateauGasConc[i]
      outputDF$meanPlatGasConcClean[i] <- outputDF$meanPlatGasConc[i]
      outputDF$sdPlatGasConcClean[i] <- outputDF$sdPlatGasConc[i]
      
      saltForBackCor <- pSaltConc
    }else{
      saltOutlierIdxs <- which(!pSaltConc %in% platSaltTest$out)
      gasOutlierIdxs <- which(!pGasConc %in% platGasTest$out)
      idxToKeep <- intersect(saltOutlierIdxs, gasOutlierIdxs)
      
      outputDF$plateauSaltConcClean[i] <- paste(pSaltConc[idxToKeep], collapse = "|")
      outputDF$meanPlatSaltConcClean[i] <- mean(pSaltConc[idxToKeep], na.rm = TRUE)
      outputDF$sdPlatSaltConcClean[i] <- stats::sd(pSaltConc[idxToKeep], na.rm = TRUE)
      
      outputDF$plateauGasConcClean[i] <- paste(pGasConc[idxToKeep], collapse = "|")
      outputDF$meanPlatGasConcClean[i] <- mean(pGasConc[idxToKeep], na.rm = TRUE)
      outputDF$sdPlatGasConcClean[i] <- stats::sd(pGasConc[idxToKeep], na.rm = TRUE)
      
      saltForBackCor <- pSaltConc[idxToKeep]
    }
    
    #Background subtract salt data
    try(outputDF$plateauSaltConcCleanCorr[i] <- paste((saltForBackCor - outputDF$backgroundSaltConc[i]), collapse = "|"))
    try(outputDF$meanPlatSaltConcCleanCorr[i] <- mean((saltForBackCor - outputDF$backgroundSaltConc[i]), na.rm = TRUE))
    try(outputDF$sdPlatSaltConcCleanCorr[i] <- stats::sd((saltForBackCor - outputDF$backgroundSaltConc[i]), na.rm = TRUE))
    
    #Salt correct the cleaned gas concentrations
    try(pGasConForCorr <-  as.numeric(strsplit(outputDF$plateauGasConcClean[i],"\\|")[[1]]))
    
    #Matched syringes
    try(pSaltConForCorr <- as.numeric(strsplit(outputDF$plateauSaltConcCleanCorr[i],"\\|")[[1]]))
    try(outputDF$corrPlatGasConcClean[i] <- paste((pGasConForCorr/pSaltConForCorr), collapse = "|"))
      
    #Using the mean salt at the station
    try(outputDF$meanCorrPlatGasConcClean[i] <- paste((pGasConForCorr/outputDF$meanPlatSaltConcCleanCorr[i]), collapse = "|"))
      
    #Flag salt data for unmixed situations
    platSaltCV <- outputDF$sdPlatSaltConcClean[i]/outputDF$meanPlatSaltConcClean[i]
    if(!is.na(platSaltCV) && platSaltCV > 0.1){
      outputDF$unmixedStationFlag[i] <- "unmixedSaltStation"
    }
    
    #Flag gas data for unmixed situations
    platGasCV <- outputDF$sdPlatGasConcClean[i]/outputDF$meanPlatGasConcClean[i]
    if(!is.na(platGasCV) && platGasCV > 0.1){
      if(is.na(outputDF$unmixedStationFlag[i])){
        outputDF$unmixedStationFlag[i] <- "unmixedGasStation"
      }else{
        outputDF$unmixedStationFlag[i] <- "unmixedSaltStation|unmixedGasStation"
      }
    }
    
    #Fill in mean wetted width
    wettedWidthVals <- rea_widthFieldData$wettedWidth[
      rea_widthFieldData$namedLocation == siteID &
        grepl(substr(startDate, 1, 10), rea_widthFieldData$collectDate)]
    
    #Remove outliers TBD
    #Calculate the mean wetted width
    outputDF$wettedWidth[i] <- ifelse(!is.nan(mean(wettedWidthVals, na.rm = T)),mean(wettedWidthVals, na.rm = T),NA)
    
    #Populate water temp
    suppressWarnings(try(outputDF$waterTemp[i] <- rea_plateauMeasurementFieldData$waterTemp[rea_plateauMeasurementFieldData$namedLocation == station &
                                                                                              rea_plateauMeasurementFieldData$eventID == outputDF$eventID[i]], silent = TRUE))
    
  }
  
  #Remove any rows where injectionType is missing
  outputDF <- outputDF[!is.na(outputDF$injectionType),]
  
  return(outputDF)
  
}

