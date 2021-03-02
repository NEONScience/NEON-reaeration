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
#' @param rea_externalLabDataSalt This dataframe contains the data for the NEON rea_externalLabDataSalt table [dataframe]
#' @param rea_externalLabDataGas This dataframe contains the data for the NEON rea_externalLabDataGas table [dataframe]
#' @param rea_widthFieldData This dataframe contains the data for the NEON rea_widthFieldData table [dataframe]
#' @param dsc_fieldData This dataframe contains the data for the NEON dsc_fieldData table [dataframe]
#' @param dsc_individualFieldData This dataframe contains the data for the NEON dsc_individualFieldData table, optional [dataframe]

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
##############################################################################################
def.format.reaeration <- function(
  rea_backgroundFieldCondData,
  rea_backgroundFieldSaltData,
  rea_fieldData,
  rea_plateauMeasurementFieldData,
  rea_externalLabDataSalt,
  rea_externalLabDataGas,
  rea_widthFieldData,
  dsc_fieldData,
  dsc_individualFieldData
) {

  if(!exists("dsc_fieldData")){
    resp <- readline("Discharge data not loaded or available. Reaeration rates cannot be determined. Do you want to continue to calculate travel time and SF6 loss rate only? y/n: ")
    if(resp %in% c("n","N")) {
      stop("Input data will not be used to make any calculations. Exiting.")
    }
    if(!(resp %in% c("y","Y"))) {
      stop("Input data will not be used to make any calculations. Exiting.")
    }
  }

  # Pull the timezone for the site(s) for making sure the eventIDs match depending on the time of day, need to convert to local time.
  allSites <- unique(rea_fieldData$siteID)

  rea_fieldData$localDate <- NA
  dsc_fieldData$localDate <- NA
  rea_plateauMeasurementFieldData$localDate <- NA
  for(currSite in allSites){
    currLocInfo <- geoNEON::getLocBySite(site = currSite)
    currTimeZone <- currLocInfo$siteTimezone

    rea_fieldData$localDate[rea_fieldData$siteID == currSite] <- format(rea_fieldData$collectDate, tz = currTimeZone, format = "%Y%m%d")

    dsc_fieldData$localDate[dsc_fieldData$siteID == currSite] <- format(dsc_fieldData$collectDate, tz = currTimeZone, format = "%Y%m%d")

    rea_plateauMeasurementFieldData$localDate[rea_plateauMeasurementFieldData$siteID == currSite] <- format(rea_plateauMeasurementFieldData$collectDate, tz = currTimeZone, format = "%Y%m%d")
  }

  # Add an eventID for later
  rea_fieldData$eventID <- paste(rea_fieldData$siteID, rea_fieldData$localDate, sep = ".")
  dsc_fieldData$eventID <- paste(dsc_fieldData$siteID, dsc_fieldData$localDate, sep = ".")
  rea_plateauMeasurementFieldData$eventID <- paste(rea_plateauMeasurementFieldData$siteID, rea_plateauMeasurementFieldData$localDate, sep = ".")

  rea_fieldData$namedLocation <- NULL #So that merge goes smoothly

  # Populate the saltBelowDetectionQF if it isn't there and remove any values with flags of 1
  rea_externalLabDataSalt$saltBelowDetectionQF[is.na(rea_externalLabDataSalt$saltBelowDetectionQF)] <- 0
  rea_externalLabDataSalt$finalConcentration[rea_externalLabDataSalt$saltBelowDetectionQF == 1] <- NA

  #Merge the rea_backgroundFieldSaltData and rea_fieldData tables
  if(exists("rea_backgroundFieldSaltData")){
    loggerSiteData <- merge(rea_backgroundFieldSaltData,
                            rea_fieldData,
                            by = c('siteID', 'collectDate'),
                            all = T)
  }else{
    loggerSiteData <- merge(rea_backgroundFieldCondData,
                            rea_fieldData,
                            by = c('siteID', 'collectDate'),
                            all = T)
  }

  #Create input file for reaeration calculations
  outputDFNames <- c(
    'siteID',
    'namedLocation', #Station at this point
    'collectDate',
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

  #Remove data for model type injections since we can't get k values from those anyway
  modelInjectionTypes <- c("model","model - slug","model - CRI")
  outputDF <- outputDF[!outputDF$injectionType%in%modelInjectionTypes & !is.na(outputDF$injectionType),]

  #Recalculate wading survey discharge using the stageQCurve package and then add to the output dataframe
  dsc_fieldData_calc <- stageQCurve::conv.calc.Q(stageData = dsc_fieldData,
                                                 dischargeData = dsc_individualFieldData)

  for(i in unique(outputDF$eventID)){
    #print(i)
    currQ <- dsc_fieldData_calc$calcQ[dsc_fieldData$eventID == i]
    try(outputDF$fieldDischarge[outputDF$eventID == i] <- currQ, silent = T)
  }

  for(i in seq(along = outputDF$siteID)){
    siteID <- outputDF$siteID[i]
    startDate <- outputDF$collectDate[i]
    station <- outputDF$namedLocation[i]
    stationType <- substr(station, 6, nchar(station))

    repRegex <- switch(stationType,
                       "AOS.reaeration.station.01" = "\\.0[12345]\\.",
                       "AOS.reaeration.station.02" = "\\.0[6789]\\.|\\.10\\.",
                       "AOS.reaeration.station.03" = "\\.1[12345]\\.",
                       "AOS.reaeration.station.04" = "\\.1[6789]\\.|\\.20\\."
    )

    #Fill in hoboSampleID from background logger table
    try(outputDF$hoboSampleID[i] <- rea_backgroundFieldCondData$hoboSampleID[
      rea_backgroundFieldCondData$namedLocation == station &
        rea_backgroundFieldCondData$startDate == startDate], silent = T)

    #Fill in background concentration data
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

    #Concatenate all values for plotting and assessment
    outputDF$plateauSaltConc[i] <- paste(pSaltConc, collapse = "|")

    #Fill in plateau gas concentration
    pGasConc <- rea_externalLabDataGas$gasTracerConcentration[
      rea_externalLabDataGas$namedLocation == station &
        rea_externalLabDataGas$startDate == startDate &
        grepl(repRegex, rea_externalLabDataGas$gasSampleID)]

    #Concatenate all values for plotting and assessment
    outputDF$plateauGasConc[i] <- paste(pGasConc, collapse = "|")

    #Fill in mean wetted width
    wettedWidthVals <- rea_widthFieldData$wettedWidth[
      rea_widthFieldData$namedLocation == siteID &
        grepl(substr(startDate, 1, 10), rea_widthFieldData$collectDate)]

    #Remove outliers TBD
    #Calculate the mean wetted width
    outputDF$wettedWidth[i] <- ifelse(!is.nan(mean(wettedWidthVals, na.rm = T)),mean(wettedWidthVals, na.rm = T),NA)

    #Populate water temp
    outputDF$waterTemp[i] <- rea_plateauMeasurementFieldData$waterTemp[rea_plateauMeasurementFieldData$namedLocation == station &
                                                                         rea_plateauMeasurementFieldData$eventID == outputDF$eventID[i]]

  }

  return(outputDF)

}

