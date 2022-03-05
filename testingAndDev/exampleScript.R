##############################################################################################
#' @title Example Script for use with neonUtilities 2.0+

#' @author
#' Kaelin M. Cawley \email{kcawley@battelleecology.org} \cr

#' @description This script downloads and formats reaeration and discharge data from the 
#' NEON data portal in order to calculate loss rate, travel time, SF6 reaeration rate, 
#' O2 gas transfer velocity, and Schmidt number 600.

#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

#' @keywords surface water, streams, rivers, reaeration, gas transfer velocity, schmidt number

# changelog and author contributions / copyrights
#   Kaelin M. Cawley (2021-02-04)
#     original creation
#	Kaelin M. Cawley (2022-03-05)
#	  updated to work with reaRate v1.0.0+
##############################################################################################

#User Inputs
siteID <- "GUIL"
plotPath <- paste0("~/reaOutputs/",siteID,"/QAQC_plots")

#String constants
reaDPID <- "DP1.20190.001"
dscDPID <- "DP1.20048.001"
wqDPID <- "DP1.20288.001"

# Download Reaeration Data (just delete the input for dates to get data for all time for a site)
reaInputList <- neonUtilities::loadByProduct(dpID = reaDPID, 
                                             site = siteID,
                                             startdate = "2019-01-01", 
                                             enddate = "2021-12-01",
                                             check.size = FALSE)

rea_backgroundFieldCondDataIn <- reaInputList$rea_backgroundFieldCondData
rea_backgroundFieldSaltDataIn <- reaInputList$rea_backgroundFieldSaltData
rea_fieldDataIn <- reaInputList$rea_fieldData
rea_plateauMeasurementFieldDataIn <- reaInputList$rea_plateauMeasurementFieldData
rea_plateauSampleFieldDataIn <- reaInputList$rea_plateauSampleFieldData
rea_externalLabDataSaltIn <- reaInputList$rea_externalLabDataSalt
rea_externalLabDataGasIn <- reaInputList$rea_externalLabDataGas
rea_widthFieldDataIn <- reaInputList$rea_widthFieldData

# Download Discharge Data
qInputList <- neonUtilities::loadByProduct(dpID = dscDPID, site = siteID, check.size = FALSE)

dsc_fieldDataIn <- qInputList$dsc_fieldData
dsc_individualFieldDataIn <- qInputList$dsc_individualFieldData
dsc_fieldDataADCPIn <- qInputList$dsc_fieldDataADCP

# Download Sensor Data
sensorData <- neonUtilities::loadByProduct(dpID = wqDPID, 
                                           site = siteID,
                                           check.size = FALSE)
waq_instantaneousIn <- sensorData$waq_instantaneous

# It can take a while to download data, so save it in case you need to go back
save.image(paste0("~/reaOutputs/",siteID,"/downloadedData.RData"))

# Format the downloaded data so everything is in one table
reaFormatted <- reaRate::def.format.reaeration(rea_backgroundFieldCondData = rea_backgroundFieldCondDataIn,
                                               rea_backgroundFieldSaltData = rea_backgroundFieldSaltDataIn,
                                               rea_fieldData = rea_fieldDataIn,
                                               rea_plateauMeasurementFieldData = rea_plateauMeasurementFieldDataIn,
                                               rea_plateauSampleFieldData = rea_plateauSampleFieldDataIn,
                                               rea_externalLabDataSalt = rea_externalLabDataSaltIn,
                                               rea_externalLabDataGas = rea_externalLabDataGasIn,
                                               rea_widthFieldData = rea_widthFieldDataIn,
                                               dsc_fieldData = dsc_fieldDataIn,
                                               dsc_individualFieldData = dsc_individualFieldDataIn,
                                               dsc_fieldDataADCP = dsc_fieldDataADCPIn,
                                               waq_instantaneous = waq_instantaneousIn)

# Example of a manual step... Remove some bad looking sensor data
reaFormatted$backgroundSensorCond[reaFormatted$backgroundSensorCond < 105] <- NA
reaFormatted$platSensorCond[reaFormatted$platSensorCond < 105] <- NA

# Calculate SF6 loss rates
plotsOut <- reaRate::gas.loss.rate.plot(inputFile = reaFormatted,
                                        savePlotPath = plotPath)

# Take a look at the background data
reaRate::bkgd.salt.conc.plot(inputFile = plotsOut,
                             savePlotPath = plotPath)

# Calculate travel times
reaRatesTrvlTime <- reaRate::def.calc.trvl.time(inputFile = plotsOut,
                                                loggerData = reaInputList$rea_conductivityFieldData,
                                                plot = TRUE,
                                                savePlotPath = plotPath)

# Example of a manual step... A few I didn't like where I picked the range, so trying again
badEventIDs <- c("GUIL.20200812", "GUIL.20201216")
reaFormattedTake2 <- reaRatesTrvlTime$inputFile[reaRatesTrvlTime$inputFile$eventID %in% badEventIDs,]
reaRatesTrvlTimeTake2 <- reaRate::def.calc.trvl.time(inputFile = reaFormattedTake2,
                                                     loggerData = reaInputList$rea_conductivityFieldData,
                                                     plot = TRUE,
                                                     savePlotPath = plotPath)
reaRatesTrvlTimeAll <- rbind(reaRatesTrvlTime$outputDF[!reaRatesTrvlTime$outputDF$eventID %in% badEventIDs,],
                             reaRatesTrvlTimeTake2$outputDF)

# Use SF6 loss rates for clean data to get k600 and K600
reaRatesCalc <- reaRate::def.calc.reaeration(inputFile = reaRatesTrvlTimeAll,
                                             lossRateSF6 = "slopeClean",
                                             outputSuffix = "clean")

# Use SF6 loss rates for background and plateaur salt corrected data to get k600 and K600
reaRatesCalc <- reaRate::def.calc.reaeration(inputFile = reaRatesCalc,
                                             lossRateSF6 = "slopeBackCorr",
                                             outputSuffix = "allCorr")

# Take a quick look at the different k660 v discharge depending on salt correction choice
plot(reaRatesCalc$meanQ,
     reaRatesCalc$k600.clean,
     pch = 19,
     xlab = "discharge",
     ylab = "k600")
points(reaRatesCalc$meanQ,
       reaRatesCalc$k600.allCorr,
       pch = 19,
       col = "blue")
legend("topright",
       legend = c("gas uncorrected","gas background salt corrected"),
       pch = 19,
       col = c("black", "blue"))


