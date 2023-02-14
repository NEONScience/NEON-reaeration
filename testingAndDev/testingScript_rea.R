##############################################################################################
#' @title Script to help with troubleshooting and updates

#' @author
#' Kaelin M. Cawley \email{kcawley@battelleecology.org} \cr

#' @description This script has all of the inputs typed out for ease of running scripts for 
#' troubleshooting and development.

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
siteID <- "COMO"
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
rea_backgroundFieldCondData = rea_backgroundFieldCondDataIn
rea_backgroundFieldSaltData = rea_backgroundFieldSaltDataIn
rea_fieldData = rea_fieldDataIn
rea_plateauMeasurementFieldData = rea_plateauMeasurementFieldDataIn
rea_plateauSampleFieldData = rea_plateauSampleFieldDataIn
rea_externalLabDataSalt = rea_externalLabDataSaltIn
rea_externalLabDataGas = rea_externalLabDataGasIn
rea_widthFieldData = rea_widthFieldDataIn
dsc_fieldData = dsc_fieldDataIn
dsc_individualFieldData = dsc_individualFieldDataIn
dsc_fieldDataADCP = dsc_fieldDataADCPIn
waq_instantaneous = waq_instantaneousIn

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

# Calculate SF6 loss rates
inputFile = reaFormatted
injectionTypeName = "injectionType"
eventID = "eventID"
stationToInjectionDistance = "stationToInjectionDistance"
plateauGasConc = "plateauGasConc"
corrPlatSaltConc = "corrPlatSaltConc"
savePlotPath = plotPathLossRates

plotsOut <- reaRate::gas.loss.rate.plot(inputFile = reaFormatted,
                                        savePlotPath = plotPath)

# Take a look at the background data
inputFile = plotsOut
savePlotPath = plotPath

reaRate::bkgd.salt.conc.plot(inputFile = plotsOut,
                             savePlotPath = plotPath)

# Calculate travel times
inputFile = plotsOut
loggerData = reaInputList$rea_conductivityFieldData
namedLocation = "namedLocation"
injectionTypeName = "injectionType"
eventID = "eventID"
slopeRaw = 'slopeRaw'
slopeClean = 'slopeClean'
slopeSaltCorr = 'slopeSaltCorr'
slopeNormClean = 'slopeNormClean'
slopeNormCorr = 'slopeNormCorr'
slopeBackCorr = 'slopeBackCorr'
stationToInjectionDistance = "stationToInjectionDistance"
hoboSampleID = "hoboSampleID"
slugPourTime = "slugPourTime"
dripStartTime = "dripStartTime"
meanBackgroundCond = "meanBackgroundCond"
discharge = "fieldDischarge"
waterTemp = "waterTemp"
wettedWidth = "wettedWidth"
plot = TRUE
savePlotPath = plotPath

reaRatesTrvlTime <- reaRate::def.calc.trvl.time(inputFile = plotsOut,
                                                loggerData = reaInputList$rea_conductivityFieldData,
                                                plot = TRUE,
                                                savePlotPath = plotPath)

# Use SF6 loss rates for clean data to get k600 and K600
inputFile = reaRatesTrvlTime$outputDF
lossRateSF6 = "lossRateSF6"
peakMaxVelocity = "peakMaxVelocity"
meanDepth = "meanDepth"
meanTemp = "meanTemp"
outputSuffix = "clean"

reaRatesCalc <- reaRate::def.calc.reaeration(inputFile = reaRatesTrvlTime$outputDF,
                                             lossRateSF6 = "slopeClean",
                                             outputSuffix = "clean")

# Take a quick look at the different k660 v discharge depending on salt correction choice
plot(reaRatesCalc$meanQ,
     reaRatesCalc$k600.clean,
     pch = 19,
     xlab = "discharge",
     ylab = "k600")
# points(reaRatesCalc$meanQ,
#        reaRatesCalc$k600.allCorr,
#        pch = 19,
#        col = "blue")
# legend("topright",
#        legend = c("gas uncorrected","gas background salt corrected"),
#        pch = 19,
#        col = c("black", "blue"))

# For rebuilding the package
# Don't forget to update the version in the description file!
setwd("~/GitHub/NEON-reaeration/reaRate")
devtools::document()
devtools::check() 
setwd("~/GitHub/NEON-reaeration/")
devtools::install("reaRate") 

