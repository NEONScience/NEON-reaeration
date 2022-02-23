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

#' @examples
#' #where the data .zip file is in the working directory and has the default name,
#' #reaFormatted <- def.format.reaeration()

#' @seealso def.calc.tracerTime.R for calculating the stream travel time,
#' def.plot.reaQcurve.R for plotting reaeration rate versusu stream flow,
#' def.format.reaeration for formatting reaeration data

# changelog and author contributions / copyrights
#   Kaelin M. Cawley (2021-02-04)
#     original creation
##############################################################################################

#User Inputs
siteID <- "CARI"
#Take a look at REDB next to see if data looks compromised from broken diffuser
#siteID <- "CARI" #Not an ADCP site (yet!), NaBr injection
#siteID <- "MAYF" #ADCP site for testing
#siteID <- "WALK" #When station 2 is used instead of station 1
#siteID <- "BLUE" #For testing to get travel-time for model experiments

#String constants
reaDPID <- "DP1.20190.001"
dscDPID <- "DP1.20048.001"
wqDPID <- "DP1.20288.001"

# Download Reaeration Data
reaInputList <- neonUtilities::loadByProduct(dpID = reaDPID, site = siteID, check.size = FALSE) #, startdate = "2019-01-01", enddate = "2021-12-01"

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

rea_backgroundFieldCondData <- rea_backgroundFieldCondDataIn
rea_backgroundFieldSaltData <- rea_backgroundFieldSaltDataIn
rea_fieldData <- rea_fieldDataIn
rea_plateauMeasurementFieldData <- rea_plateauMeasurementFieldDataIn
rea_plateauSampleFieldData <- rea_plateauSampleFieldDataIn
rea_externalLabDataSalt <- rea_externalLabDataSaltIn
rea_externalLabDataGas <- rea_externalLabDataGasIn
rea_widthFieldData <- rea_widthFieldDataIn
dsc_fieldData <- dsc_fieldDataIn
dsc_individualFieldData <- dsc_individualFieldDataIn
dsc_fieldDataADCP <- dsc_fieldDataADCPIn
waq_instantaneous <- waq_instantaneousIn

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

inputFile = reaFormatted
injectionTypeName = "injectionType"
eventID = "eventID"
stationToInjectionDistance = "stationToInjectionDistance"
plateauGasConc = "plateauGasConc"
corrPlatSaltConc = "corrPlatSaltConc"
savePlotPath = "C:/Users/kcawley/Desktop/reaTesting/GUIL_20220211"

plotsOut <- reaRate::gas.loss.rate.plot(inputFile = reaFormatted,
                                        injectionTypeName = "injectionType",
                                        eventID = "eventID",
                                        stationToInjectionDistance = "stationToInjectionDistance",
                                        plateauGasConc = "plateauGasConc",
                                        corrPlatSaltConc = "corrPlatSaltConc",
                                        savePlotPath = "C:/Users/kcawley/Desktop/reaTesting/GUIL_20220211")

inputFile = reaFormatted
loggerData = reaInputList$rea_conductivityFieldData
namedLocation = "namedLocation"
injectionTypeName = "injectionType"
eventID = "eventID"
stationToInjectionDistance = "stationToInjectionDistance"
plateauGasConc = "plateauGasConc"
corrPlatSaltConc = "corrPlatSaltConc"
hoboSampleID = "hoboSampleID"
discharge = "fieldDischarge"
waterTemp = "waterTemp"
wettedWidth = "wettedWidth"
plot = TRUE
savePlotPath = "C:/Users/kcawley/Desktop/reaTesting/"
processingInfo = NULL

reaRatesCalc <- reaRate::def.calc.reaeration(inputFile = reaFormatted,
                                             loggerData = reaInputList$rea_conductivityFieldData,
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
                                             processingInfo = NULL)

outputDF <- reaRatesCalc$outputDF
outputDFClean <- outputDF[outputDF$k600 > 0,]
plot(outputDFClean$meanQ, outputDFClean$k600)

#Testing when not using loadByProduct
rea_backgroundFieldCondDataIn <- read.table("C:/Users/kcawley/Desktop/NEON.D12.BLDE.DP1.20190.001/NEON.D12.BLDE.DP1.20190.001.rea_backgroundFieldCondData.2018-06.expanded.20201218T123834Z.csv",
                                            sep = ",",
                                            header = TRUE)
rea_backgroundFieldSaltDataIn <- read.table("C:/Users/kcawley/Desktop/NEON.D12.BLDE.DP1.20190.001/NEON.D12.BLDE.DP1.20190.001.rea_backgroundFieldSaltData.2018-06.expanded.20201218T123834Z.csv",
                                            sep = ",",
                                            header = TRUE)
rea_fieldDataIn <- read.table("C:/Users/kcawley/Desktop/NEON.D12.BLDE.DP1.20190.001/NEON.D12.BLDE.DP1.20190.001.rea_fieldData.2018-06.expanded.20201218T123834Z.csv",
                              sep = ",",
                              header = TRUE)
rea_plateauMeasurementFieldDataIn <- read.table("C:/Users/kcawley/Desktop/NEON.D12.BLDE.DP1.20190.001/NEON.D12.BLDE.DP1.20190.001.rea_plateauMeasurementFieldData.2018-06.expanded.20201218T123834Z.csv",
                                                sep = ",",
                                                header = TRUE)
rea_externalLabDataSaltIn <- read.table("C:/Users/kcawley/Desktop/NEON.D12.BLDE.DP1.20190.001/NEON.D12.BLDE.DP1.20190.001.rea_externalLabDataSalt.2018-06.expanded.20201218T123834Z.csv",
                                        sep = ",",
                                        header = TRUE)
rea_externalLabDataGasIn <- read.table("C:/Users/kcawley/Desktop/NEON.D12.BLDE.DP1.20190.001/NEON.D12.BLDE.DP1.20190.001.rea_externalLabDataGas.2018-06.expanded.20201218T123834Z.csv",
                                       sep = ",",
                                       header = TRUE)
rea_widthFieldDataIn <- read.table("C:/Users/kcawley/Desktop/NEON.D12.BLDE.DP1.20190.001/NEON.D12.BLDE.DP1.20190.001.rea_widthFieldData.2018-06.expanded.20201218T123834Z.csv",
                                   sep = ",",
                                   header = TRUE)
#Discharge data download
dsc_fieldDataIn <- read.table("C:/Users/kcawley/Desktop/NEON.D12.BLDE.DP1.20048.001/NEON.D12.BLDE.DP1.20048.001.dsc_fieldData.2018-06.basic.20210112T153513Z.csv",
                              sep = ",",
                              header = TRUE)
dsc_individualFieldDataIn <- read.table("C:/Users/kcawley/Desktop/NEON.D12.BLDE.DP1.20048.001/NEON.D12.BLDE.DP1.20048.001.dsc_individualFieldData.2018-06.basic.20210112T153513Z.csv",
                                        sep = ",",
                                        header = TRUE)

reaFormatted <- reaRate::def.format.reaeration(rea_backgroundFieldCondData = rea_backgroundFieldCondDataIn,
                                               rea_backgroundFieldSaltData = rea_backgroundFieldSaltDataIn,
                                               rea_fieldData = rea_fieldDataIn,
                                               rea_plateauMeasurementFieldData = rea_plateauMeasurementFieldDataIn,
                                               rea_externalLabDataSalt = rea_externalLabDataSaltIn,
                                               rea_externalLabDataGas = rea_externalLabDataGasIn,
                                               rea_widthFieldData = rea_widthFieldDataIn,
                                               dsc_fieldData = dsc_fieldDataIn,
                                               dsc_individualFieldData = dsc_individualFieldDataIn)



#Testing some ideas with plotly

library(plotly)
library(crosstalk)
library(DT)


sd <- SharedData$new(iris)

a <- plot_ly(sd, x = ~Sepal.Width, y = ~Sepal.Length) %>% 
  add_markers(alpha = 0.5) %>%
  highlight("plotly_selected", dynamic = TRUE)


options(persistent = TRUE)

p <- datatable(sd)

bscols(widths = c(6, 4), a, p)

test <- as.data.frame(p)


library(plotly)
library(shiny)

mtcars$key <- row.names(mtcars)
mtcars$col <- "black"

ui <- fluidPage(
  plotlyOutput("plot")
)

server <- function(input, output, session) {
  output$plot <- renderPlotly({
    click_data <- event_data("plotly_click", priority   = "event")
    print(click_data)
    select_data <- event_data("plotly_selected", priority   = "event")
    if (!is.null(select_data)) {
      mtcars[mtcars$key %in% select_data$customdata, "col"] <- "blue"
    }
    if (!is.null(click_data)) {
      mtcars[mtcars$key %in% click_data$customdata, "col"] <- "red"
    }
    p <- plot_ly(mtcars, x = ~mpg, y=~wt, colors = ~sort(unique(col)), color = ~col, customdata = ~key, type = "scatter", mode = "markers") %>% layout(dragmode = "lasso")
  })
}

shinyApp(ui, server)
