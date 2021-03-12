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
#siteID <- "BLDE" #Not an ADCP site (yet!), NaBr injection
siteID <- "MAYF" #ADCP site for testing

#String constants
reaDPID <- "DP1.20190.001"
dscDPID <- "DP1.20048.001"

# Download Reaeration Data
reaInputList <- neonUtilities::loadByProduct(dpID = reaDPID, site = siteID, check.size = FALSE)

rea_backgroundFieldCondDataIn <- reaInputList$rea_backgroundFieldCondData
rea_backgroundFieldSaltDataIn <- reaInputList$rea_backgroundFieldSaltData
rea_fieldDataIn <- reaInputList$rea_fieldData
rea_plateauMeasurementFieldDataIn <- reaInputList$rea_plateauMeasurementFieldData
rea_externalLabDataSaltIn <- reaInputList$rea_externalLabDataSalt
rea_externalLabDataGasIn <- reaInputList$rea_externalLabDataGas
rea_widthFieldDataIn <- reaInputList$rea_widthFieldData

# Download Discharge Data
qInputList <- neonUtilities::loadByProduct(dpID = dscDPID, site = siteID, check.size = FALSE)

dsc_fieldDataIn <- qInputList$dsc_fieldData
dsc_individualFieldDataIn <- qInputList$dsc_individualFieldData
dsc_fieldDataADCPIn <- qInputList$dsc_fieldDataADCP

rea_backgroundFieldCondData <- rea_backgroundFieldCondDataIn
rea_backgroundFieldSaltData <- rea_backgroundFieldSaltDataIn
rea_fieldData <- rea_fieldDataIn
rea_plateauMeasurementFieldData <- rea_plateauMeasurementFieldDataIn
rea_externalLabDataSalt <- rea_externalLabDataSaltIn
rea_externalLabDataGas <- rea_externalLabDataGasIn
rea_widthFieldData <- rea_widthFieldDataIn
dsc_fieldData <- dsc_fieldDataIn
dsc_individualFieldData <- dsc_individualFieldDataIn
dsc_fieldDataADCP <- dsc_fieldDataADCPIn

reaFormatted <- reaRate::def.format.reaeration(rea_backgroundFieldCondData = rea_backgroundFieldCondDataIn,
                                      rea_backgroundFieldSaltData = rea_backgroundFieldSaltDataIn,
                                      rea_fieldData = rea_fieldDataIn,
                                      rea_plateauMeasurementFieldData = rea_plateauMeasurementFieldDataIn,
                                      rea_externalLabDataSalt = rea_externalLabDataSaltIn,
                                      rea_externalLabDataGas = rea_externalLabDataGasIn,
                                      rea_widthFieldData = rea_widthFieldDataIn,
                                      dsc_fieldData = dsc_fieldDataIn,
                                      dsc_individualFieldData = dsc_individualFieldDataIn,
                                      dsc_fieldDataADCP = dsc_fieldDataADCPIn)

# inputFile = reaFormatted
# loggerData = reaInputList$rea_conductivityFieldData
# namedLocation = "namedLocation"
# injectionTypeName = "injectionType"
# eventID = "eventID"
# stationToInjectionDistance = "stationToInjectionDistance"
# plateauGasConc = "plateauGasConc"
# corrPlatSaltConc = "corrPlatSaltConc"
# hoboSampleID = "hoboSampleID"
# discharge = "fieldDischarge"
# waterTemp = "waterTemp"
# wettedWidth = "wettedWidth"
# plot = TRUE
# savePlotPath = NULL
# processingInfo = NULL

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
