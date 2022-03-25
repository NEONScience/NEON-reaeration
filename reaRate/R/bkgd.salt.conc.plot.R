##############################################################################################
#' @title Background Salt Plot for QAQC Decisions

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
#' @importFrom grDevices colorRampPalette
#' @importFrom stats lm
#' @importFrom stats lsfit
#' @importFrom methods is

#' @param inputFile Name of the data frame containing the information needed to plot background 
#' salt information for QAQC. [string]
#' @param savePlotPath If a user specifies a path the plots will be saved to this location [string]

#' @return This function makes three plots, background salt concentrations color coded for 
#' date and discharge and a plot of sensor specific conductance 20 minutes before the 
#' injection time versus grab sample concentration.

#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

#' @keywords surface water, streams, rivers, reaeration, deaeration, SF6, metabolism, tracer

#' @examples
#' #TBD

#' @seealso def.calc.peakTime for calculating travel times and def.format.reaeration for
#' formatting reaeration data

#' @export

# changelog and author contributions / copyrights
#   Kaelin M. Cawley (2022-02-26)
#     original creation
##############################################################################################
#This code is for calculating reaeration rates and Schmidt numbers
bkgd.salt.conc.plot <- function(
  inputFile = NULL,
  savePlotPath = NULL
){

  if(length(unique(inputFile$siteID)) > 1){
    stop("This function is still in development and can only handle one site at a time in the input file.")
  }
  
  # Make it cyberpunk
  colfunc <- colorRampPalette(c("red","orange","yellow","green","blue","purple"))
  theColors <- colfunc(length(unique(inputFile$eventID)))[order(unique(inputFile$eventID))]
  
  # Set up the plot
  yMin <- min(inputFile$backgroundSaltConc, na.rm = TRUE)
  yMax <- max(inputFile$backgroundSaltConc, na.rm = TRUE)
  xMin <- min(inputFile$stationToInjectionDistance, na.rm = TRUE)
  xMax <- max(inputFile$stationToInjectionDistance, na.rm = TRUE)

  colIdx <- 1
  xData <- inputFile$stationToInjectionDistance[inputFile$eventID == inputFile$eventID[1]]
  yData <- inputFile$backgroundSaltConc[inputFile$eventID == inputFile$eventID[1]]
  plot(xData[order(xData)],
       yData[order(xData)],
       xlab = "Distance from Injection (m)",
       ylab = "Background Salt Concentration (ppm)",
       main = paste0(unique(inputFile$siteID), " - Date"),
       ylim = c(-1, yMax),
       xlim = c(xMin, xMax),
       col = theColors[colIdx],
       type = "b",
       pch = 19)
  
  for(currEventID in unique(inputFile$eventID)[2:length(unique(inputFile$eventID))]){
    colIdx <- colIdx + 1
    xData <- inputFile$stationToInjectionDistance[inputFile$eventID == currEventID]
    yData <- inputFile$backgroundSaltConc[inputFile$eventID == currEventID]
    lines(xData[order(xData)],
           yData[order(xData)],
          col = theColors[colIdx])
    points(xData[order(xData)],
          yData[order(xData)],
          col = theColors[colIdx],
          pch = 19)
  }
  graphics::legend(x = "bottomright",
                   legend = c("Oldest","Most Recent"),
                   col = c("red","purple"),
                   pch = c(19, 19),
                   cex = c(0.9, 0.9))
  
  if(!is.null(savePlotPath)){
    dev.copy(png,paste0(savePlotPath,"/",unique(inputFile$siteID),"_backSaltByEventID.png"))
    dev.off()
  }
  
  #Plot related to discharge
  dscInpout <- inputFile <- inputFile[!is.na(inputFile$fieldDischarge),]
  theColors <- colfunc(length(unique(dscInpout$eventID)))[order(unique(dscInpout$eventID))]
  
  # Order eventIDs by discharge
  dscInpout <- dscInpout[order(dscInpout$fieldDischarge),]
  
  # Set up the plot
  yMin <- min(dscInpout$backgroundSaltConc, na.rm = TRUE)
  yMax <- max(dscInpout$backgroundSaltConc, na.rm = TRUE)
  xMin <- min(dscInpout$stationToInjectionDistance, na.rm = TRUE)
  xMax <- max(dscInpout$stationToInjectionDistance, na.rm = TRUE)
  
  colIdx <- 1
  xData <- dscInpout$stationToInjectionDistance[dscInpout$eventID == dscInpout$eventID[1]]
  yData <- dscInpout$backgroundSaltConc[dscInpout$eventID == dscInpout$eventID[1]]
  plot(xData[order(xData)],
       yData[order(xData)],
       xlab = "Distance from Injection (m)",
       ylab = "Background Salt Concentration (ppm)",
       main = paste0(unique(dscInpout$siteID), " - Discharge"),
       ylim = c(0, yMax),
       xlim = c(xMin, xMax),
       type = "b",
       pch = 19,
       col = theColors[colIdx])
  
  for(currEventID in unique(dscInpout$eventID)[2:length(unique(dscInpout$eventID))]){
    colIdx <- colIdx + 1
    xData <- dscInpout$stationToInjectionDistance[dscInpout$eventID == currEventID]
    yData <- dscInpout$backgroundSaltConc[dscInpout$eventID == currEventID]
    points(xData[order(xData)],
           yData[order(xData)],
           pch = 19,
           col = theColors[colIdx])
    lines(xData[order(xData)],
          yData[order(xData)],
           col = theColors[colIdx])
  }
  graphics::legend(x = "bottomright",
                   legend = c(paste0("Low Flow", " (",round(min(dscInpout$fieldDischarge, na.rm = TRUE), digits = 0)," lps)"),
                              paste0("High Flow", " (",round(max(dscInpout$fieldDischarge, na.rm = TRUE), digits = 0)," lps)")),
                   col = c("red","purple"),
                   pch = c(19, 19),
                   cex = c(0.9, 0.9))
  
  if(!is.null(savePlotPath)){
    dev.copy(png,paste0(savePlotPath,"/",unique(inputFile$siteID),"_backSaltByDischarge.png"))
    dev.off()
  }
  
  # Plot background concentration to sensor specific conductance
  if(unique(inputFile$siteID) == "GUIL"){ 
    inputFile$backgroundSensorCond[inputFile$backgroundSensorCond < 105] <- NA
  }
  thing3 <- plot(inputFile$backgroundSaltConc, 
       inputFile$backgroundSensorCond,
       ylab = "Sensor Specific Conductance (uS/cm at 25 C)",
       xlab = "Background Salt Concentration (ppm)",
       main = "Sensor versus grab sample background concentrations",
       col = "blue")
  try(sensorFit <- lsfit(inputFile$backgroundSaltConc, 
                         inputFile$backgroundSensorCond), silent = T)
  abline(a = sensorFit$coefficients[["Intercept"]], b = sensorFit$coefficients[["X"]], col = "blue")
  
  if(!is.null(savePlotPath)){
    dev.copy(png,paste0(savePlotPath,"/",unique(inputFile$siteID),"_sensorVersusGrab.png"))
    dev.off()
  }
  
}
