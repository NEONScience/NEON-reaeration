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

#' @param inputFile Name of the data frame containing the information needed to calculate the
#' reaeration parameters. If the headers are named: "injectionType", "eventID",
#' "stationToInjectionDistance", "plateauGasConc", "corrPlatSaltConc", "hoboSampleID",
#' "wettedWidth", respectively, no other inputs are required. Otherwise, the names of the
#' columns need to be input for the function to work. [string]
#' @param savePlotPath If a user specifies a path the plots will be saved to this location [string]

#' @return This function returns a list of two dataframes, the input dataframe of data for up to
#' 4 stations per site per date and an output dataframe appended with loss rate, travel time,
#' SF6 reaeration rate, O2 gas transfer velocity, and Schmidt number 600 for a given site and date

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
  colfunc <- colorRampPalette(c("cyan","deeppink"))
  theColors <- colfunc(length(unique(inputFile$eventID)))[order(unique(inputFile$eventID))]
  
  # Set up the plot
  yMin <- min(inputFile$backgroundSaltConc, na.rm = TRUE)
  yMax <- max(inputFile$backgroundSaltConc, na.rm = TRUE)
  xMin <- min(inputFile$stationToInjectionDistance, na.rm = TRUE)
  xMax <- max(inputFile$stationToInjectionDistance, na.rm = TRUE)

  colIdx <- 1
  plot(inputFile$stationToInjectionDistance[inputFile$eventID == inputFile$eventID[1]],
       inputFile$backgroundSaltConc[inputFile$eventID == inputFile$eventID[1]],
       xlab = "Distance from Injection (m)",
       ylab = "Background Salt Concentration (ppm)",
       main = paste0(unique(inputFile$siteID)),
       ylim = c(yMin, yMax),
       xlim = c(xMin, xMax),
       col = theColors[colIdx])
  
  for(currEventID in unique(inputFile$eventID)[2:length(unique(inputFile$eventID))]){
    colIdx <- colIdx + 1
    points(inputFile$stationToInjectionDistance[inputFile$eventID == currEventID],
          inputFile$backgroundSaltConc[inputFile$eventID == currEventID],
          col = theColors[colIdx])
  }
  graphics::legend(x = "bottomright",
                   legend = c("Oldest","Most Recent"),
                   col = c("cyan","magenta"),
                   pch = c(19, 19),
                   cex = c(0.9, 0.9))
  
  #Plot related to discharge
  dscInpout <- inputFile <- inputFile[!is.na(inputFile$fieldDischarge),]
  # Make it cyberpunk
  colfunc <- colorRampPalette(c("cyan","deeppink"))
  theColors <- colfunc(length(unique(dscInpout$eventID)))[order(unique(dscInpout$eventID))]
  
  # Order eventIDs by discharge
  dscInpout <- dscInpout[order(dscInpout$fieldDischarge),]
  
  # Set up the plot
  yMin <- min(dscInpout$backgroundSaltConc, na.rm = TRUE)
  yMax <- max(dscInpout$backgroundSaltConc, na.rm = TRUE)
  xMin <- min(dscInpout$stationToInjectionDistance, na.rm = TRUE)
  xMax <- max(dscInpout$stationToInjectionDistance, na.rm = TRUE)
  
  colIdx <- 1
  plot(dscInpout$stationToInjectionDistance[dscInpout$eventID == dscInpout$eventID[1]],
       dscInpout$backgroundSaltConc[dscInpout$eventID == dscInpout$eventID[1]],
       xlab = "Distance from Injection (m)",
       ylab = "Background Salt Concentration (ppm)",
       main = paste0(unique(dscInpout$siteID)),
       ylim = c(yMin, yMax),
       xlim = c(xMin, xMax),
       pch = 19,
       col = theColors[colIdx])
  
  for(currEventID in unique(dscInpout$eventID)[2:length(unique(dscInpout$eventID))]){
    colIdx <- colIdx + 1
    points(dscInpout$stationToInjectionDistance[dscInpout$eventID == currEventID],
           dscInpout$backgroundSaltConc[dscInpout$eventID == currEventID],
           col = theColors[colIdx])
  }
  graphics::legend(x = "bottomright",
                   legend = c("Low Flow","High Flow"),
                   col = c("cyan","magenta"),
                   pch = c(19, 19),
                   cex = c(0.9, 0.9))
  
  # Plot background concentration to sensor specific conductance
  if(unique(inputFile$siteID) == "GUIL"){ 
    inputFile$backgroundSensorCond[inputFile$backgroundSensorCond < 105] <- NA
  }
  plot(inputFile$backgroundSaltConc, 
       inputFile$backgroundSensorCond,
       ylab = "Sensor Specific Conductance (uS/cm at 25 C)",
       xlab = "Background Salt Concentration (ppm)",
       main = "Sensor versus grab sample background concentrations",
       col = "blue")
  try(sensorFit <- lsfit(inputFile$backgroundSaltConc, 
                         inputFile$backgroundSensorCond), silent = T)
  abline(a = sensorFit$coefficients[["Intercept"]], b = sensorFit$coefficients[["X"]], col = "blue")
  
  return(inputFile)
}
