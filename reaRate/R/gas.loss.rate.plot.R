##############################################################################################
#' @title Loss rate plots for QAQC

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
#' @importFrom stats lm
#' @importFrom stats lsfit
#' @importFrom methods is

#' @param inputFile Name of the data frame containing the information needed to calculate the
#' reaeration parameters. If the headers are named: "injectionType", "eventID",
#' "stationToInjectionDistance", "plateauGasConc", "corrPlatSaltConc", "hoboSampleID",
#' "wettedWidth", respectively, no other inputs are required. Otherwise, the names of the
#' columns need to be input for the function to work. [string]
#' @param injectionTypeName Either constant rate or slug [string]
#' @param eventID A string identifier to link records collected as part of the same experiment,
#' SITE.YYYYMMDD for NEON [string]
#' @param stationToInjectionDistance Dataframe column name for distance from station to
#' injection [string]
#' @param plateauGasConc Dataframe column name for natural log of gas concentration normalized to
#' background corrected salt concentration [string]
#' @param corrPlatSaltConc Dataframe column name for natural log of gas concentration normalized to
#' background corrected salt concentration [string]
#' @param savePlotPath If a user specifies a path the plots will be saved to this location [string]

#' @return This function returns a list of two dataframes, the input dataframe of data for up to
#' 4 stations per site per date and an output dataframe appended with loss rate, travel time,
#' SF6 reaeration rate, O2 gas transfer velocity, and Schmidt number 600 for a given site and date

#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

#' @keywords surface water, streams, rivers, reaeration, deaeration, SF6, metabolism, tracer

#' @examples
#' #where the data frame "reaFormatted" is already read in
#' #reaRatesCalc <- def.calc.reaeration(inputFile = reaFormatted,
#' #dataDir = paste(path.package("reaRate"),"inst\\extdata", sep = "\\"), plot = TRUE)
#' #where the data is read in from a file in the working directory (also works with a full path)
#' #reaRatesCalc <- def.calc.reaeration(inputFile =
#' #system.file("extdata", "reaTestData.csv", package = "reaRate"))

#' @seealso def.calc.peakTime for calculating travel times and def.format.reaeration for
#' formatting reaeration data

#' @export

# changelog and author contributions / copyrights
#   Kaelin M. Cawley (2021-09-06)
#     original creation
#   Kaelin M. Cawley (2022-02-11)
#     updated to work with cleaned data columns
##############################################################################################
#This code is for calculating reaeration rates and Schmidt numbers
gas.loss.rate.plot <- function(
  inputFile = NULL,
  injectionTypeName = "injectionType",
  eventID = "eventID",
  stationToInjectionDistance = "stationToInjectionDistance",
  plateauGasConc = "plateauGasConc",
  corrPlatSaltConc = "corrPlatSaltConc",
  savePlotPath = NULL
){

  #Remove all model injections since they don't have gas tracer data
  inputFile <- inputFile[inputFile$injectionType %in% c("NaCl", "NaBr"),]
  inputFile$slopeRaw <- NA
  inputFile$slopeClean <- NA
  inputFile$slopeSaltCorr <- NA
  inputFile$slopeBackCorr <- NA

  for(currEventID in unique(inputFile$eventID)){

    #Other handy things
    currQ <- unique(inputFile$fieldDischarge[inputFile$eventID == currEventID])

    #Background correct salt samples, normalize gas concentration, and natural log transform the plateau gas concentrations
    stations <- inputFile$namedLocation[inputFile$eventID == currEventID]

    backSaltY <- inputFile$backgroundSaltConc[inputFile$eventID == currEventID]
    backSaltY[is.na(backSaltY)] <- 0
    backSaltX <- inputFile$stationToInjectionDistance[inputFile$eventID == currEventID]

    platSalt <- as.character(inputFile$plateauSaltConc[inputFile$eventID == currEventID])
    platSaltY <- NULL
    platSaltX <- NULL
    platSaltYCorr <- NULL
    
    mixStatus <- as.character(inputFile$unmixedStationFlag[inputFile$eventID == currEventID])
    saltColClean <- NULL
    gasColClean <- NULL
    
    platSaltClean <- as.character(inputFile$plateauSaltConcClean[inputFile$eventID == currEventID])
    platSaltYClean <- NULL
    platSaltXClean <- NULL
    platSaltYCorrClean <- NULL

    platGas <- as.character(inputFile$plateauGasConc[inputFile$eventID == currEventID])
    platGasY <- NULL
    platGasX <- NULL
    platGasSaltCorrY <- NULL
    platGasBackSaltCorrY <- NULL
    
    platGasClean <- as.character(inputFile$plateauGasConcClean[inputFile$eventID == currEventID])
    platGasYClean <- NULL
    platGasXClean <- NULL
    platGasSaltCorrYClean <- NULL
    platGasBackSaltCorrYClean <- NULL
    
    

    #Create the data for plotting
    for(currStat in seq(along = stations)){
      platSaltToAdd <- as.numeric(strsplit(platSalt[currStat],"\\|")[[1]])
      statDistToAddSalt <- rep(backSaltX[currStat],length(platSaltToAdd))
      corrSaltToAdd <- platSaltToAdd - backSaltY[currStat]
      platSaltY <- c(platSaltY, platSaltToAdd)
      platSaltX <- c(platSaltX, statDistToAddSalt)
      platSaltYCorr <- c(platSaltYCorr, corrSaltToAdd)
      
      platSaltToAddClean <- as.numeric(strsplit(platSaltClean[currStat],"\\|")[[1]])
      statDistToAddSaltClean <- rep(backSaltX[currStat],length(platSaltToAddClean))
      corrSaltToAddClean <- platSaltToAddClean - backSaltY[currStat]
      platSaltYClean <- c(platSaltYClean, platSaltToAddClean)
      platSaltXClean <- c(platSaltXClean, statDistToAddSaltClean)
      platSaltYCorrClean <- c(platSaltYCorrClean, corrSaltToAddClean)

      platGasToAdd <- as.numeric(strsplit(platGas[currStat],"\\|")[[1]])
      statDistToAddGas <- rep(backSaltX[currStat],length(platGasToAdd))
      platGasY <- c(platGasY, platGasToAdd)
      platGasX <- c(platGasX, statDistToAddGas)
      
      platGasToAddClean <- as.numeric(strsplit(platGasClean[currStat],"\\|")[[1]])
      statDistToAddGasClean <- rep(backSaltX[currStat],length(platGasToAddClean))
      platGasYClean <- c(platGasYClean, platGasToAddClean)
      platGasXClean <- c(platGasXClean, statDistToAddGasClean)
      
      platGasSaltCorrToAdd <- platGasToAddClean/mean(platSaltToAddClean, na.rm = TRUE)
      platGasBackSaltCorrToAdd <- platGasToAddClean/(mean(platSaltToAddClean, na.rm = TRUE) - backSaltY[currStat])
      platGasSaltCorrY <- c(platGasSaltCorrY, platGasSaltCorrToAdd)
      platGasBackSaltCorrY <- c(platGasBackSaltCorrY, platGasBackSaltCorrToAdd)
      
      # if(grepl("unmixedSaltStation",mixStatus[currStat])){
      #   saltCol <- c(saltCol,rep("grey21",length(platSaltToAdd)))
      # }else{
      #   saltCol <- c(saltCol,rep("grey61",length(platSaltToAdd)))
      # }
      
      if(grepl("unmixedSaltStation",mixStatus[currStat])){
        saltColClean <- c(saltColClean,rep("darkblue",length(platSaltToAddClean)))
      }else{
        saltColClean <- c(saltColClean,rep("deepskyblue3",length(platSaltToAddClean)))
      }
      
      if(grepl("unmixedGasStation",mixStatus[currStat])){
        gasColClean <- c(gasColClean,rep("purple4",length(platGasToAddClean)))
      }else{
        gasColClean <- c(gasColClean,rep("purple",length(platGasToAddClean)))
      }
    }

    #Skip the site if there isn't any salt data
    if(length(platSaltY) < 1 | length(platGasY) < 1){
      next
    }

    #Plot the salt background and plateau data
    # backRange <- c(min(backSaltY, na.rm = TRUE),
    #                max(backSaltY, na.rm = TRUE) * 5)
    # minPlat <- min(platSaltYCorr,platSaltY, na.rm = TRUE)
    # maxPlat <- max(platSaltYCorr,platSaltY, na.rm = TRUE)
    saltRange <- c(min(backSaltY, platSaltY, platSaltYCorr, na.rm = TRUE),
                   max(backSaltY, platSaltY, platSaltYCorr, na.rm = TRUE))
    #It just doesn't seem to look right in the rstudio plots
    # par(mar = c(4, 4, 7, 4) + 0.55)
    # plot(backSaltX,
    #      backSaltY,
    #      xlab = "Distance from Injection (m)",
    #      ylab = "Background Salt Concentration (ppm)",
    #      main = paste0(currEventID, " (", currQ, " lps)"),
    #      ylim = saltRange,
    #      pch = 19,
    #      col = c("magenta"))
    # par(new = TRUE)
    # plot(platSaltX,
    #      platSaltY,
    #      axes = FALSE,
    #      xlab = "",
    #      ylab = "",
    #      col = "red",
    #      pch = 4,
    #      cex = 1.5,
    #      ylim = saltRange)
    # axis(side = 4, at = pretty(saltRange))
    # mtext("Plateau Salt Concentrations (ppm)", side = 4, line = 3)
    # points(platSaltXClean,platSaltYClean,col = saltColClean, pch = 18, cex = 2)
    # points(platSaltXClean, platSaltYCorrClean, col = "green", pch = 19)
    # 
    # #Add lengend at the top
    # graphics::legend(x = "topright",
    #                  inset = c(-0.17,-0.28),
    #                  xpd = TRUE,
    #                  legend = c("Background","Plateau - all","Plateau Clean - mixed","Plateau Clean - unmixed","Corrected"),
    #                  col = c("magenta","red","deepskyblue3","darkblue","green"),
    #                  pch = c(19,4,18,18,19),
    #                  cex = c(0.9,0.9,0.9,0.9,0.9))

    if(!is.null(savePlotPath)){
      png(paste0(savePlotPath,"/salt_",gsub("\\.","_",currEventID),".png"))
      par(mar = c(4, 4, 7, 4) + 0.55)
      plot(backSaltX,
           backSaltY,
           xlab = "Distance from Injection (m)",
           ylab = "Background Salt Concentration (ppm)",
           main = paste0(currEventID, " (", currQ, " lps)"),
           ylim = saltRange,
           pch = 19,
           col = c("magenta"))
      par(new = TRUE)
      plot(platSaltX,
           platSaltY,
           axes = FALSE,
           xlab = "",
           ylab = "",
           col = "red",
           pch = 4,
           cex = 1.5,
           ylim = saltRange)
      axis(side = 4, at = pretty(saltRange))
      mtext("Plateau Salt Concentrations (ppm)", side = 4, line = 3)
      points(platSaltXClean,platSaltYClean,col = saltColClean, pch = 18, cex = 2)
      points(platSaltXClean, platSaltYCorrClean, col = "green", pch = 19)

      #Add lengend at the top
      graphics::legend(x = "topright",
                       inset = c(-0.17,-0.28),
                       xpd = TRUE,
                       legend = c("Background","Plateau - all","Plateau Clean - mixed","Plateau Clean - unmixed","Corrected"),
                       col = c("magenta","red","deepskyblue3","darkblue","green"),
                       pch = c(19,4,18,18,19),
                       cex = c(0.9,0.9,0.9,0.9,0.9))
      dev.off()
    }


    #Plot the gas data with different corrections
    logPlatGasY <- log(platGasY)
    platGasXforY <- platGasX[!is.infinite(logPlatGasY)]
    logPlatGasY <- logPlatGasY[!is.infinite(logPlatGasY)]
    
    logPlatGasYClean <- log(platGasYClean)
    platGasXforYClean <- platGasXClean[!is.infinite(logPlatGasYClean)]
    logPlatGasYClean <- logPlatGasYClean[!is.infinite(logPlatGasYClean)]

    logPlatGasSaltCorrY <- log(platGasSaltCorrY)
    platGasXforSaltCorrY <- platGasXClean[!is.infinite(logPlatGasSaltCorrY)]
    logPlatGasSaltCorrY <- logPlatGasSaltCorrY[!is.infinite(logPlatGasSaltCorrY)]
    
    logPlatGasBackSaltCorrY <- log(platGasBackSaltCorrY)
    platGasXforBackSaltCorrY <- platGasXClean[!is.infinite(logPlatGasBackSaltCorrY)]
    logPlatGasBackSaltCorrY <- logPlatGasBackSaltCorrY[!is.infinite(logPlatGasBackSaltCorrY)]

    gasRange <- c(min(logPlatGasY,logPlatGasYClean,logPlatGasSaltCorrY,logPlatGasBackSaltCorrY, na.rm = TRUE),
                  max(logPlatGasY,logPlatGasYClean,logPlatGasSaltCorrY,logPlatGasBackSaltCorrY, na.rm = TRUE))
    # par(mar = c(7, 4, 4, 4) + 0.3)
    # plot(platGasXforY,
    #      logPlatGasY,
    #      xlab = "Distance from Injection (m)",
    #      ylab = "Gas Concentration log(ppmv)",
    #      ylim = gasRange,
    #      main = paste0(currEventID, " (", currQ, " lps)"),
    #      pch = 4,
    #      col = "red")
    # par(new = TRUE)
    # plot(platGasXforSaltCorrY,
    #      logPlatGasSaltCorrY,
    #      axes = FALSE,
    #      xlab = "",
    #      ylab = "",
    #      ylim = gasRange,
    #      col = "cyan",
    #      pch = 19)
    # axis(side = 4, at = pretty(gasRange))
    # mtext("Salt Corrected Gas log(ppmv/ppm)", side = 4, line = 3)
    # points(platGasXforYClean, logPlatGasYClean, col = "purple", pch = 18)
    # points(platGasXforBackSaltCorrY, logPlatGasBackSaltCorrY, col = "orange", pch = 1)
    # 
    # try(rawLineFit <- lsfit(platGasX,logPlatGasY), silent = T)
    # rawSlope <- rawLineFit$coefficients[["X"]]
    # try(cleanLineFit <- lsfit(platGasXClean,logPlatGasYClean), silent = T)
    # cleanSlope <- cleanLineFit$coefficients[["X"]]
    # try(corrLineFit <- lsfit(platGasX,logPlatGasSaltCorrY), silent = T)
    # corrSlope <- corrLineFit$coefficients[["X"]]
    # try(backLineFit <- lsfit(platGasX,logPlatGasBackSaltCorrY), silent = T)
    # backSlope <- backLineFit$coefficients[["X"]]
    # 
    # abline(a = rawLineFit$coefficients[["Intercept"]], b = rawLineFit$coefficients[["X"]], col = "red")
    # abline(a = cleanLineFit$coefficients[["Intercept"]], b = cleanLineFit$coefficients[["X"]], col = "purple")
    # abline(a = corrLineFit$coefficients[["Intercept"]], b = corrLineFit$coefficients[["X"]], col = "cyan")
    # abline(a = backLineFit$coefficients[["Intercept"]], b = backLineFit$coefficients[["X"]], col = "orange")
    # 
    # #Add legend at the top
    # graphics::legend(x = "bottom",
    #        inset = c(-0.1,-0.35),
    #        xpd = TRUE,
    #        legend = c(paste0("Raw Gas ", signif(rawLineFit$coefficients[["X"]], 3)),
    #                   paste0("Clean gas ", signif(corrLineFit$coefficients[["X"]], 3)),
    #                   paste0("Plateau Corrected ", signif(corrLineFit$coefficients[["X"]], 3)),
    #                   paste0("Background Corrected ", signif(backLineFit$coefficients[["X"]], 3))),
    #        col = c("red","purple","cyan","orange"),
    #        pch = c(4,19,18,1),
    #        cex = c(0.9,0.9,0.9,0.9),
    #        horiz = TRUE)

    if(!is.null(savePlotPath)){
      png(paste0(savePlotPath,"/gas_",gsub("\\.","_",currEventID),".png"))
      par(mar = c(4, 4, 5, 9) + 0.3)
      plot(platGasXforY,
           logPlatGasY,
           xlab = "Distance from Injection (m)",
           ylab = "Gas Concentration log(ppmv)",
           ylim = gasRange,
           main = paste0(currEventID, " (", currQ, " lps)"),
           pch = 4,
           col = "red")
      par(new = TRUE)
      plot(platGasXforSaltCorrY,
           logPlatGasSaltCorrY,
           axes = FALSE,
           xlab = "",
           ylab = "",
           ylim = gasRange,
           col = "cyan",
           pch = 19)
      axis(side = 4, at = pretty(gasRange))
      mtext("Salt Corrected Gas log(ppmv/ppm)", side = 4, line = 3)
      points(platGasXforYClean, logPlatGasYClean, col = "purple", pch = 18)
      points(platGasXforBackSaltCorrY, logPlatGasBackSaltCorrY, col = "orange", pch = 1)
      
      try(rawLineFit <- lsfit(platGasX,logPlatGasY), silent = T)
      rawSlope <- rawLineFit$coefficients[["X"]]
      try(cleanLineFit <- lsfit(platGasXClean,logPlatGasYClean), silent = T)
      cleanSlope <- cleanLineFit$coefficients[["X"]]
      try(corrLineFit <- lsfit(platGasXClean,logPlatGasSaltCorrY), silent = T)
      corrSlope <- corrLineFit$coefficients[["X"]]
      try(backLineFit <- lsfit(platGasXClean,logPlatGasBackSaltCorrY), silent = T)
      backSlope <- backLineFit$coefficients[["X"]]
      
      abline(a = rawLineFit$coefficients[["Intercept"]], b = rawLineFit$coefficients[["X"]], col = "red")
      abline(a = cleanLineFit$coefficients[["Intercept"]], b = cleanLineFit$coefficients[["X"]], col = "purple")
      abline(a = corrLineFit$coefficients[["Intercept"]], b = corrLineFit$coefficients[["X"]], col = "cyan")
      abline(a = backLineFit$coefficients[["Intercept"]], b = backLineFit$coefficients[["X"]], col = "orange")
      
      #Add legend at the top
      graphics::legend(x = "topright",
                       inset = c(-0.45,-0.2),
                       xpd = TRUE,
                       legend = c(paste0("Raw Gas ", signif(rawLineFit$coefficients[["X"]], 3)),
                                  paste0("Clean gas ", signif(corrLineFit$coefficients[["X"]], 3)),
                                  paste0("Plateau Corrected ", signif(corrLineFit$coefficients[["X"]], 3)),
                                  paste0("Background Corrected ", signif(backLineFit$coefficients[["X"]], 3))),
                       col = c("red","purple","cyan","orange"),
                       pch = c(4,19,18,1),
                       cex = c(0.9,0.9,0.9,0.9))
      dev.off()
    }

    inputFile$slopeRaw[inputFile$eventID == currEventID] <- rawLineFit$coefficients[["X"]]
    inputFile$slopeClean[inputFile$eventID == currEventID] <- cleanLineFit$coefficients[["X"]]
    inputFile$slopeSaltCorr[inputFile$eventID == currEventID] <- corrLineFit$coefficients[["X"]]
    inputFile$slopeBackCorr[inputFile$eventID == currEventID] <- backLineFit$coefficients[["X"]]
  }
  return(inputFile)
}
