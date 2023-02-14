##############################################################################################
#' @title MART Multilevel model for estimate gas exchange in streams

#' @author
#' Kaelin M. Cawley \email{kcawley@battelleecology.org} \cr
#' Robert O. Hall \cr

#' @description This script retrieves, joins, manually cleans, and applies a stan model to 
#' determine gas exchange rates for NEON data.

#' @references
#' License: CC0 1.0 Universal

#' @keywords surface water, streams, gas exchange, reaeration, deaeration, metabolism, tracer

# changelog and author contributions / copyrights
#   Kaelin M. Cawley (2023-02-14)
#     original creation
#########################################################################################

#User Inputs
siteID <- "MART"
workingPath <- "SET TO YOUR PATH/"
stanPlots <- paste0(workingPath,siteID)
plotPathLossRates <- paste0(workingPath,siteID,"/saltAndGasLossPlots")
plotPathBackground <- paste0(workingPath,siteID,"/backgroundPlots")
plotPathTravelTime <- paste0(workingPath,siteID,"/travelTimePlots")

#String constants
reaDPID <- "DP1.20190.001"
dscDPID <- "DP1.20048.001"
wqDPID <- "DP1.20288.001"

setwd(stanPlots)

# Download Reaeration Data
reaInputList <- neonUtilities::loadByProduct(dpID = reaDPID, site = siteID, check.size = FALSE)

# Download Discharge Data
qInputList <- neonUtilities::loadByProduct(dpID = dscDPID, site = siteID, check.size = FALSE)

# Download Sensor Data
sensorData <- neonUtilities::loadByProduct(dpID = wqDPID,
                                           site = siteID,
                                           check.size = FALSE)

save.image(file = paste0(stanPlots,"/data_", Sys.Date(),".RData"))
# load(paste0(stanPlots,"/data_2022-03-08.RData"))

# Format all of the tables into one
reaFormatted <- reaRate::def.format.reaeration(rea_backgroundFieldCondData = reaInputList$rea_backgroundFieldCondData,
                                               rea_backgroundFieldSaltData = reaInputList$rea_backgroundFieldSaltData,
                                               rea_fieldData = reaInputList$rea_fieldData,
                                               rea_plateauMeasurementFieldData = reaInputList$rea_plateauMeasurementFieldData,
                                               rea_plateauSampleFieldData = reaInputList$rea_plateauSampleFieldData,
                                               rea_externalLabDataSalt = reaInputList$rea_externalLabDataSalt,
                                               rea_externalLabDataGas = reaInputList$rea_externalLabDataGas,
                                               rea_widthFieldData = reaInputList$rea_widthFieldData,
                                               dsc_fieldData = qInputList$dsc_fieldData,
                                               dsc_individualFieldData = qInputList$dsc_individualFieldData,
                                               dsc_fieldDataADCP = qInputList$dsc_fieldDataADCP,
                                               waq_instantaneous = sensorData$waq_instantaneous)

# Write out csv of reaFormatted
write.csv(reaFormatted, 
          paste0(siteID, "_reaFormatted.csv"),
          row.names = FALSE)
reaFormatted <- read.csv(paste0(siteID, "_reaFormatted.csv"))

# Plot the salt and gas loss over distance and calculate slopes
plotsOut <- reaRate::gas.loss.rate.plot(inputFile = reaFormatted,
                                        injectionTypeName = "injectionType",
                                        eventID = "eventID",
                                        stationToInjectionDistance = "stationToInjectionDistance",
                                        plateauGasConc = "plateauGasConc",
                                        corrPlatSaltConc = "corrPlatSaltConc",
                                        savePlotPath = plotPathLossRates)

reaRate::bkgd.salt.conc.plot(inputFile = plotsOut,
                             savePlotPath = plotPathBackground)

reaRatesTrvlTime <- reaRate::def.calc.trvl.time(inputFile = plotsOut,
                                                loggerData = reaInputList$rea_conductivityFieldData,
                                                plot = TRUE,
                                                savePlotPath = plotPathTravelTime)

write.csv(reaRatesTrvlTime$outputDF, 
          paste0(siteID, "_outputDF.csv"),
          row.names = FALSE)



##############################################################################################
########## Prep and format things for stan and manual method ##########
##############################################################################################
#Starting from here for re-loading data to use new stan model as of 2/26/22
siteID <- "MART"
localPath <- "SET TO YOUR PATH/"
plotPathLossRates <- paste0(localPath,siteID,"/saltAndGasLossPlots")
plotPathBackground <- paste0(localPath,siteID,"/backgroundPlots")
plotPathTravelTime <- paste0(localPath,siteID,"/travelTimePlots")
stanPlots <- paste0(localPath,siteID)
load(paste0(localPath,siteID,"/MART_workspace.RData"))

plotsOut <- reaRate::gas.loss.rate.plot(inputFile = reaFormatted,
                                        injectionTypeName = "injectionType",
                                        eventID = "eventID",
                                        stationToInjectionDistance = "stationToInjectionDistance",
                                        plateauGasConc = "plateauGasConc",
                                        corrPlatSaltConc = "corrPlatSaltConc",
                                        savePlotPath = plotPathLossRates)

#Combine travel time info and new slopes
reaRateManual <- reaRatesTrvlTime$outputDF
for(currEvt in reaRateManual$eventID){
  try(reaRateManual$slopeNormClean[reaRateManual$eventID == currEvt] <- unique(plotsOut$slopeNormClean[plotsOut$eventID == currEvt]))
  try(reaRateManual$slopeNormCorr[reaRateManual$eventID == currEvt] <- unique(plotsOut$slopeNormCorr[plotsOut$eventID == currEvt]))
}


###Plots for the TWG Meeting
library(plotly)

###### k600 and K600 plots #####

`%nin%` <- Negate(`%in%`)

badExperiments.clean <- c("")

badExperiments.salt <- c("MART.20170228","MART.20170822","MART.20171204","MART.20180215","MART.20181228","MART.20190103","MART.20190226","MART.20201026","MART.20210120")

reaRatesCalc.normClean <- reaRate::def.calc.reaeration(inputFile = reaRateManual,
                                                       lossRateSF6 = "slopeNormClean",
                                                       peakMaxVelocity = "peakMaxVelocity",
                                                       meanDepth = "meanDepth",
                                                       meanTemp = "meanTemp",
                                                       outputSuffix = "normClean")

reaRatesCalc.normClean$OK[reaRatesCalc.normClean$eventID %nin% badExperiments.clean] <- "yes"
reaRatesCalc.normClean$OK[reaRatesCalc.normClean$eventID %in% badExperiments.clean] <- "no"

reaRatesCalc.normCorr <- reaRate::def.calc.reaeration(inputFile = reaRateManual,
                                                      lossRateSF6 = "slopeNormCorr",
                                                      peakMaxVelocity = "peakMaxVelocity",
                                                      meanDepth = "meanDepth",
                                                      meanTemp = "meanTemp",
                                                      outputSuffix = "normCorr")

reaRatesCalc.normCorr$OK[reaRatesCalc.normCorr$eventID %nin% badExperiments.salt] <- "yes"
reaRatesCalc.normCorr$OK[reaRatesCalc.normCorr$eventID %in% badExperiments.salt] <- "no"

reaRates <- merge(reaRatesCalc.normClean,reaRatesCalc.normCorr, by="eventID")
k600.ymin <- min(reaRates$k600.normClean, reaRates$k600.normCorr, na.rm=T)
k600.ymax <- max(reaRates$k600.normClean, reaRates$k600.normCorr, na.rm=T)
K600.ymin <- min(reaRates$K600.normClean, reaRates$K600.normCorr, na.rm=T)
K600.ymax <- max(reaRates$K600.normClean, reaRates$K600.normCorr, na.rm=T)

fig3a <- plot_ly(data = reaRatesCalc.normClean,
                 x = ~meanQ_lps,
                 y = ~k600.normClean,
                 type = 'scatter',
                 mode = 'markers',
                 #symbol = ~OK,
                 #symbols = c('x','circle'),
                 #color = ~OK,
                 #colors = c("red", "blue"),
                 marker = list(color = "blue", symbol = "circle"),
                 text = ~eventID) %>%
  layout(
    title = paste0(siteID," Manual Analysis - normalized, no salt corr"), 
    yaxis = list(title="k600 (m d<sup>-1</sup>)", 
                 color='black',
                 range=c(k600.ymin-5, k600.ymax+5)),
    xaxis = list(title="Q (lps)"), 
    font=list(size=18),
    margin=list(t=75, r =75),
    legend=list(x=100,y=1.4),
    showlegend = FALSE
  )
fig3b <- plot_ly(data = reaRatesCalc.normClean,
                 x = ~meanQ_lps,
                 y = ~K600.normClean,
                 type = 'scatter',
                 mode = 'markers',
                 #symbol = ~OK,
                 #symbols = c('x','circle'),
                 #color = ~OK,
                 #colors = c("red", "blue"),
                 marker = list(color = "blue", symbol = "circle"),
                 text = ~eventID) %>%
  layout(
    title = paste0(siteID," Manual Analysis - normalized, no salt corr"), 
    yaxis = list(title="K600 (d<sup>-1</sup>)", 
                 color='black',
                 range=c(K600.ymin-5, K600.ymax+5)),
    xaxis = list(title="Q (lps)"), 
    font=list(size=12),
    margin=list(t=75, r =75),
    legend=list(x=100,y=1.4),
    showlegend = FALSE
  )


fig4a <- plot_ly(data = reaRatesCalc.normCorr,
                 x = ~meanQ_lps,
                 y = ~k600.normCorr,
                 type = 'scatter', 
                 mode = 'markers',
                 symbol = ~OK, 
                 symbols = c('x','circle'),
                 color = ~OK,
                 colors = c("red", "blue"),
                 text = ~eventID) %>%
  layout(
    title = paste0(siteID," Manual Analysis - normalized, salt corr"), 
    yaxis = list(title="k600 (m d<sup>-1</sup>)", 
                 color='black',
                 range=c(k600.ymin-5, k600.ymax+5)),
    xaxis = list(title="Q (lps)"), 
    font=list(size=18),
    margin=list(t=75, r =75),
    legend=list(x=100,y=1.4),
    showlegend = FALSE
  )
fig4b <- plot_ly(data = reaRatesCalc.normCorr,
                 x = ~meanQ_lps,
                 y = ~K600.normCorr,
                 type = 'scatter', 
                 mode = 'markers',
                 symbol = ~OK, 
                 symbols = c('x','circle'),
                 color = ~OK,
                 colors = c("red", "blue"),
                 text = ~eventID) %>%
  layout(
    title = paste0(siteID," Manual Analysis - normalized, salt corr"), 
    yaxis = list(title="K600 (d<sup>-1</sup>)", 
                 color='black',
                 range=c(K600.ymin-5, K600.ymax+5)),
    xaxis = list(title="Q (lps)"), 
    font=list(size=12),
    margin=list(t=75, r =75),
    legend=list(x=100,y=1.4),
    showlegend = FALSE
  )



fig.K600 <- subplot(fig3a, fig4a, fig3b, fig4b, nrows=2,
                    titleY = TRUE, titleX = TRUE, margin = 0.1) %>% 
  layout(title = paste0(siteID,' k600 and K600'))

annotations = list( 
  list( 
    x = 0.2,  
    y = 1.0,  
    text = "normalized, no salt corr",  
    xref = "paper",  
    yref = "paper",  
    xanchor = "center",  
    yanchor = "bottom",  
    showarrow = FALSE 
  ),  
  list( 
    x = 0.8,  
    y = 1,  
    text = "normalized, salt corr",  
    xref = "paper",  
    yref = "paper",  
    xanchor = "center",  
    yanchor = "bottom",  
    showarrow = FALSE 
  ),  
  list( 
    x = 0.2,  
    y = 0.4,  
    text = "normalized, no salt corr",  
    xref = "paper",  
    yref = "paper",  
    xanchor = "center",  
    yanchor = "bottom",  
    showarrow = FALSE 
  ),
  list( 
    x = 0.8,  
    y = 0.4,  
    text = "normalized, salt corr",  
    xref = "paper",  
    yref = "paper",  
    xanchor = "center",  
    yanchor = "bottom",  
    showarrow = FALSE 
  ))

fig.K600 <- fig.K600 %>%layout(annotations = annotations)
fig.K600 #Save as Web Page to retain hovering labels

htmltools::save_html(
  html = htmltools::as.tags(
    x = plotly::toWebGL(fig.K600),
    standalone = TRUE),
  file = paste0(localPath,"plotsForTWG/",siteID,"_K600.html"))

###### Q plots #####

measuredDischargeDPID <- "DP4.00133.001"
fieldQ <- neonUtilities::loadByProduct(dpID = measuredDischargeDPID,site = siteID, check.size = FALSE)

figQ <- plot_ly(y = ~reaRatesCalc.normClean$meanQ_lps[reaRatesCalc.normClean$OK == "yes"],
                type = "box",
                boxpoints = "all",
                name = "rea exp - no salt corr")%>%
  add_trace(y = ~reaRatesCalc.normCorr$meanQ_lps[reaRatesCalc.normCorr$OK == "yes"], 
            name = "rea exp - salt corr") %>%
  add_trace(y = ~fieldQ$sdrc_gaugeDischargeMeas$streamDischarge, 
            name = "field Qs") %>%
  layout(
    title = paste0(siteID," Q comparison"), 
    yaxis = list(title="Q (lps)", 
                 color='black',
                 type = "log"),
    #xaxis = list(title="Q (lps)"), 
    font=list(size=18),
    margin=list(t=75, r =75),
    legend=list(x=100,y=1.4),
    showlegend = FALSE
  )
figQ #Save as Web Page to retain hovering labels

htmltools::save_html(
  html = htmltools::as.tags(
    x = plotly::toWebGL(figQ),
    standalone = TRUE),
  file = paste0(localPath,"plotsForTWG/",siteID,"_discharge.html"))

##############################################################################################
########## Prep and format things for stan and manual method ##########
##############################################################################################
#Starting from here for re-loading data to use new stan model as of 2/26/22
siteID <- "MART"
localPath <- "SET TO YOUR PATH/"
plotPathLossRates <- paste0(localPath,siteID,"/saltAndGasLossPlots")
plotPathBackground <- paste0(localPath,siteID,"/backgroundPlots")
plotPathTravelTime <- paste0(localPath,siteID,"/travelTimePlots")
stanPlots <- paste0(localPath,siteID)
load(paste0(localPath,siteID,"/MART_workspace.RData"))
localPath <- "SET TO YOUR PATH/"
plotPathLossRates <- paste0(localPath,siteID,"/saltAndGasLossPlots")
plotPathBackground <- paste0(localPath,siteID,"/backgroundPlots")
plotPathTravelTime <- paste0(localPath,siteID,"/travelTimePlots")
stanPlots <- paste0(localPath,siteID)

plotsOut <- reaRate::gas.loss.rate.plot(inputFile = reaFormatted,
                                        injectionTypeName = "injectionType",
                                        eventID = "eventID",
                                        stationToInjectionDistance = "stationToInjectionDistance",
                                        plateauGasConc = "plateauGasConc",
                                        corrPlatSaltConc = "corrPlatSaltConc",
                                        savePlotPath = plotPathLossRates)

#Combine travel time info and new slopes
reaRateManual <- reaRatesTrvlTime$outputDF
for(currEvt in reaRateManual$eventID){
  try(reaRateManual$slopeNormClean[reaRateManual$eventID == currEvt] <- unique(plotsOut$slopeNormClean[plotsOut$eventID == currEvt]))
  try(reaRateManual$slopeNormCorr[reaRateManual$eventID == currEvt] <- unique(plotsOut$slopeNormCorr[plotsOut$eventID == currEvt]))
}

reaRatesCalc <- reaRate::def.calc.reaeration(inputFile = reaRateManual,
                                             lossRateSF6 = "slopeNormClean",
                                             peakMaxVelocity = "peakMaxVelocity",
                                             meanDepth = "meanDepth",
                                             meanTemp = "meanTemp",
                                             outputSuffix = "clean")

reaRatesCalc <- reaRate::def.calc.reaeration(inputFile = reaRatesCalc,
                                             lossRateSF6 = "slopeNormCorr",
                                             peakMaxVelocity = "peakMaxVelocity",
                                             meanDepth = "meanDepth",
                                             meanTemp = "meanTemp",
                                             outputSuffix = "corr")

##############################################################################################
########## Calculate using cleaned gas data only (no salt or background correction) ##########
##############################################################################################

#Remove things that are messed up
plotsOutClean <- plotsOut[is.na(plotsOut$highDwnStrConcFlgClean),]
plotsOutClean <- plotsOutClean[!(!is.na(plotsOutClean$lowStOneFlagClean) & grepl("station\\.01",plotsOutClean$namedLocation)),]

#Create input file for reaeration calculations
KdfNamesClean <- c('eventID','stationID','fieldDischarge','originalDist','stationOneDist',
                   'meanRelativeCleanSF6','velocity','eventNum', 'temp',
                   'manualCleanBigK','manualCorrBigK', 'manualCleanLittlek','manualCorrLittlek')

KdfClean <- data.frame(matrix(data=NA, ncol=length(KdfNamesClean), nrow=nrow(plotsOutClean)))
names(KdfClean) <- KdfNamesClean

KdfClean$eventID <- plotsOutClean$eventID
KdfClean$stationID <- plotsOutClean$namedLocation
KdfClean$fieldDischarge <- plotsOutClean$fieldDischarge_lps * 60 / 1000 # Convert discharge from L/s to m3/min
KdfClean$meanRelativeCleanSF6 <- log(as.numeric(plotsOutClean$normGasClean))
KdfClean$originalDist <- plotsOutClean$stationToInjectionDistance
KdfClean$temp <- plotsOutClean$waterTemp

#Remove stations that don't have SF6 data.
eventIDToExclude <- KdfClean$eventID[is.na(KdfClean$meanRelativeCleanSF6)]
KdfClean <- KdfClean[!KdfClean$eventID %in% eventIDToExclude,]

# for(i in seq(along = KdfClean$relativeCleanSF6)){
#   KdfClean$meanRelativeCleanSF6[i] <-  mean(as.numeric(strsplit(KdfClean$relativeCleanSF6[i],"\\|")[[1]]), na.rm = TRUE)
#   KdfClean$meanRelativeCorrSF6[i] <-  mean(as.numeric(strsplit(KdfClean$relativeCorrSF6[i],"\\|")[[1]]), na.rm = TRUE)
# }

#Got to loop through to re-calculate downstream distances, relative SF6, and velocity
for(id in unique(KdfClean$eventID)){
  #Calculate distance relative to the first station
  minDist <- min(KdfClean$originalDist[KdfClean$eventID == id])
  KdfClean$stationOneDist[KdfClean$eventID == id] <- KdfClean$originalDist[KdfClean$eventID == id] - minDist
  
  #Populate velocity
  maxDist <- max(KdfClean$stationOneDist[KdfClean$eventID == id])
  KdfClean$velocity[KdfClean$eventID == id] <- reaRatesTrvlTime$outputDF$peakMaxVelocity[reaRatesTrvlTime$outputDF$eventID == id]
  
  #Populate the manual calculations as well
  KdfClean$manualCleanBigK[KdfClean$eventID == id] <- reaRatesCalc$K600.clean[reaRatesCalc$eventID == id]
  KdfClean$manualCorrBigK[KdfClean$eventID == id] <- reaRatesCalc$K600.corr[reaRatesCalc$eventID == id]
  KdfClean$manualCleanLittlek[KdfClean$eventID == id] <- reaRatesCalc$k600.clean[reaRatesCalc$eventID == id]
  KdfClean$manualCorrLittlek[KdfClean$eventID == id] <- reaRatesCalc$k600.corr[reaRatesCalc$eventID == id]
}
KdfClean$velocity <- KdfClean$velocity * 60*60*24 # Convert velocity from m/s to m/day

# Remove data missing travel time or discharge
cat("removing events with NA discharge: ", paste(unique(KdfClean$eventID[is.na(KdfClean$fieldDischarge)]), collapse = " & "), "\n")
cat("removing events with NA velocity: ", paste(unique(KdfClean$eventID[is.na(KdfClean$velocity)]), collapse = " & "), "\n")
cat("removing events with negative velocity: ", paste(unique(KdfClean$eventID[KdfClean$velocity <= 0]), collapse = " & "), "\n")
KdfClean <- KdfClean[!is.na(KdfClean$fieldDischarge),]
KdfClean <- KdfClean[!is.na(KdfClean$velocity),]
KdfClean <- KdfClean[KdfClean$velocity > 0,]

#Remove data with NA for temp
KdfClean <- KdfClean[!is.na(KdfClean$temp),]

Q <- rep(NA, length(unique(KdfClean$eventID)))
V <- rep(NA, length(unique(KdfClean$eventID)))
temp <- rep(NA,length(unique(KdfClean$eventID)))
for(num in seq(along = unique(KdfClean$eventID))){
  currEventID <- unique(KdfClean$eventID)[num]
  KdfClean$eventNum[KdfClean$eventID == currEventID] <- num
  Q[num] <- unique(KdfClean$fieldDischarge[KdfClean$eventID == currEventID])
  V[num] <- unique(KdfClean$velocity[KdfClean$eventID == currEventID])
  temp[num] <- mean(KdfClean$temp[KdfClean$eventID == currEventID], na.rm = TRUE)
}



#now for the Stan model on SF6 data only
sf6.datClean <- list( N= length(KdfClean$eventID), 
                      nexpt=length(Q), 
                      exptID=KdfClean$eventNum, 
                      dist=KdfClean$stationOneDist, 
                      logSf6=KdfClean$meanRelativeCleanSF6,
                      temp=temp,
                      Q=Q, 
                      V=V)

fitClean <- rstan::stan(paste0(localPath,"/NEON_multi_sf6.stan"), data = sf6.datClean,  iter = 5000, chains = 4)

print(fitClean,  digits_summary = 3)

rstan::traceplot(fitClean)

samples_fitClean<-rstan::extract(fitClean)

logk_sumClean <- rstan::summary(fitClean, pars = c("logk"), probs = c(0.5))$summary
Kd_sumClean <- rstan::summary(fitClean, pars = c("Kd"), probs = c(0.5))$summary
a_b_sumClean <- rstan::summary(fitClean, pars = c("a", "b"), probs = c(0.5))$summary
sigmaClean <- rstan::summary(fitClean, pars = c("sigma"), probs = c(0.5))$summary
sigma_exptClean <- rstan::summary(fitClean, pars = c("sigma_expt"), probs = c(0.5))$summary

pred_fullClean <- data.frame(eventNum = KdfClean$eventNum, dist=KdfClean$stationOneDist)
pred_fullClean$sf6Pred <- -Kd_sumClean[pred_fullClean$eventNum,1]*pred_fullClean$dist

eventNum.labs <- unique(KdfClean$eventID)
names(eventNum.labs) <- unique(KdfClean$eventNum)

ggplot2::ggplot(data= KdfClean, ggplot2::aes(y=meanRelativeCleanSF6, x=stationOneDist) ) + ggplot2::geom_point() +
  ggplot2::geom_line(data = pred_fullClean, ggplot2::aes(y=sf6Pred, x=dist)) +
  ggplot2::facet_wrap(~eventNum, labeller = ggplot2::labeller(eventNum = eventNum.labs))

dev.copy(jpeg,paste0(stanPlots,"/",siteID,'_clean_modelFits.jpg'), width = 16, height = 9, unit = "in", res = 1080)
dev.off()

# Generate plots of stan model and manual model outputs
lowY <- min(c(log(KdfClean$manualCleanBigK[is.finite(log(KdfClean$manualCleanBigK))]),logk_sumClean[,1]), na.rm = TRUE)
highY <- max(c(log(KdfClean$manualCleanBigK[is.finite(log(KdfClean$manualCleanBigK))]),logk_sumClean[,1]), na.rm = TRUE)

dev.new(width=16, height=9, unit="in", noRStudioGD = TRUE)
par(mfrow = c(1, 2))

plot (log(Q),
      logk_sumClean[,1],
      ylim = c(lowY*0.9, highY*1.1),
      xlab = "ln(discharge)",
      ylab = "ln(K600.clean)",
      main = paste0(siteID," Stan Model (clean)"))

for (i in 1:1000){
  lines(log(Q), samples_fitClean$a[i] + samples_fitClean$b[i]*log(Q), col='light gray')
}
points (log(Q), logk_sumClean[,1], col='blue', pch=16)
lines(log(Q), a_b_sumClean[1,1] +a_b_sumClean[2,1]*log(Q) )

plot(log(KdfClean$fieldDischarge),
     log(KdfClean$manualCleanBigK),
     ylim = c(lowY*0.9, highY*1.1),
     xlab = "ln(discharge)",
     ylab = "ln(K600.clean)",
     main = paste0(siteID," Manual Analysis (clean)"))

dev.copy(jpeg,paste0(stanPlots,"/",siteID,'_Outputs_clean.jpg'), width = 16, height = 9, unit = "in", res = 1080)
dev.off()

##################################################################################################
#now for the Stan model on Salt corrected SF6 data
##################################################################################################

#Remove things that are messed up
plotsOutCorr <- plotsOut[is.na(plotsOut$highDwnStrConcFlgCorr),]
plotsOutCorr <- plotsOutCorr[!(!is.na(plotsOutCorr$lowStOneFlagCorr) & grepl("station\\.01",plotsOutCorr$namedLocation)),]

#Create input file for reaeration calculations
KdfNamesCorr <- c('eventID','stationID','fieldDischarge','originalDist','stationOneDist',
                  'meanRelativeCorrSF6','velocity','eventNum', 'temp',
                  'manualCleanBigK','manualCorrBigK', 'manualCleanLittlek','manualCorrLittlek')

KdfCorr <- data.frame(matrix(data=NA, ncol=length(KdfNamesCorr), nrow=nrow(plotsOutCorr)))
names(KdfCorr) <- KdfNamesCorr

KdfCorr$eventID <- plotsOutCorr$eventID
KdfCorr$stationID <- plotsOutCorr$namedLocation
KdfCorr$fieldDischarge <- plotsOutCorr$fieldDischarge_lps * 60 / 1000 # Convert discharge from L/s to m3/min
KdfCorr$meanRelativeCorrSF6 <- log(as.numeric(plotsOutCorr$normGasCorr))
KdfCorr$originalDist <- plotsOutCorr$stationToInjectionDistance
KdfCorr$temp <- plotsOutCorr$waterTemp

#Remove stations that don't have SF6 data.
eventIDToExclude <- KdfCorr$eventID[is.na(KdfCorr$meanRelativeCorrSF6)]
KdfCorr <- KdfCorr[!KdfCorr$eventID %in% eventIDToExclude,]

#Got to loop through to re-calculate downstream distances, relative SF6, and velocity
for(id in unique(KdfCorr$eventID)){
  #Calculate distance relative to the first station
  minDist <- min(KdfCorr$originalDist[KdfCorr$eventID == id])
  KdfCorr$stationOneDist[KdfCorr$eventID == id] <- KdfCorr$originalDist[KdfCorr$eventID == id] - minDist
  
  #Populate velocity
  maxDist <- max(KdfCorr$stationOneDist[KdfCorr$eventID == id])
  KdfCorr$velocity[KdfCorr$eventID == id] <- reaRatesTrvlTime$outputDF$peakMaxVelocity[reaRatesTrvlTime$outputDF$eventID == id]
  
  #Populate the manual calculations as well
  KdfCorr$manualCleanBigK[KdfCorr$eventID == id] <- reaRatesCalc$K600.clean[reaRatesCalc$eventID == id]
  KdfCorr$manualCorrBigK[KdfCorr$eventID == id] <- reaRatesCalc$K600.corr[reaRatesCalc$eventID == id]
  KdfCorr$manualCleanLittlek[KdfCorr$eventID == id] <- reaRatesCalc$k600.clean[reaRatesCalc$eventID == id]
  KdfCorr$manualCorrLittlek[KdfCorr$eventID == id] <- reaRatesCalc$k600.corr[reaRatesCalc$eventID == id]
}
KdfCorr$velocity <- KdfCorr$velocity * 60*60*24 # Convert velocity from m/s to m/day

# Remove data missing travel time or discharge
cat("removing events with NA discharge: ", paste(unique(KdfCorr$eventID[is.na(KdfCorr$fieldDischarge)]), collapse = " & "), "\n")
cat("removing events with NA velocity: ", paste(unique(KdfCorr$eventID[is.na(KdfCorr$velocity)]), collapse = " & "), "\n")
cat("removing events with negative velocity: ", paste(unique(KdfCorr$eventID[KdfCorr$velocity <= 0]), collapse = " & "), "\n")
KdfCorr <- KdfCorr[!is.na(KdfCorr$fieldDischarge),]
KdfCorr <- KdfCorr[!is.na(KdfCorr$velocity),]
KdfCorr <- KdfCorr[KdfCorr$velocity > 0,]

#Remove data with NA for temp
KdfCorr <- KdfCorr[!is.na(KdfCorr$temp),]

Q <- rep(NA, length(unique(KdfCorr$eventID)))
V <- rep(NA, length(unique(KdfCorr$eventID)))
temp <- rep(NA,length(unique(KdfCorr$eventID)))
for(num in seq(along = unique(KdfCorr$eventID))){
  currEventID <- unique(KdfCorr$eventID)[num]
  KdfCorr$eventNum[KdfCorr$eventID == currEventID] <- num
  Q[num] <- unique(KdfCorr$fieldDischarge[KdfCorr$eventID == currEventID])
  V[num] <- unique(KdfCorr$velocity[KdfCorr$eventID == currEventID])
  temp[num] <- mean(KdfCorr$temp[KdfCorr$eventID == currEventID], na.rm = TRUE)
}

sf6.datCorr<- list( N= length(KdfCorr$eventID),
                    nexpt=length(Q),
                    exptID=KdfCorr$eventNum,
                    dist=KdfCorr$stationOneDist,
                    logSf6=KdfCorr$meanRelativeCorrSF6,
                    temp=temp,
                    Q=Q,
                    V=V)

fitCorr <- rstan::stan(paste0(localPath,"/NEON_multi_sf6.stan"), data = sf6.datCorr,  iter = 5000, chains = 4)

print(fitCorr,  digits_summary = 3)

rstan::traceplot(fitCorr)

samples_fitCorr<-rstan::extract(fitCorr)

logk_sumCorr <- rstan::summary(fitCorr, pars = c("logk"), probs = c(0.5))$summary
Kd_sumCorr <- rstan::summary(fitCorr, pars = c("Kd"), probs = c(0.5))$summary
a_b_sumCorr <- rstan::summary(fitCorr, pars = c("a", "b"), probs = c(0.5))$summary

pred_fullCorr <- data.frame(eventNum = KdfCorr$eventNum, dist=KdfCorr$stationOneDist)
pred_fullCorr$sf6Pred <- -Kd_sumCorr[pred_fullCorr$eventNum,1]*pred_fullCorr$dist

eventNumCorr.labs <- unique(KdfCorr$eventID)
names(eventNumCorr.labs) <- unique(KdfCorr$eventNum)

ggplot2::ggplot(data= KdfCorr, ggplot2::aes(y=meanRelativeCorrSF6, x=stationOneDist) ) + ggplot2::geom_point() +
  ggplot2::geom_line(data = pred_fullCorr, ggplot2::aes(y=sf6Pred, x=dist)) +
  ggplot2::facet_wrap(~eventNum, labeller = ggplot2::labeller(eventNum = eventNumCorr.labs))

dev.copy(jpeg,paste0(stanPlots,"/",siteID,'_corr_modelFits.jpg'), width = 16, height = 9, unit = "in", res = 1080)
dev.off()

# Generate plots of stan model and manual model outputs
lowY <- min(c(log(KdfCorr$manualCorrBigK[is.finite(log(KdfCorr$manualCorrBigK))]),log(logk_sumCorr[,1])), na.rm = TRUE)
highY <- max(c(log(KdfCorr$manualCorrBigK[is.finite(log(KdfCorr$manualCorrBigK))]),log(logk_sumCorr[,1])), na.rm = TRUE)

dev.new(width=16, height=9, unit="in", noRStudioGD = TRUE)
par(mfrow = c(1, 2))

plot (log(Q),
      logk_sumCorr[,1],
      ylim = c(lowY*0.9, highY*1.1),
      xlab = "ln(discharge)",
      ylab = "ln(K600.corr)",
      main = paste0(siteID," Stan Model (corr)"))

for (i in 1:1000){
  lines(log(Q), samples_fitCorr$a[i] + samples_fitCorr$b[i]*log(Q), col='light gray')
}
points (log(Q), logk_sumCorr[,1], col='blue', pch=16)
lines(log(Q), a_b_sumCorr[1,1] +a_b_sumCorr[2,1]*log(Q) )

plot(log(KdfCorr$fieldDischarge),
     log(KdfCorr$manualCorrBigK),
     ylim = c(lowY*0.9, highY*1.1),
     xlab = "ln(discharge)",
     ylab = "ln(K600.corr)",
     main = paste0(siteID," Manual Analysis (corr)"))

dev.copy(jpeg,paste0(stanPlots,"/",siteID,'_Outputs_corr.jpg'), width = 16, height = 9, unit = "in", res = 1080)
dev.off()


save.image(paste0(stanPlots,"/",siteID, "_withStan.RData"))

#JASM Plots
siteID <- "MART"
localPath <- "SET TO YOUR PATH/"
plotPathLossRates <- paste0(localPath,siteID,"/saltAndGasLossPlots")
plotPathBackground <- paste0(localPath,siteID,"/backgroundPlots")
plotPathTravelTime <- paste0(localPath,siteID,"/travelTimePlots")
stanPlots <- paste0(localPath,siteID)
load(paste0(stanPlots,"/",siteID, "_withStan.RData"))

plot (log(Q),
      logk_sumCorr[,1],
      xlab = "ln(discharge, lps)",
      ylab = "ln(K600, d-1)",
      main = paste0(siteID," Stan Model"))

for (i in 1:1000){
  lines(log(Q), samples_fitCorr$a[i] + samples_fitCorr$b[i]*log(Q), col='light gray')
}
points (log(Q), logk_sumCorr[,1], col='blue', pch=16)
lines(log(Q), a_b_sumCorr[1,1] +a_b_sumCorr[2,1]*log(Q) )

min(exp(logk_sumCorr[,1]))
max(exp(logk_sumCorr[,1]))
