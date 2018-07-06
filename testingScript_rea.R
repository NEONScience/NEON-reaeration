
library(devtools)
library(roxygen2)

setwd("C:/Users/kcawley/Documents/GitHub/NEON-reaeration")
#setwd("C:/Users/Kaelin/Documents/GitHub/biogeochemistryIPT/reaeration/Science Only/rCodeForRelease")
#create("reaRate")
#install_github("NEONScience/NEON-utilities/neonUtilities", force = T)
install("reaRate")
library(reaRate)

dataDir <- "C:/Users/kcawley/Downloads/NEON_reaeration.zip"
#dataDir <- "C:/Users/Kaelin/Downloads/NEON_reaeration.zip"
#dataDir <- "C:/Users/Kcawley/Desktop/NEON_reaeration.zip"

#For use with the API functionality
dataDir <- "API"
site <- "KING"
#site <- "all"

reaFormatted <- def.format.reaeration(dataDir = dataDir, site = site, fieldQ = TRUE)
#write.csv(reaFormatted, "C:/Users/kcawley/Documents/GitHub/biogeochemistryIPT/reaeration/Science Only/rCodeForRelease/reaRate/inst/extdata/reaTestData.csv", row.names = F)
#write.csv(condDataS1, "C:/Users/kcawley/Documents/GitHub/biogeochemistryIPT/reaeration/Science Only/rCodeForRelease/reaRate/inst/extdata/condDataS1.csv", row.names = F)

#reaRatesCalc <- def.calc.reaeration(inputFile = reaFormatted, loggerFile = , dataDir = dataDir, plot = TRUE)
reaRatesCalc <- def.calc.reaeration(inputFile = reaFormatted, 
                                    dataDir = "C:/Users/kcawley/Documents/GitHub/NEON-reaeration/filesToStack20190/stackedFiles/", 
                                    loggerFile = "rea_conductivityFieldData.csv", 
                                    plot = TRUE,
                                    savePlotPath = "H:/Operations Optimization/savedPlots")

outputDF <- reaRatesCalc$outputDF
inputFile <- reaRatesCalc$inputFile

setwd("C:/Users/kcawley/Documents/GitHub/NEON-reaeration/reaRate")
#setwd("C:/Users/Kaelin/Documents/GitHub/biogeochemistryIPT/reaeration/Science Only/rCodeForRelease/reaRate")
document()
#devtools::use_data(reaFormatted, reaFormatted)
#devtools::use_data(condDataS1, condDataS1)
devtools::check()
