##############################################################################################
#' @title Constant rate or slug travel time between two stations

#' @author 
#' Kaelin M. Cawley \email{kcawley@battelleecology.org} \cr

#' @description This function determines the timestamp of the peak tracer conductance for use 
#' calculating travel time between an upstream and downstream sensor set.

#' @importFrom neonUtilities stackByTable

#' @param dataDir User input of the directory holding the logger files [string]
#' @param currEventID User input of the eventID of the tracer experiment [string]
#' @param injectionType User input of the injection type either "constant" or "slug" [string]
#' @param bPlot User input to display a visual check, default is FALSE [boolean]

#' @return This function returns the peak tracer timestamp [dateTime]

#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

#' @keywords surface water, streams, rivers, velocity, travel time, reaeration, metabolism

#' @examples
#' #Using an example file
#' travelTime <- def.calc.travelTime(
#' dataDir = paste(path.package("reaRate"),"inst\\extdata", sep = "\\"), 
#' currEventID = "GUIL.20150129", injectionType = "constant", bPlot = T)

#' @seealso def.calc.reaeration.R for calculating reaeration rates

#' @export

# changelog and author contributions / copyrights
#   Kaelin M. Cawley (2017-08-03)
#     original creation
##############################################################################################
def.calc.travelTime <- function(
  dataDir,
  currEventID,
  injectionType,
  bPlot = FALSE
  ){
  
  #Stack field and external lab data if needed
  if(!dir.exists(substr(dataDir, 1, (nchar(dataDir)-4)))){
    neonUtilities::stackByTable(dpID="DP1.20190.001",filepath=dataDir)
  }
  
  #Read in stacked logger data
  #Allows for using the reaeration tables in addition to the salt-based discharge tables
  allFolders <- list.files(paste(gsub("\\.zip","",dataDir), sep = "/"))
  eventFolder <- allFolders[grepl(paste0(substr(currEventID,1,4), ".*", substr(currEventID,6,9),"\\-", substr(currEventID,10,11)), allFolders)]
  allFiles <- list.files(paste(gsub("\\.zip","",dataDir), eventFolder, sep = "/"))
  
  if(length(allFiles[grepl("conductivityFieldData", allFiles)]) == 0){
    
    print(paste0("Conductivity logger data not available for ", currEventID))
    
  }else{
    
    loggerFile <- allFiles[grepl("conductivityFieldData", allFiles)]
    loggerData <- read.csv(
      paste(gsub("\\.zip","",dataDir), eventFolder, loggerFile, sep = "/"), 
      stringsAsFactors = F)
  
    s1LoggerData <- loggerData[loggerData$hoboSampleID == paste0(substr(currEventID, 1, 4), "_S1_", substr(currEventID, 6, 13)),]
    s1LoggerData <- s1LoggerData[order(s1LoggerData$measurementNumber),]
    
    s4LoggerData <- loggerData[loggerData$hoboSampleID == paste0(substr(currEventID, 1, 4), "_S4_", substr(currEventID, 6, 13)),]
    s4LoggerData <- s4LoggerData[order(s4LoggerData$measurementNumber),]
    
    if(length(s1LoggerData[[1]]) <= 0){
      print(paste0("Conductivity logger data not available for ", currEventID, ", station S1"))
    }else if(length(s4LoggerData[[1]]) <= 0){
      print(paste0("Conductivity logger data not available for ", currEventID, ", station S4"))
    }else{
      #If low range isn't collected use the full range
      ifelse(all(is.na(s1LoggerData$lowRangeSpCondNonlinear)),
             condDataS1 <- s1LoggerData$fullRangeSpCondNonlinear,
             condDataS1 <- s1LoggerData$lowRangeSpCondNonlinear)
    
      ifelse(all(is.na(s4LoggerData$lowRangeSpCondNonlinear)),
             condDataS4 <- s4LoggerData$fullRangeSpCondNonlinear,
             condDataS4 <- s4LoggerData$lowRangeSpCondNonlinear)
  
      #Slug injections, find the index/measurement number of the tracer peak
      if(injectionType == "slug"){
      #After checking for spikes, return peak value
        s1peakLoc <- slugTime(condDataS1)
        s1PkTime <- strptime(s1LoggerData$dateTimeLogger[s1peakLoc], tz = "UTC", format = "%Y-%m-%dT%H:%M:%S")
      
        s4peakLoc <- slugTime(condDataS4)
        s4PkTime <- strptime(s4LoggerData$dateTimeLogger[s4peakLoc], tz = "UTC", format = "%Y-%m-%dT%H:%M:%S")
      
        travelTime <- difftime(s4PkTime, s1PkTime, units = "secs")
      
      }else if(injectionType == "constant"){
        #Determine the data frequency, usually every 10 seconds for NEON data
        #Look a few readings in since sometimes there are weird timestamps when it's first turned on
        #freq_S1 <- difftime(strptime(s1LoggerData$dateTimeLogger[6], tz = "UTC", format = "%Y-%m-%dT%H:%M:%S"),strptime(s1LoggerData$dateTimeLogger[5], tz = "UTC", format = "%Y-%m-%dT%H:%M:%S"), units = "secs")
        #freq_S4 <- difftime(strptime(s4LoggerData$dateTimeLogger[6], tz = "UTC", format = "%Y-%m-%dT%H:%M:%S"),strptime(s4LoggerData$dateTimeLogger[5], tz = "UTC", format = "%Y-%m-%dT%H:%M:%S"), units = "secs")
      
        s1peakLoc <- constTime(reaCond = condDataS1, 
                               plotName = paste0(currEventID, "_S1"),
                               bPlot = bPlot)
        if(s1peakLoc == -9999){
          print(paste0("Conductivity logger data did not contain a peak/plateau ", currEventID, ", station S1"))
          return(-9999)
        }else{
          s1PkTime <- strptime(s1LoggerData$dateTimeLogger[s1peakLoc], tz = "UTC", format = "%Y-%m-%dT%H:%M:%S")
      
          s4peakLoc <- constTime(reaCond = condDataS4, 
                                 plotName = paste0(currEventID, "_S4"),
                                 bPlot = bPlot)
          if(s4peakLoc == -9999){
            print(paste0("Conductivity logger data did not contain a peak/plateau ", currEventID, ", station S4"))
            return(-9999)
          }else{
            s4PkTime <- strptime(s4LoggerData$dateTimeLogger[s4peakLoc], tz = "UTC", format = "%Y-%m-%dT%H:%M:%S")
      
            travelTime <- difftime(s4PkTime, s1PkTime, units = "secs")
          }
        }
        
      }else{
        print(paste0("Injection Type was not specified for ", currEventID, "and travel time could not be determined."))
      }
      if(travelTime <= 0){
        print(paste0("Conductivity logger data travel time calculated less than 0 for ", currEventID))
        #return(-9999)
      }else{
        return(travelTime)
      }
    }
  }
}


                                           