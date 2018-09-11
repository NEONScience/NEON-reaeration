
library(devtools)
library(roxygen2)

setwd("C:/Users/kcawley/Documents/GitHub/how-to-make-a-data-product/REST_R")
install("restR")
library(restR)

setwd("C:/Users/kcawley/Documents/GitHub/how-to-make-a-data-product/REST_R/restR")
document()
devtools::check()
