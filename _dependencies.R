# Development and report packages that renv does not detect on its own.
# This file only exists so that renv records them in renv.lock; it is never run.
library(devtools)
library(usethis)
library(testthat)
library(roxygen2)
library(knitr)
library(rmarkdown)
library(quarto)
library(ggplot2)
# Suggests in DESCRIPTION (comparison methods and spatial data)
library(SpatialBSS)
library(sf)
