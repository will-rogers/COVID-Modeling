### Running files 
setwd("~/Downloads/Research/COVID-19/COVIDstuff-master/LAMP-Project")

# Always rerun these if you edit the source file
source("sir_lamp_function.R")
source("sir_simple_step_function.R")
source("uni_sim_function.R")
source("uni_sim_par_function.R")
source("plotting_functions.R")

# These are good to have on hand
library(ggplot2)
library(dplyr)
library(data.table)
library(mc2d)
library(abind)

output <- uni_sims_par(tests=c(1,500,1000), compliance=1, introductions = 25, 
                       ppn_sympt = 0.8, care.seeking = 1, R0.on = 2.5, R0.off = 2.5, 
                       test.scenario = c("No Delay"),
                       sensitivity = c(.8, .9, 1), specificity = 1, times = c(1),
                       ncores=2)

out <- running_tots(output)
head(out)
