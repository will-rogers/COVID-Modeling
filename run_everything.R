### Running files 
setwd("~/Downloads/Research/COVID-19/LAMP/COVID-Modeling")

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

set.seed(12345)
output <- uni_sims_par(tests=c(1000), compliance=.8, introductions = 20, 
                       ppn_sympt = 0.4, care.seeking = .5, R0.on = 3, R0.off = 3, 
                       test.scenario = c("No Delay"), sens.pcr = c(.99), spec.pcr = .99,
                       sens.lamp = c(1), spec.lamp = .99, 
                       ncores=NULL)

out <- running_tots(output)
ggplot(out, aes(x = day, y = cum_all.symptomatics.on+cum_all.symptomatics.off, color = factor(sens.pcr))) +
  geom_smooth() + 
  facet_grid(tests~.)+
  theme_classic()





