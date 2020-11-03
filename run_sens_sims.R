### Running files 
setwd("~/Downloads/Research/COVID-19/LAMP/COVID-Modeling")

# Always rerun these if you edit the source file
source("sir_sens_function.R")
source("sens_sim_function.R")
source("sens_sims_par_function.R")
source("plotting_functions.R")

# These are good to have on hand
library(ggplot2)
library(dplyr)
library(data.table)
library(mc2d)
library(abind)

set.seed(12345)

output <- sens_sims_par(tests = c(0,500,1000), compliance = 1, introduction = 25,
                        ppn_sympt = 0.2, care.seeking = 1,
                        R0 = 2.5, sensitivity = .5, times = 1, 
                        form = c("uniform-max","uniform-mean","uniform-min",
                                 "normal", "early beta", "late beta"), ncores = 2)

out <- output %>% 
  group_by(group) %>% 
  mutate(cum.cases.on = cumsum(new.cases))


ggplot(out, aes(x = day, y = cum.cases.on, color = form)) +
  geom_smooth() +
  facet_grid(.~tests)


output$active.inf

out <- running_tots(output)
head(out)

ggplot(out, aes(x = day, y = cum_cases.on+cum_cases.off, color = factor(sensitivity))) +
  geom_smooth() +
  facet_grid(tests~testing.per.patient)+
  theme_classic()