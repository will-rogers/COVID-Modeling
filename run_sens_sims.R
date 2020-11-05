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

output <- sens_sims_par(tests = c(100,150,500), compliance = 1, introduction = 5,
                        ppn_sympt = 0.2, care.seeking = 1,
                        R0 = 2.5, sensitivity = 1, times = 1, 
                        form = c("uniform-max",
                                 "normal", "late beta"), ncores = 2)

out <- output %>% 
  group_by(group) %>% 
  mutate(cum.cases = cumsum(new.cases))

ggplot(out, aes(x = day, y = cum.cases, color = form)) +
  geom_smooth() +
  facet_grid(.~tests) + 
  labs(y = "Cumulative Cases") +
  theme_classic()


output$active.inf

out <- running_tots(output)
head(out)

ggplot(out, aes(x = day, y = cum_cases.on+cum_cases.off, color = factor(sensitivity))) +
  geom_smooth() +
  facet_grid(tests~testing.per.patient)+
  theme_classic()



hist(rpois(1000, 2))
1000*.995





