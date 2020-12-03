### Running files 
setwd("~/Downloads/Research/COVID-19/LAMP/COVID-Modeling")

# These are good to have on hand
library(ggplot2)
library(dplyr)
library(data.table)
library(mc2d)
library(abind)
library(pbapply)

# Always rerun these if you edit the source file
source("sir_lamp_function.R")
source("sir_simple_step_function.R")
source("uni_sim_function.R")
source("uni_sim_par_function.R")
source("plotting_functions.R")

set.seed(12345)
output <- uni_sims_par(tst = c(0), 
                       test.timeline = c("Sustained"),
                       compliance = c(0.5), # worst-best
                       init.prev = c(0.002), # best-mid-worst
                       ppn_sympt = c(0.2), # best-worst
                       care.seeking = c(0.3), # worst-best
                       R0.on = c(2.4),  # best-mid-worst
                       R0.off = c(2.4), # best-mid-worst
                       test.scenario = c("No Delay"),# best-worst
                       sens.pcr = 1, # for all intents this is reasonable
                       spec.pcr = 1, # for all intents this is reasonable
                       sens.lamp = c(.9), # bad day/lab to good day/lab
                       spec.lamp = .99, # not dynamic yet
                       lamp.diagnostic = c(F), # will affect PCR demand
                       community.intro.daily.on = 1, 
                       community.prob.daily.on = c(0.1), # every ten, 2, or 1 days
                       community.intro.daily.off = 1, 
                       community.prob.daily.off = c(0.1), # every ten, 2, or 1 days
                       immunity = c(0.0), # based on Fall 20 symptm, extrapolation to total
                       N0 = 16750, #campus pop
                       on.campus.prop = .25, #on/off division, only matters if we change on/off characteristics
                       contact.tracing.limit = c(0), # limit on number of contact traces per day
                       pooling = c(1), 
                       pooling.multi = c(1), #is the effect of pooling on accuracy 1:1, or do added pools only reduce sensitivity slightly?
                       days = 150, #simulation days
                       sims = 50, # number of simulations
                       ncores=NULL)

output$lamp.diagnostic.f <- factor(output$lamp.diagnostic, levels = c("FALSE", "TRUE"), labels = c("PCR Conf.", "LAMP Diag."))
out <- output %>% 
  group_by(group) %>% 
  mutate(cum.cases.on = cumsum(new.cases.on),
         cum.reporting.symptoms.on = cumsum(reporting.symptoms.on),
         cum.all.symptomatics.on = cumsum(all.symptomatics.on),
         cum.all.asymptomatics.on = cumsum(positive.asympt.on),
         cum.cases.off = cumsum(new.cases.off),
         cum.reporting.symptoms.off = cumsum(reporting.symptoms.off),
         cum.all.symptomatics.off = cumsum(all.symptomatics.off),
         cum.all.asymptomatics.off = cumsum(positive.asympt.off),
         cum.sum.missed = cumsum(missed.pcr),
         pcr.demand = cumsum(symp.pcr) + cumsum(asymp.pcr)
  )

# out %>% 
#   filter(day == max(day)) %>% 
#   group_by(test.timeline, sens.lamp, init.prev) %>% 
#   summarize(mean = mean(cum.cases.on+cum.cases.off),
#                      sd = sd(cum.cases.on+cum.cases.off),
#                      upper = mean + 1.96*sd,
#                      lower = mean - 1.96*sd) %>% 
#   ggplot(aes(x = test.timeline, y = mean, color = test.timeline)) +
#   geom_point(stat="identity") +
#   geom_errorbar(aes(ymin = lower, ymax = upper)) +
#   facet_grid(sens.lamp~init.prev) +
#   theme_classic()


# write.csv(out, "out.csv")
# 
# plot <- out %>% 
#   group_by(day, lamp.diagnostic.f) %>% 
#   summarize(mean = mean(symp.pcr+asymp.pcr),
#          sd = sd(symp.pcr+asymp.pcr),
#          upper = mean + 1.96*sd,
#          lower = mean - 1.96*sd) 
# ggplot(plot, aes(x = day, y = mean, color = lamp.diagnostic.f, fill = lamp.diagnostic.f)) +
#   geom_line() + 
#   geom_ribbon(aes(ymax = upper, ymin = lower), alpha = .5) +
#   theme_classic()
# 
# plot <- out %>% 
#   filter(day == max(day)) %>% 
#   group_by(lamp.diagnostic.f) %>% 
#   summarize(mean = mean(cum.sum.missed),
#             sd = sd(cum.sum.missed),
#             upper = mean + 1.96*sd,
#             lower = mean - 1.96*sd) 
# ggplot(plot, aes(x = lamp.diagnostic.f, fill = lamp.diagnostic.f)) +
#   geom_bar(aes(y = mean), stat = "identity") + 
#   geom_errorbar(aes(ymax = upper, ymin = lower)) +
#   theme_classic()
# 
df <- data.frame(day = seq(17, by = 7, length.out = 13),
                 new_MSU_cases = c(NA,3,7,66,43,60,65,99,132,212,265,203,105),
                 cumulative_MSU_cases = c(38,41,48,114,157,217,282,381,513,725,990,1193,1298))
out %>%
  filter(tests == 0) %>%
  group_by(day) %>%
  summarize(mean = mean(cum.reporting.symptoms.on+cum.reporting.symptoms.off)+0.0025*16750,
                     sd = sd(cum.reporting.symptoms.on+cum.reporting.symptoms.off),
                     upper = mean + 1.96*sd,
                     lower = mean - 1.96*sd) %>%
  ggplot(aes(x = day, y = mean)) +
  geom_line() +
  geom_ribbon(aes(ymax = upper, ymin = lower), alpha = 0.25) +
  geom_point(data = df, aes(x = day, y = cumulative_MSU_cases, color = "Empirical")) +
  labs(y="Estimated Cumulative Cases \n (95% Of Simulations)",
       x = "Day of Semester") +
  scale_color_discrete(name = "",
                       labels = "Weekly Cumulative \n Case Counts") +
  theme_classic()
# ggsave("ComparisontoReal.png")
# 
# out %>% 
#   group_by(tests, sens.lamp, day) %>%
#   summarize(mean = mean(cum.cases.on+cum.cases.off),
#                        sd = sd(cum.cases.on+cum.cases.off),
#                        upper = mean + 1.96*sd,
#                        lower = mean - 1.96*sd) %>%
#   ggplot(aes(x = day, y = mean, color = factor(sens.lamp), 
#              fill = factor(sens.lamp), group = factor(sens.lamp))) +
#   geom_line() +
#   geom_ribbon(aes(ymax = upper, ymin = lower), alpha = 0.25, colour = NA) +
#   facet_grid(.~tests) +
#   labs(y="Estimated Cumulative Cases \n (95% Of Simulations)",
#         x = "Day of Semester") +
#   scale_color_discrete(name = "LAMP \nSensitivity",
#                        labels = c("80%","90%","100%")) +
#   scale_fill_discrete(name = "LAMP \nSensitivity",
#                        labels = c("80%","90%","100%")) +
#   theme_classic()
# ggsave("LAMPSensvtime.png")
# 
# out %>% 
#   filter(day == max(day)) %>% 
#   group_by(tests, sens.lamp) %>%
#   summarize(mean = mean(cum.cases.on+cum.cases.off),
#             sd = sd(cum.cases.on+cum.cases.off),
#             upper = mean + 1.96*sd,
#             lower = mean - 1.96*sd) %>%
#   ggplot(aes(x = sens.lamp, y = mean, color = factor(sens.lamp), 
#              fill = factor(sens.lamp), group = factor(sens.lamp))) +
#   geom_bar(stat = "identity") +
#   geom_errorbar(aes(ymax = upper, ymin = lower),color="black") +
#   facet_grid(.~tests) +
#   labs(y="Estimated Cumulative Cases \n (95% Of Simulations)",
#        x = "Day of Semester") +
#   scale_color_discrete(name = "LAMP \nSensitivity",
#                        labels = c("80%","90%","100%")) +
#   scale_fill_discrete(name = "LAMP \nSensitivity",
#                       labels = c("80%","90%","100%")) +
#   theme_classic()
# ggsave("LAMPSens.png")
# 
# out %>% 
#   group_by(tests, sens.lamp, day) %>%
#   summarize(mean = mean((1/.05)*symp.pcr+asymp.pcr),
#             sd = sd((1/.05)*symp.pcr+asymp.pcr),
#             upper = mean + 1.96*sd,
#             lower = mean - 1.96*sd) %>%
#   ggplot(aes(x = day, y = mean, color = factor(sens.lamp), 
#              fill = factor(sens.lamp), group = factor(sens.lamp))) +
#   geom_line() +
#   geom_ribbon(aes(ymax = upper, ymin = lower), alpha = 0.25, colour = NA) +
#   facet_grid(.~tests) +
#   labs(y="Estimated Daily PCRs Required \n (95% Of Simulations, 10% Positivity)",
#        x = "Day of Semester") +
#   scale_color_discrete(name = "LAMP \nSensitivity",
#                        labels = c("80%","90%","100%")) +
#   scale_fill_discrete(name = "LAMP \nSensitivity",
#                       labels = c("80%","90%","100%")) +
#   theme_classic()
# ggsave("LAMP_aids_PCR_demand.png")
# 
# out %>% 
#   filter(tests > 0) %>% 
#   group_by(tests, sens.lamp, day) %>%
#   summarize(mean = mean(missed.pcr),
#             sd = sd(missed.pcr),
#             upper = mean + 1.96*sd,
#             lower = mean - 1.96*sd) %>%
#   ggplot(aes(x = day, y = mean, color = factor(sens.lamp), 
#              fill = factor(sens.lamp), group = factor(sens.lamp))) +
#   geom_line() +
#   geom_ribbon(aes(ymax = upper, ymin = lower), alpha = 0.25, colour = NA) +
#   facet_grid(.~tests) +
#   labs(y="Estimated Daily PCRs Required \n (95% Of Simulations, 10% Positivity)",
#        x = "Day of Semester") +
#   scale_color_discrete(name = "LAMP \nSensitivity",
#                        labels = c("80%","90%","100%")) +
#   scale_fill_discrete(name = "LAMP \nSensitivity",
#                       labels = c("80%","90%","100%")) +
#   theme_classic()
# 
# out %>% 
#   filter(day == max(day)) %>% 
#   group_by(tests, sens.lamp, day) %>%
#   summarize(mean = mean(cum.sum.missed),
#             sd = sd(cum.sum.missed),
#             upper = mean + 1.96*sd,
#             lower = mean - 1.96*sd) %>%
#   ggplot(aes(x = factor(sens.lamp), y = mean, color = factor(sens.lamp), 
#              fill = factor(sens.lamp), group = factor(sens.lamp))) +
#   geom_point() +
#   geom_errorbar(aes(ymax = upper, ymin = lower)) +
#   facet_grid(.~tests) +
#   labs(y="Cumulative Missed Cases \n (95% Of Simulations)",
#        x = "Sensitivity") +
#   scale_color_discrete(name = "LAMP \nSensitivity",
#                        labels = c("80%","85%","90%","95%")) +
#   scale_fill_discrete(name = "LAMP \nSensitivity",
#                       labels = c("80%","85%","90%","95%")) +
#   theme_classic() +
#   theme(axis.ticks.x = element_blank(),
#         axis.text.x = element_blank())
# ggsave("LAMP_missed_cases.png")


