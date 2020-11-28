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
output <- uni_sims_par(tst = 500, 
                       compliance = 0.75, 
                       introductions = 5, 
                       ppn_sympt = .8, 
                       care.seeking = 0.5, 
                       R0.on = 3,  R0.off = 1.5, 
                       test.scenario = c("No Delay"),
                       sens.pcr = .99, spec.pcr = .99, 
                       sens.lamp = c(1), spec.lamp = .99, 
                       lamp.diagnostic = F, 
                       size.intro.on = 1, prob.into.on =0.1,
                       size.intro.off = 1, prob.into.off =0.1,
                       immunity = 0.1, 
                       N0 = 16750, 
                       on.campus.prop = .25, 
                       contact.tracing.limit = 100,
                       pooling = 4, pooling.multi = 1,
                       days = 300, sims = 5,
                       engage.lamp = 25,
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

write.csv(out, "comparison.csv")

set.seed(12345)
output <- uni_sims_par(tests=c(0,100,250,500,1000,2000,3000,4000), compliance=1, introductions = 50, 
                       ppn_sympt = 0.2, care.seeking = 1, R0.on = 2.6, R0.off = 2.6, 
                       test.scenario = c("1 Day"), sens.pcr = c(1), spec.pcr = 1,
                       sens.lamp = c(.8,.85,.9,.95,1), spec.lamp = 1, lamp.diagnostic = c(F),
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

write.csv(out, "out.csv")

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
df <- data.frame(day = seq(17, by = 7, length.out = 11),
                 new_MSU_cases = c(NA,3,7,66,43,60,65,99,132,212,265),
                 cumulative_MSU_cases = c(38,41,48,114,157,217,282,381,513,725,990))
out %>%
  filter(tests == 0) %>%
  group_by(day) %>%
  summarize(mean = mean(cum.reporting.symptoms.on+cum.reporting.symptoms.off),
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

out %>% 
  filter(day == max(day)) %>% 
  group_by(tests, sens.lamp, day) %>%
  summarize(mean = mean(cum.sum.missed),
            sd = sd(cum.sum.missed),
            upper = mean + 1.96*sd,
            lower = mean - 1.96*sd) %>%
  ggplot(aes(x = factor(sens.lamp), y = mean, color = factor(sens.lamp), 
             fill = factor(sens.lamp), group = factor(sens.lamp))) +
  geom_point() +
  geom_errorbar(aes(ymax = upper, ymin = lower)) +
  facet_grid(.~tests) +
  labs(y="Cumulative Missed Cases \n (95% Of Simulations)",
       x = "Sensitivity") +
  scale_color_discrete(name = "LAMP \nSensitivity",
                       labels = c("80%","85%","90%","95%")) +
  scale_fill_discrete(name = "LAMP \nSensitivity",
                      labels = c("80%","85%","90%","95%")) +
  theme_classic() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
ggsave("LAMP_missed_cases.png")


