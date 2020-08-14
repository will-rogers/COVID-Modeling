library(ggplot2)
library(dplyr)
library(data.table)
library(reshape2)
library(mc2d)
source('sim_functions_onoffbins_wr.R')


# simulate epidemics ------------------------------------------------------
set.seed(12345)
output <- sims_par(tests=c(250, 500, 1000),
                       compliance=c(1,.75,.50), 
                       introductions = c(25),
                       ppn_sympt = .2, 
                       care.seeking = .5,
                       R0.on = c(2,2.5,3),
                       R0.off = c(2,2.5,3),
                       test.scenario = c("5 Days", "4 Days", "3 Days", "2 Days", "1 Day","No Delay"),
                       ncores=NULL) ## ncores not needed for this script at the moment 
                                  # just a demo in case we want to scale up

# output <- output %>% 
#   filter(R0.on == R0.off)

output <- running_tots(output)

# summarize variables -----------------------------------------------------
output$compliance <- factor(output$compliance, levels = c("50%","75%", "100%"), labels = c("50% Compliance","75% Compliance", "100% Compliance"), ordered = T)
output$compliance <- factor(output$compliance, levels=rev(levels(output$compliance)))
# 
output$introductions <- factor(output$introductions)
output$introductions <- factor(output$introductions, levels = c(25, 50, 100), labels = c("25 Initial Cases","50 Initial Cases", "100 Initial Cases"), ordered = T)
output$introductions <- factor(output$introductions, levels=rev(levels(output$introductions)))

output$test.scenario <- factor(output$test.scenario, levels = c("5 Days", "4 Days","3 Days","2 Days","1 Day","No Delay"), labels = c("5 Days", "4 Days","3 Days","2 Days","1 Day","No Delay"), ordered = T)
output$test.scenario <- factor(output$test.scenario, levels=rev(levels(output$test.scenario)))

# output$ppn_sympt <- factor(output$ppn_sympt, levels = c(0.2,0.4,0.6), labels = c("20% Symptoms","40% Symptoms", "60% Symptoms"), ordered = T)
output$R0 <- factor(output$R0.on, levels = c(2,2.5,3), labels = c("R0 = 2","R0 = 2.5", "R0 = 3"), ordered = T)

output$totalrooms <- output$isolation.complying.on + output$quarantine.complying.on
output$new.cases <- output$new.cases.on + output$new.cases.off
output$tot.cases <- output$cum_cases.on + output$cum_cases.off

plotting.final <- output %>% 
  group_by(test.scenario, compliance, R0, tests, day) %>% 
  mutate(mean = mean(totalrooms),
         sd = sd(totalrooms),
         upper = mean + 1.96*sd,
         lower = mean - 1.96*sd)

ggplot(plotting.final %>% filter(tests == 500)  , 
       aes(x=day, y = mean, fill = test.scenario, color = test.scenario)) +
  # geom_bar(position = "dodge", stat = "identity") +
  # geom_errorbar(aes(ymin=lower, ymax=upper))+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, colour = NA) +
  geom_smooth() +
  geom_hline(yintercept = 92, linetype = "dotted") +
  geom_hline(yintercept = 220, linetype = "dotted") +
  facet_grid(R0~compliance) +
  scale_fill_ordinal(name = "Testing Delay", 
                     # guide = F
                     ) +
  scale_color_ordinal(name = "Testing Delay", 
                      # guide = F
                      ) +
  labs(y = "On-Campus Isolation/Quarantine Rooms Per Day",
       x = "Day of School Year") +
  theme_classic() +
  scale_y_continuous(n.breaks = 12) 
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
 # output
# 
# output$new.cases.off
# 
# 
# 
# ggplot(output,aes(x = day, y = cum_cases, color = factor(tests), linetype = factor(R0))) +
#   # geom_line(aes(group = group), alpha = 0.025) +
#   geom_smooth() +
#   facet_grid(compliance~ppn_sympt) +
#   labs(
#     x = "Day of School Year") +
#   # scale_y_log10(n.breaks = 10)+
#   theme_classic()
# 

plotting.final. <- plotting.final %>% 
  filter(test.scenario == "No Delay",
         tests == 500) %>% 
  group_by(compliance, R0) %>%
  summarize(baseline = mean) %>% 
  distinct()

output. <- output %>% 
  filter(tests == 500) 
  

prop.dif <- merge(output., plotting.final., by = c("compliance", "R0"))

a <- prop.dif %>% 
  filter(day == max(day)) %>% 
  group_by(R0, compliance) %>% 
  mutate(prop.dif = (tot.cases-baseline)) %>% 
  filter(test.scenario != "No Delay") %>% 
  group_by(test.scenario, compliance, R0) %>% 
  summarize(mean = mean(prop.dif),
            sd = sd(prop.dif)) %>% 
  mutate(upper = mean + 1.96*sd,
         lower = mean - 1.96*sd)
  
ggplot(a %>% 
         filter(compliance == "100% Compliance"), aes(x = test.scenario, y = mean, fill = test.scenario,
                                                      color = test.scenario)) +
  geom_point(stat = "identity", position = "dodge", size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = "dodge") +
  facet_grid(.~R0) +
  scale_fill_ordinal(guide = F) +
  scale_color_ordinal(guide = F) +
  labs(x = "Testing Delay",
       y =  "# More Cases than 'No Delay'") +
  scale_y_continuous(#labels = scales::percent_format(accuracy = 1),
                     n.breaks = 10) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# ggplot(output, aes(x = day, y = totalrooms, color = factor(R0), linetype=factor(tests))) +
#   geom_smooth() + 
#   facet_grid(compliance~ppn_sympt) +
#   scale_color_discrete(name = bquote(R[0])) +
#   scale_linetype_discrete(name = "Available Tests") +
#   labs(x = "Available Tests Per Day",
#        y = "Isolation & Quarantine Rooms Per Day") +
#   scale_y_log10(n.breaks = 12) +
#   theme_classic()
# 
ggsave("RoomsbyDayR0CompAblines.JPG", dpi=300)
# 
# 
# # plot by compliance ------------------------------------------------------
# 
# plotting.function(output, "cum_cases", log.10 = T)
# plotting.function(output, "cum_reporting.symptoms", log.10 = T)
# plotting.function(output, "cum_all.symptomatics", log.10 = T)
# plotting.function(output, "cum_all.asymptomatics", log.10 = T)
# plotting.function(output, "isolation.complying", log.10 = T)
# plotting.function(output, "quarantine.complying", log.10 = T)
# plotting.function(output, "total_traces", log.10 = T)
# plotting.function(output, "new.cases", log.10 = T)
# 
# plot(ggplot(output,aes(x = day, y = cum_cases, color = factor(tests))) +
#        geom_line(aes(group = group), alpha = 0.025) +
#        geom_smooth() +
#        # facet_grid(compliance~.) +
#        labs(y = "% reduction of Cumulative Cases",
#             x = "Day of School Year",
#             title = "Number of Symptomatic Cases") +
#        scale_y_log10(n.breaks=10)+
#        theme_classic())


