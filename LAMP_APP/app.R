library(shiny)
library(shinydashboard)
library(scales)
library(tidyverse)
library(data.table)
library(combinat)
library(splines)
source("sir_lamp_function.R")
source("sir_simple_step_function.R")
source("uni_sim_function.R")
source("uni_sim_par_function.R")
source("plotting_functions.R")


body <- dashboardBody(
    fluidRow(
        box(title = "Part 1: Set Model Parameters for Simulated Epidemic", 
            status = "primary", 
            width = 12)
    ),
    fluidRow(
        box(
            title = "Epidemic Parameters: R0 on campus",
            sliderInput('R0.on', "R0 On", value = 3, min = 0, max=8, step = .25),
            solidHeader = TRUE,
            width = 3
        ),
        box(
            title = "Epidemic Parameters: R0 off campus",
            sliderInput('R0.off', "R0 Off", value = 3, min = 0, max=8, step = .25),
            solidHeader = TRUE,
            width = 3
        ),
        box(
            title = "Epidemic Parameters: Exposure Period",
            sliderInput('exposure.days', "Exposure Period", value = 5, min = 0, max=10, step = 1),
            solidHeader = TRUE,
            width = 3
        ),
        box(
            title = "Epidemic Parameters: Presymptomatic Period",
            sliderInput('presymptom.days', "Presymptom Period", value = 2, min = 0, max=10, step = 1),
            solidHeader = TRUE,
            width = 3
        ),
        box(
            title = "Epidemic Parameters: Postsymptomatic Period",
            sliderInput('postsymptom.days', "Postsymptom Period", value = 7, min = 0, max=15, step = 1),
            solidHeader = TRUE,
            width = 3
        ),
        box(
            title = "Epidemic Parameters: Initial Prevalence",
            sliderInput('init.prev', "Initial Prevalence", value = .03, min = .01, max=0.5, step = .01),
            solidHeader = TRUE,
            width = 3
        ),
        box(
            title = "Epidemic Parameters: Immunity",
            sliderInput('immunity', "Proportion immune or recovered", value = .15, min = 0, max=1, step = .01),
            solidHeader = TRUE,
            width = 3
        ),
        box(
            title = "Epidemic Parameters: Proportion Symptomatic ",
            sliderInput('ppn_sympt', "Population Size", value = 0.65, min = 0, max=1, step = 0.01),
            solidHeader = TRUE,
            width = 3
        ),
        box(
            title = "Epidemic Parameters: Proportion of Symptomatic Seeking Care",
            sliderInput('care.seeking', "Symptomatic Seeking Care", value = 1, min = 0, max=1, step = 0.01),
            solidHeader = TRUE,
            width = 3
        ),
        box(
            title = "Epidemic Parameters: Population Size",
            sliderInput('N0', "Population Size", value = 16500, min = 0, max=30000, step = 500),
            solidHeader = TRUE,
            width = 3
        ),
        box(
            title = "Epidemic Parameters: Introduction Rate On-campus",
            sliderInput('community.prob.daily.on', "Intro. Prob", value = 0.1, min = 0, max=1, step = 0.01),
            solidHeader = TRUE,
            width = 3
        ),
        box(
            title = "Epidemic Parameters: Introduction Size On-campus",
            sliderInput('community.intro.daily.on', "Population Size", value = 1, min = 0, max=50, step = 1),
            solidHeader = TRUE,
            width = 3
        ),
        box(
            title = "Epidemic Parameters: Introduction Rate Off-campus",
            sliderInput('community.prob.daily.off', "Intro. Size", value = 0.1, min = 0, max=1, step = 0.01),
            solidHeader = TRUE,
            width = 3
        ),
        box(
            title = "Epidemic Parameters: Introduction Size Off-campus",
            sliderInput('community.intro.daily.off', "Intro. Size", value = 1, min = 0, max=50, step = 1),
            solidHeader = TRUE,
            width = 3
        ),
        box(
            title = "Epidemic Parameters: Proportion on Campus",
            sliderInput('on.campus.prop', "Prop. on Campus", value = 0.25, min = 0, max=1, step = 0.01),
            solidHeader = TRUE,
            width = 3
        ),
        box(
            title = "Epidemic Parameters: Days to simulate",
            sliderInput('days', "Days", value = 100, min = 2, max=200, step = 1),
            solidHeader = TRUE,
            width = 3
        ),
        box(
            title = "Epidemic Parameters: Simulations",
            sliderInput('sims', "Simulations", value = 25, min = 2, max=200, step = 1),
            solidHeader = TRUE,
            width = 3
        )
    ),
    fluidRow(
        box(title = "Simulated Epidemic with no Surveillance Testing", 
            
            width = 12)
        
    ),
    fluidRow(
        box(title = 'Cumulative Cases ',
            plotOutput("epi_curve"),
            width = 3),
        box(title = 'On Campus Isolation/Quarantine Demand',
            plotOutput("iso"),
            width = 3),
        box(title = 'PCR Demand',
            plotOutput("pcr"),
            width = 3),
        box(title = 'Cumulative Class Days Missed',
            plotOutput("class_days"),
            width = 3)
    ),
    fluidRow(
        box(title = "Part 2: Choose Surveillance Intervention", 
            status = "primary", 
            width = 12)
    ),
    fluidRow(
        box(
            title = "Intervention Parameters: Screening Tests",
            sliderInput('tst', "Number of Screening Tests", value = 2500, min = 0, max=5000, step = 50),
            solidHeader = TRUE,
            width = 4
        ),
        box(
            title = "Intervention Parameters: Screening Test Sensitivity",
            sliderInput('sens.lamp', "Screening Sensitivity", value = 0.75, min = 0, max=1, step = 0.01),
            solidHeader = TRUE,
            width = 4
        ),
        box(
            title = "Intervention Parameters: Screening Test Specificity",
            sliderInput('spec.lamp', "Screening Specificity", value = 0.98, min = 0, max=1, step = 0.01),
            solidHeader = TRUE,
            width = 4
        ),
        box(
            title = "Intervention Parameters: PCR Sensitivity",
            sliderInput('sens.pcr', "PCR Sensitivity", value = 0.99, min = 0, max=1, step = 0.01),
            solidHeader = TRUE,
            width = 4
        ),
        box(
            title = "Intervention Parameters: PCR Specificity",
            sliderInput('spec.pcr', "PCR Specificity", value = 0.99, min = 0, max=1, step = 0.01),
            solidHeader = TRUE,
            width = 4
        ),
        box(
            title = "Intervention Parameters: Test Timeline",
            
            selectInput('test.timeline', "Test Timeline", choices = c("Initial", "Sustained", "Both"), selected = "Initial"),
            solidHeader = TRUE,
            width = 4
        ),
        box(
            title = "Intervention Parameters: Test Delay",
            
            selectInput('test.scenario', "Test Delay", choices = c("2 Days","1 Day","No Delay"), selected = "1 Day"),
            solidHeader = TRUE,
            width = 4
        ),
        box(
            title = "Intervention Parameters: Does the screening test serve as a diagnostic?",
            
            selectInput('lamp.diagnostic', "Screening Test Legally Appropriate", choices = c(TRUE, FALSE), selected = FALSE),
            solidHeader = TRUE,
            width = 4
        ),
        box(
            title = "Intervention Parameters: Compliance",
            sliderInput('compliance', "Compliance", value = 1, min = 0, max=1, step = 0.01),
            solidHeader = TRUE,
            width = 4
        ),
        box(
            title = "Intervention Parameters: Contact Tracing Limit",
            sliderInput('contact.tracing.limit', "Contact Tracing Limit", value = 25, min = 0, max=500, step = 1),
            solidHeader = TRUE,
            width = 4
        ),
        box(
            title = "Intervention Parameters: Days to Isolate",
            sliderInput('days.to.isolate', "Days to Isolate", value = 25, min = 0, max=500, step = 1),
            solidHeader = TRUE,
            width = 4
        ),
        box(
            title = "Intervention Parameters: Days to Quarantine",
            sliderInput('days.to.quarantine', "Days to Quarantine", value = 25, min = 0, max=500, step = 1),
            solidHeader = TRUE,
            width = 4
        )
    ),
    fluidRow(
        box(
            title = "LAMP Pool Size: ",
            sliderInput('pooling', "LAMP Pool Size", value = 1, min = 1, max=10, step = 1),
            solidHeader = TRUE,
            width = 6
        ),
        box(
            title = "Pooled Sensitivity Reduction: ",
            sliderInput('pooling.multi', "Pooled Sensitivity Reduction", value = 0, min = 0, max=.1, step = .01),
            solidHeader = TRUE,
            width = 6
        )
    ),
    fluidRow(
        box(title = 'Cumulative Cases ',
            plotOutput("epi_curve2"),
            width = 3),
        box(title = 'On Campus Isolation/Quarantine Demand',
            plotOutput("iso2"),
            width = 3),
        box(title = 'PCR Demand',
            plotOutput("pcr2"),
            width = 3),
        box(title = 'Cumulative Class Days Missed',
            plotOutput("class_days2"),
            width = 3)
    ),
    fluidRow(
        box(title = "Part 3: Strategy Comparison", 
            status = "primary", 
            width = 12)
    ),
    fluidRow(
        box(
            title = "LAMP Test Cost: ",
            sliderInput('lamp_cost', "LAMP Test Cost", value = 3.5, min = 0, max=20, step = .25),
            solidHeader = TRUE,
            width = 4
        ),
        box(
            title = "PCR Test Cost:",
            sliderInput('pcr_cost', "PCR Test Cost", value = 12.5, min = 0, max=25, step = .25),
            solidHeader = TRUE,
            width = 4
        ),
        box(
            title = "Positivity Rate:",
            sliderInput('pos_rate', "Symptomatic PCR Test positivity rate", value = .2, min = 0, max=.50, step = .01),
            solidHeader = TRUE,
            width = 4,
            'note this is used to derive the total number of PCR tests'
        )
    ),
    fluidRow(
        box(
            title = "Quarantine / Isolation Cost: ",
            sliderInput('iso_cost', "Quarantine / Isolation Cost", value = 50, min = 0, max=200, step = 10),
            solidHeader = TRUE,
            width = 6,
            'defined as the cost of housing a student in quarantine/isolation for a day'
            
        ),
        box(
            title = "Missed Class Cost: ",
            sliderInput('missed_cost', "Cost (Value) of a missed class", value = 50, min = 0, max=200, step = 10),
            solidHeader = TRUE,
            width = 6,
            'defined as a day student is unable to attend class in person'
        )
    ),
    fluidRow(
        box(
            title = "No surveillance intervention: ",
            textOutput("no_strategy_size"),
            textOutput("no_strategy"),
            textOutput("nostrategy_iso"),
            textOutput("nostrategy_iso_cost"),
            textOutput("nostrategy_class_cost"),
            textOutput("nostrategy_total_cost"),
            width = 6
        ),
        box(
            title = "Cost with intervention: ",
            textOutput("strategy_size"),
            textOutput("strategy"),
            textOutput("strategy_iso"),
            textOutput("strategy_iso_cost"),
            textOutput("strategy_class_cost"),
            textOutput("strategy_total_cost"),
            width = 6
        )
    )
)

sidebar <- dashboardSidebar(collapsed = T)

ui <- dashboardPage(
    dashboardHeader(title = "Bobcat LAMP: A Safe, Economically Viable Return to Campus",
                    titleWidth = 800),
    sidebar,
    # tabPanel("ModelDescription",
    #          withMathJax(includeHTML("description_combo.html"))
    # ),
    body
)

server <- function(input, output){
    no_intervention <- reactive({
        uni_sims_par(tst = 0, 
                     test.timeline = input$test.timeline,
                     compliance = input$compliance, 
                     init.prev = input$init.prev, 
                     ppn_sympt = input$ppn_sympt, 
                     care.seeking = input$care.seeking, 
                     R0.on = input$R0.on,  
                     R0.off = input$R0.off, 
                     test.scenario = input$test.scenario,
                     sens.pcr = input$sens.pcr, 
                     spec.pcr = input$spec.pcr, 
                     sens.lamp = input$sens.lamp, 
                     spec.lamp = input$spec.lamp, 
                     lamp.diagnostic = input$lamp.diagnostic, 
                     community.intro.daily.on = input$community.intro.daily.on, 
                     community.prob.daily.on = input$community.prob.daily.on,
                     community.intro.daily.off = input$community.intro.daily.off, 
                     community.prob.daily.off = input$community.prob.daily.off,
                     immunity = input$immunity, 
                     N0 = input$N0, 
                     on.campus.prop = input$on.campus.prop, 
                     contact.tracing.limit = input$contact.tracing.limit,
                     pooling = input$pooling, 
                     pooling.multi = input$pooling.multi,
                     days = input$days, 
                     sims = input$sims,
                     days.to.isolate = input$days.to.isolate,
                     days.to.quarantine = input$days.to.quarantine,
                     exposure.days = input$exposure.days, 
                     presymptom.days = input$presymptom.days, 
                     postsymptom.days = input$postsymptom.days,
                     ncores=1)
    })
    
    intervention <- reactive({
        uni_sims_par(tst = input$tst, 
                     test.timeline = input$test.timeline,
                     compliance = input$compliance, 
                     init.prev = input$init.prev, 
                     ppn_sympt = input$ppn_sympt, 
                     care.seeking = input$care.seeking, 
                     R0.on = input$R0.on,  
                     R0.off = input$R0.off, 
                     test.scenario = input$test.scenario,
                     sens.pcr = input$sens.pcr, 
                     spec.pcr = input$spec.pcr, 
                     sens.lamp = input$sens.lamp, 
                     spec.lamp = input$spec.lamp, 
                     lamp.diagnostic = input$lamp.diagnostic, 
                     community.intro.daily.on = input$community.intro.daily.on, 
                     community.prob.daily.on = input$community.prob.daily.on,
                     community.intro.daily.off = input$community.intro.daily.off, 
                     community.prob.daily.off = input$community.prob.daily.off,
                     immunity = input$immunity, 
                     N0 = input$N0, 
                     on.campus.prop = input$on.campus.prop, 
                     contact.tracing.limit = input$contact.tracing.limit,
                     pooling = input$pooling, 
                     pooling.multi = input$pooling.multi,
                     days = input$days, 
                     sims = input$sims,
                     days.to.isolate = input$days.to.isolate,
                     days.to.quarantine = input$days.to.quarantine,
                     exposure.days = input$exposure.days, 
                     presymptom.days = input$presymptom.days, 
                     postsymptom.days = input$postsymptom.days,
                     ncores=1)
    })
    
    output$epi_curve <- renderPlot(
        
        no_intervention() %>%
            group_by(group) %>%
            mutate(cum.cases.on = cumsum(new.cases.on),
                   cum.cases.off = cumsum(new.cases.off)
            ) %>%
            ggplot(aes(x = day, y = cum.cases.on + cum.cases.off, group = factor(group))) +
            geom_line(alpha = .4) + theme_bw() + ylab("Number of Students") + ylim(NA, input$N0)
        
    )
    output$iso <- renderPlot(
        
        no_intervention() %>%
            ggplot(aes(x = day, y = isolation.complying.on + quarantine.complying.on, group = factor(group))) +
            geom_line(alpha = .4) + theme_bw() + ylab("Number of Students")
        
    )
    output$pcr <- renderPlot(
        no_intervention() %>% 
            ggplot(aes(x = day, y = (symp.pcr + asymp.pcr) * (1 / input$pos_rate), group = factor(group))) +
            geom_line(alpha = .4) + theme_bw() + ylab("Number of Tests") + labs(caption = 'total derived with percent positivity')
    )
    output$class_days <- renderPlot(
        no_intervention() %>%
            group_by(group) %>%
            mutate(cum.iso.on = cumsum(isolation.complying.on),
                   cum.iso.off = cumsum(isolation.complying.off),
                   cum.qua.on = cumsum(quarantine.complying.on),
                   cum.qua.off = cumsum(quarantine.complying.off)
            ) %>%
            ggplot(aes(x = day, y = cum.iso.on + cum.iso.off + cum.qua.on + cum.qua.off, group = factor(group))) +
            geom_line(alpha = .4) + theme_bw() + ylab("Number of Student Class Days Missed \n Due to Isolation or Quarantine") +
            scale_y_continuous(sec.axis = sec_axis(~ . / input$N0, name = "Average Number of Class Days Missed per Student \n Due to Isolation or Quarantine" ))
    )
    
    output$epi_curve2 <- renderPlot(
        
        intervention() %>%
            group_by(group) %>%
            mutate(cum.cases.on = cumsum(new.cases.on),
                   cum.cases.off = cumsum(new.cases.off)
            ) %>%
            ggplot(aes(x = day, y = cum.cases.on + cum.cases.off, group = factor(group))) +
            geom_line(alpha = .4) + theme_bw() + ylab("Number of Students") + ylim(NA, input$N0)
        
    )
    output$iso2 <- renderPlot(
        
        intervention() %>%
            ggplot(aes(x = day, y = isolation.complying.on + quarantine.complying.on, group = factor(group))) +
            geom_line(alpha = .4) + theme_bw() + ylab("Number of Students") + 
            ylim(0, no_intervention() %>% mutate(y = isolation.complying.on + quarantine.complying.on ) %>% summarize(max_y = max(y)) %>% select(max_y) %>% pull())
        
    )
    output$pcr2 <- renderPlot(
        intervention() %>% 
            ggplot(aes(x = day, y = (symp.pcr + asymp.pcr) * (1 / input$pos_rate), group = factor(group))) +
            geom_line(alpha = .4) + theme_bw() + ylab("Number of Tests") + labs(caption = 'total derived with percent positivity') +
            ylim(0, no_intervention() %>% mutate(y = (symp.pcr + asymp.pcr) * (1 / input$pos_rate) ) %>% summarize(max_y = max(y)) %>% select(max_y) %>% pull())
    )
    output$class_days2 <- renderPlot(
        intervention() %>%
            group_by(group) %>%
            mutate(cum.iso.on = cumsum(isolation.complying.on),
                   cum.iso.off = cumsum(isolation.complying.off),
                   cum.qua.on = cumsum(quarantine.complying.on),
                   cum.qua.off = cumsum(quarantine.complying.off),
                   type = 'intervention'
            )  %>%
            ggplot(aes(x = day, y = cum.iso.on + cum.iso.off + cum.qua.on + cum.qua.off, group = factor(group))) +
            geom_path(alpha = .4) + theme_bw() + ylab("Number of Student Class Days Missed \n Due to Isolation or Quarantine") +
            
            scale_y_continuous(limits = c(0,no_intervention() %>%
                                              group_by(group) %>%
                                              mutate(cum.iso.on = cumsum(isolation.complying.on),
                                                     cum.iso.off = cumsum(isolation.complying.off),
                                                     cum.qua.on = cumsum(quarantine.complying.on),
                                                     cum.qua.off = cumsum(quarantine.complying.off),
                                                     y = cum.iso.on + cum.iso.off + cum.qua.on + cum.qua.off) %>%
                                              ungroup() %>% summarise(y = max(y)) %>% select(y) %>% pull()), sec.axis = sec_axis(~ . / input$N0, name = "Average Number of Class Days Missed per Student \n Due to Isolation or Quarantine" ))
    )
    output$no_strategy <- renderText(
        paste("The total cost of PCR testing is",
              no_intervention() %>%
                  group_by(group) %>%
                  mutate(pcr.demand.sym = cumsum(symp.pcr), 
                         pcr.demand.asym = cumsum(asymp.pcr),
                         total.pcr = pcr.demand.sym * (1 / input$pos_rate) + pcr.demand.asym
                  ) %>% group_by(group) %>%
                  summarise(value=last(total.pcr)) %>% ungroup() %>%
                  summarise(tests = mean(value) ) %>% 
                  mutate(total_cost = tests * input$pcr_cost) %>% pull()  %>% dollar_format()()
              , "on average"))
    output$strategy <- renderText(
        paste("The total cost of PCR testing + LAMP is",
              intervention() %>%
                  group_by(group) %>%
                  mutate(pcr.demand.sym = cumsum(symp.pcr), 
                         pcr.demand.asym = cumsum(asymp.pcr),
                         total.pcr = pcr.demand.sym * (1 / input$pos_rate) + pcr.demand.asym
                  ) %>% group_by(group) %>%
                  summarise(value=last(total.pcr)) %>% ungroup() %>%
                  summarise(tests = mean(value)) %>% 
                  mutate(total_cost = tests * input$pcr_cost + input$lamp_cost * 100 * input$tst) %>% pull()  %>% dollar_format()()
              , "on average"))
    output$no_strategy_size <- renderText(
        paste("The total number of positive COVID-19 cases is",
              no_intervention() %>%
                  group_by(group) %>%
                  mutate(cum.cases.on = cumsum(new.cases.on),
                         cum.cases.off = cumsum(new.cases.off),
                         cases = cum.cases.on + cum.cases.off
                  ) %>% group_by(group) %>%
                  summarise(value=last(cases)) %>% ungroup() %>%
                  summarise(tests = mean(value)) %>% pull() %>% round()
              , "on average"))
    output$strategy_size <- renderText(
        paste("The total number of positive COVID-19 cases is",
              intervention() %>%
                  group_by(group) %>%
                  mutate(cum.cases.on = cumsum(new.cases.on),
                         cum.cases.off = cumsum(new.cases.off),
                         cases = cum.cases.on + cum.cases.off
                  ) %>% group_by(group) %>%
                  summarise(value=last(cases)) %>% ungroup() %>%
                  summarise(tests = mean(value)) %>% pull() %>% round()
              , "on average"))
    output$strategy_iso <- renderText(
        paste("The total number of on campus student days spent in quarantine or isolation is",
              round(intervention() %>%
                        group_by(group) %>%
                        mutate(cum.iso.on = cumsum(isolation.complying.on),
                               cum.qua.on = cumsum(quarantine.complying.on),
                               iso = cum.iso.on + cum.qua.on
                        ) %>% group_by(group) %>%
                        summarise(value=last(iso)) %>% ungroup() %>%
                        summarise(tests = mean(value)) %>% pull()  / input$N0 / .25,1),
              'per on campus student'
        )
    )
    output$nostrategy_iso <- renderText(
        paste("The total number of on campus student days spent in quarantine or isolation is",
              round(no_intervention() %>%
                        group_by(group) %>%
                        mutate(cum.iso.on = cumsum(isolation.complying.on),
                               cum.qua.on = cumsum(quarantine.complying.on),
                               iso = cum.iso.on + cum.qua.on
                        ) %>% group_by(group) %>%
                        summarise(value=last(iso)) %>% ungroup() %>%
                        summarise(tests = mean(value)) %>% pull()  / input$N0 / .25,1),
              'per on campus student'
        )
    )
    output$strategy_iso_cost <- renderText(
        paste("The total cost of on campus student days spent in quarantine or isolation is",
              intervention() %>%
                  group_by(group) %>%
                  mutate(cum.iso.on = cumsum(isolation.complying.on),
                         cum.qua.on = cumsum(quarantine.complying.on),
                         iso = cum.iso.on + cum.qua.on
                  ) %>% group_by(group) %>%
                  summarise(value=last(iso)) %>% ungroup() %>%
                  summarise(tests = mean(value)* input$iso_cost) %>% pull() %>% round() %>% dollar_format()()
        )
    )
    output$nostrategy_iso_cost <- renderText(
        paste("The total cost of on campus student days spent in quarantine or isolation is",
              no_intervention() %>%
                  group_by(group) %>%
                  mutate(cum.iso.on = cumsum(isolation.complying.on),
                         cum.qua.on = cumsum(quarantine.complying.on),
                         iso = cum.iso.on + cum.qua.on
                  ) %>% group_by(group) %>%
                  summarise(value=last(iso)) %>% ungroup() %>%
                  summarise(tests = mean(value)* input$iso_cost) %>% pull() %>% round() %>% dollar_format()()
        )
    )
    output$nostrategy_class_cost <- renderText(
        paste("The total cost of missed class days is",
              no_intervention() %>%
                  group_by(group) %>%
                  mutate(cum.iso.on = cumsum(isolation.complying.on),
                         cum.iso.off = cumsum(isolation.complying.off),
                         cum.qua.on = cumsum(quarantine.complying.on),
                         cum.qua.off = cumsum(quarantine.complying.off),
                         missed = cum.iso.on + cum.qua.on + cum.iso.off + cum.qua.off
                  ) %>% group_by(group) %>%
                  summarise(value=last(missed)) %>% ungroup() %>%
                  summarise(tests = mean(value) * input$missed_cost) %>% pull() %>% round() %>% dollar_format()()
        )
    )
    output$strategy_class_cost <- renderText(
        paste("The total cost of missed class days is",
              intervention() %>%
                  group_by(group) %>%
                  mutate(cum.iso.on = cumsum(isolation.complying.on),
                         cum.iso.off = cumsum(isolation.complying.off),
                         cum.qua.on = cumsum(quarantine.complying.on),
                         cum.qua.off = cumsum(quarantine.complying.off),
                         missed = cum.iso.on + cum.qua.on + cum.iso.off + cum.qua.off
                  ) %>% group_by(group) %>%
                  summarise(value=last(missed)) %>% ungroup() %>%
                  summarise(tests = mean(value) * input$missed_cost) %>% pull() %>% round() %>% dollar_format()()
        )
    )
    output$strategy_total_cost <- renderText(
        paste("The implied total cost of testing, quarantine, missed class days is",
              (intervention() %>%
                   group_by(group) %>%
                   mutate(cum.iso.on = cumsum(isolation.complying.on),
                          cum.iso.off = cumsum(isolation.complying.off),
                          cum.qua.on = cumsum(quarantine.complying.on),
                          cum.qua.off = cumsum(quarantine.complying.off),
                          missed = cum.iso.on + cum.qua.on + cum.iso.off + cum.qua.off
                   ) %>% group_by(group) %>%
                   summarise(value=last(missed)) %>% ungroup() %>%
                   summarise(tests = mean(value) * input$missed_cost) %>% pull() %>% round()  +
                   intervention() %>%
                   group_by(group) %>%
                   mutate(cum.iso.on = cumsum(isolation.complying.on),
                          cum.qua.on = cumsum(quarantine.complying.on),
                          iso = cum.iso.on + cum.qua.on
                   ) %>% group_by(group) %>%
                   summarise(value=last(iso)) %>% ungroup() %>%
                   summarise(tests = mean(value)* input$iso_cost) %>% pull() %>% round() + 
                   intervention() %>%
                   group_by(group) %>%
                   mutate(pcr.demand.sym = cumsum(symp.pcr), 
                          pcr.demand.asym = cumsum(asymp.pcr),
                          total.pcr = pcr.demand.sym * (1 / input$pos_rate) + pcr.demand.asym
                   ) %>% group_by(group) %>%
                   summarise(value=last(total.pcr)) %>% ungroup() %>%
                   summarise(tests = mean(value)) %>% 
                   mutate(total_cost = tests * input$pcr_cost + input$lamp_cost * 100 * input$tst) %>% pull())
              %>% dollar_format()()
        )
    )
    output$nostrategy_total_cost <- renderText(
        paste("The implied total cost of testing, quarantine, missed class days is",
              (no_intervention() %>%
                   group_by(group) %>%
                   mutate(cum.iso.on = cumsum(isolation.complying.on),
                          cum.iso.off = cumsum(isolation.complying.off),
                          cum.qua.on = cumsum(quarantine.complying.on),
                          cum.qua.off = cumsum(quarantine.complying.off),
                          missed = cum.iso.on + cum.qua.on + cum.iso.off + cum.qua.off
                   ) %>% group_by(group) %>%
                   summarise(value=last(missed)) %>% ungroup() %>%
                   summarise(tests = mean(value) * input$missed_cost) %>% pull() %>% round()  +
                   no_intervention() %>%
                   group_by(group) %>%
                   mutate(cum.iso.on = cumsum(isolation.complying.on),
                          cum.qua.on = cumsum(quarantine.complying.on),
                          iso = cum.iso.on + cum.qua.on
                   ) %>% group_by(group) %>%
                   summarise(value=last(iso)) %>% ungroup() %>%
                   summarise(tests = mean(value)* input$iso_cost) %>% pull() %>% round() + 
                   no_intervention() %>%
                   group_by(group) %>%
                   mutate(pcr.demand.sym = cumsum(symp.pcr), 
                          pcr.demand.asym = cumsum(asymp.pcr),
                          total.pcr = pcr.demand.sym * (1 / input$pos_rate) + pcr.demand.asym
                   ) %>% group_by(group) %>%
                   summarise(value=last(total.pcr)) %>% ungroup() %>%
                   summarise(tests = mean(value)) %>% 
                   mutate(total_cost = tests * input$pcr_cost ) %>% pull())
              %>% dollar_format()()
        )
    )
}


##########################
# Run App
##########################
shinyApp(ui = ui, server = server)

