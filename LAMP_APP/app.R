library(shiny)
library(shinydashboard)
library(scales)
library(tidyverse)
library(data.table)
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
            title = "Epidemic Parameters: R0",
            sliderInput('R0', "R0", value = 3, min = 0, max=5, step = .25),
            solidHeader = TRUE,
            width = 3
        ),
        box(
            title = "Epidemic Parameters: Initial Prevalence",
            sliderInput('init_prev', "Initial Prevalence", value = .005, min = .005, max=.1, step = .005),
            solidHeader = TRUE,
            width = 3
        ),
        box(
            title = "Epidemic Parameters: Immunity",
            sliderInput('immune', "Proportion immune or recovered", value = .05, min = 0, max=.2, step = .025),
            solidHeader = TRUE,
            width = 3
        ),
        box(
            title = "Epidemic Parameters: Population Size",
            sliderInput('N0', "Population Size", value = 16500, min = 0, max=30000, step = 500),
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
        box(title = 'Isolation/Quarantine Demand',
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
            title = "Intervention Parameters: # of tests",
            sliderInput('test_numb', "Number of Tests", value = 2000, min = 0, max=2500, step = 250),
            solidHeader = TRUE,
            width = 4
        ),
        box(
            title = "Intervention Parameters: Test Sensitivity",
            sliderInput('sens', "Test Sensitivity", value = .77, min = .6, max=1, step = .03),
            solidHeader = TRUE,
            width = 4
        ),
        box(
            title = "Intervention Parameters: Test Specificity",
            sliderInput('spec', "Test Specificity", value = .98, min = .9, max=1, step = .01),
            solidHeader = TRUE,
            width = 4
        )
    ),
    fluidRow(
        box(title = 'Cumulative Cases ',
            plotOutput("epi_curve2"),
            width = 3),
        box(title = 'Isolation/Quarantine Demand',
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
            width = 6
        ),
        box(
            title = "PCR Test Cost:",
            sliderInput('pcr_cost', "PCR Test Cost", value = 20, min = 0, max=25, step = .25),
            solidHeader = TRUE,
            width = 6
        )
    ),
    fluidRow(
        box(
            title = "Cost with no intervention: ",
            textOutput("no_strategy"),
            width = 6
        ),
        box(
            title = "Cost with intervention: ",
            textOutput("strategy"),
            width = 6
        )
    )
)

sidebar <- dashboardSidebar(collapsed = T)

ui <- dashboardPage(
    dashboardHeader(title = "Bobcat LAMP: A Safe, Economically Viable Return to Campus",
                    titleWidth = 800),
    sidebar,
    body
)

server <- function(input, output){
    no_intervention <- reactive({
       uni_sims_par(tst = c(0),
                             test.timeline = c("Sustained"),
                             compliance = c(1),
                             init.prev = input$init_prev,
                             ppn_sympt = c(0.35),
                             care.seeking = c(1),
                             R0.on = input$R0,
                             R0.off = input$R0,
                             test.scenario = c("1 Day"),
                             sens.pcr = 0.99,
                             spec.pcr = 0.99,
                             sens.lamp = c(.925),
                             spec.lamp = 0.98,
                             lamp.diagnostic = c(F),
                             community.intro.daily.on = 1,
                             community.prob.daily.on = c(0.1),
                             community.intro.daily.off = 1,
                             community.prob.daily.off = c(0.1),
                             immunity = input$immune,
                             N0 = input$N0,
                             on.campus.prop = .25,
                             contact.tracing.limit = c(25),
                             pooling = c(1),
                             pooling.multi = c(.5),
                             days = 100,
                             sims = 25,
                             ncores=1)
    })
    
    intervention <- reactive({
        uni_sims_par(tst = input$test_numb,
                     test.timeline = c("Sustained"),
                     compliance = c(1),
                     init.prev = input$init_prev,
                     ppn_sympt = c(0.35),
                     care.seeking = c(1),
                     R0.on = input$R0,
                     R0.off = input$R0,
                     test.scenario = c("1 Day"),
                     sens.pcr = 0.99,
                     spec.pcr = 0.99,
                     sens.lamp = input$sens,
                     spec.lamp = input$spec,
                     lamp.diagnostic = c(F),
                     community.intro.daily.on = 1,
                     community.prob.daily.on = c(0.1),
                     community.intro.daily.off = 1,
                     community.prob.daily.off = c(0.1),
                     immunity = input$immune,
                     N0 = input$N0,
                     on.campus.prop = .25,
                     contact.tracing.limit = c(25),
                     pooling = c(1),
                     pooling.multi = c(.5),
                     days = 100,
                     sims = 25,
                     ncores=1)
    })
    
    output$epi_curve <- renderPlot(
       
        no_intervention() %>%
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
                   cum.iso.on = cumsum(isolation.complying.on),
                   cum.iso.off = cumsum(isolation.complying.off),
                   cum.qua.on = cumsum(quarantine.complying.on),
                   cum.qua.off = cumsum(quarantine.complying.off),
                   pcr.demand.sym = cumsum(symp.pcr), 
                   pcr.demand.asym = cumsum(asymp.pcr)
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
            ggplot(aes(x = day, y = (symp.pcr + asymp.pcr) * 20, group = factor(group))) +
            geom_line(alpha = .4) + theme_bw() + ylab("Number of Students") + labs(caption = 'assumes 5% positivity')
    )
    output$class_days <- renderPlot(
        no_intervention() %>%
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
                   cum.iso.on = cumsum(isolation.complying.on),
                   cum.iso.off = cumsum(isolation.complying.off),
                   cum.qua.on = cumsum(quarantine.complying.on),
                   cum.qua.off = cumsum(quarantine.complying.off),
                   pcr.demand.sym = cumsum(symp.pcr), 
                   pcr.demand.asym = cumsum(asymp.pcr)
            ) %>%
            ggplot(aes(x = day, y = cum.iso.on + cum.iso.off + cum.qua.on + cum.qua.off, group = factor(group))) +
            geom_line(alpha = .4) + theme_bw() + ylab("Number of Student Class Days Missed \n Due to Isolation or Quarantine")
    )
    
    output$epi_curve2 <- renderPlot(
        
        intervention() %>%
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
                   cum.iso.on = cumsum(isolation.complying.on),
                   cum.iso.off = cumsum(isolation.complying.off),
                   cum.qua.on = cumsum(quarantine.complying.on),
                   cum.qua.off = cumsum(quarantine.complying.off),
                   pcr.demand.sym = cumsum(symp.pcr), 
                   pcr.demand.asym = cumsum(asymp.pcr)
            ) %>%
            ggplot(aes(x = day, y = cum.cases.on + cum.cases.off, group = factor(group))) +
            geom_line(alpha = .4) + theme_bw() + ylab("Number of Students") + ylim(NA, input$N0)
        
    )
    output$iso2 <- renderPlot(
        
        intervention() %>%
            ggplot(aes(x = day, y = isolation.complying.on + quarantine.complying.on, group = factor(group))) +
            geom_line(alpha = .4) + theme_bw() + ylab("Number of Students")
        
    )
    output$pcr2 <- renderPlot(
        intervention() %>% 
            ggplot(aes(x = day, y = (symp.pcr + asymp.pcr) * 20, group = factor(group))) +
            geom_line(alpha = .4) + theme_bw() + ylab("Number of Students") + labs(caption = 'assumes 5% positivity')
    )
    output$class_days2 <- renderPlot(
        intervention() %>%
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
                   cum.iso.on = cumsum(isolation.complying.on),
                   cum.iso.off = cumsum(isolation.complying.off),
                   cum.qua.on = cumsum(quarantine.complying.on),
                   cum.qua.off = cumsum(quarantine.complying.off),
                   pcr.demand.sym = cumsum(symp.pcr), 
                   pcr.demand.asym = cumsum(asymp.pcr)
            ) %>%
            ggplot(aes(x = day, y = cum.iso.on + cum.iso.off + cum.qua.on + cum.qua.off, group = factor(group))) +
            geom_line(alpha = .4) + theme_bw() + ylab("Number of Student Class Days Missed \n Due to Isolation or Quarantine")
    )
    output$no_strategy <- renderText(
        paste("The total cost of PCR testing is",
        no_intervention() %>%
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
                   cum.iso.on = cumsum(isolation.complying.on),
                   cum.iso.off = cumsum(isolation.complying.off),
                   cum.qua.on = cumsum(quarantine.complying.on),
                   cum.qua.off = cumsum(quarantine.complying.off),
                   pcr.demand.sym = cumsum(symp.pcr), 
                   pcr.demand.asym = cumsum(asymp.pcr),
                   total.pcr = pcr.demand.sym * 20 + pcr.demand.asym
            ) %>% group_by(group) %>%
            summarise(value=last(total.pcr)) %>% ungroup() %>%
            summarise(tests = mean(value) ) %>% 
            mutate(total_cost = tests * input$pcr_cost) %>% pull()  %>% dollar_format()()
    , "on average"))
    output$strategy <- renderText(
        paste("The total cost of PCR testing + LAMP is",
              intervention() %>%
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
                         cum.iso.on = cumsum(isolation.complying.on),
                         cum.iso.off = cumsum(isolation.complying.off),
                         cum.qua.on = cumsum(quarantine.complying.on),
                         cum.qua.off = cumsum(quarantine.complying.off),
                         pcr.demand.sym = cumsum(symp.pcr), 
                         pcr.demand.asym = cumsum(asymp.pcr),
                         total.pcr = pcr.demand.sym * 20 + pcr.demand.asym
                  ) %>% group_by(group) %>%
                  summarise(value=last(total.pcr)) %>% ungroup() %>%
                  summarise(tests = mean(value)) %>% 
                  mutate(total_cost = tests * input$pcr_cost + input$lamp_cost * 100 * input$test_numb) %>% pull()  %>% dollar_format()()
        , "on average"))
}


##########################
# Run App
##########################
shinyApp(ui = ui, server = server)
