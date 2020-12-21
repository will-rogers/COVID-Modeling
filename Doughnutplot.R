library(readr)
univ <- read_csv("univ.csv")
Sheet <- read_csv("University Testing Programs and Results - Sheet2.csv")
View(univ)

Sheet$Surveillance <- ifelse(str_detect(Sheet$'Asymptomatic Screening?',"es"), "Surveillance Testing", "No Surveillance Testing")

Sheet$Surveillance <- ifelse(Sheet$Online == "Yes", 
                             paste(Sheet$Surveillance, "& Online"),
                             paste(Sheet$Surveillance, "& In-Person"))
Sheet$Surveillance <- factor(Sheet$Surveillance, levels = c("No Surveillance Testing & Online",
                                      "No Surveillance Testing & In-Person", 
                                      "Surveillance Testing & In-Person",     
                                      NA  ), 
       labels = c("No Surveillance Testing & Online",
                  "No Surveillance Testing & In-Person",
                  "Surveillance Testing & In-Person"))

Sheet %>% 
  group_by(Surveillance,Group,Online) %>% 
  summarize(total = n()) %>% 
  mutate(prop = ifelse(Group == "Big Sky", 
                       total/14, NA),
         prop = ifelse(Group == "Land Grant", 
                       total/73, prop)) %>% 
  ggplot(aes(x = 2, y = prop, fill = factor(Surveillance))) +
  geom_bar(stat = "identity", color = "white") +
  scale_fill_discrete(name = "",
                       labels = c("No Surveillance Testing & Online",
                                  "No Surveillance Testing & In-Person",
                                  "Surveillance Testing & In-Person",
                                  "Unspecified"),
                       na.value = "grey") +
  facet_grid(.~Group) +
  xlim(0.5, 2.5) +
  theme(legend.position = "right",
        text = element_text(size = 15))


?scale_fill_manual()


