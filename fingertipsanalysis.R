library(tidyverse)

trimres <- read_csv("data/fingertips_trimres.csv")
trimres <- select(trimres, Value, `Time period`, Count, Denominator) %>%
    setNames(., c("perc", "time", "count", "denom")) %>%
    mutate(.,
        prop = count/denom
    )
lstimes <- str_split_fixed(trimres$time, " ", 2)
quarters <- lstimes[,2]
month <- case_when(
    quarters == "Q1" ~ "03",
    quarters == "Q2" ~ "06",
    quarters == "Q3" ~ "09",
    quarters == "Q4" ~ "12",
    .default = as.character(quarters)
)
times <- paste0(lstimes[,1], "/", month, "/", "01")
trimres$time <- as.Date(times, format = "%Y/%m/%d")

png("TrimResistance.png", width = 15, height = 8, units = "cm", res = 300)
p <- ggplot(trimres, aes(time, prop)) + 
    geom_point() + 
    theme_bw() + 
    labs(y = "Proportion of urine samples of E. coli\nresistant to trimethoprim",
        x = "Month")
print(p)
dev.off()

#need a left join not a bindrow
trimressel <- trimres %>%
    ungroup() %>%
    select(time, prop)
#needs to be in days from 0
trimressel$timeindays <- as.numeric(abs(trimressel$time[1] - trimressel$time))

modeldf <- read_csv("modeloutput.csv") %>%
    mutate(., time = trimressel$time[1] + timeindays,
        rep = factor(proputi))

png("TrimResistanceModvsData.png", width = 25, height = 8, units = "cm", res = 300)
p <-ggplot(modeldf, aes(time, rfreq, group = proputi)) + 
    geom_line(col = "darkred", alpha = 0.5) + 
    geom_point(data = trimressel, aes(x = time, y = prop, group = NA, col = NA), col = "black")+
    theme_bw() +
    labs(x = "Time since 2015 in days", 
        y = "Proportion of urine samples of E. coli\nresistant to trimethoprim")
print(p)
dev.off()

average_modresults <- modeldf %>%
    group_by(time) %>%
    summarise(.,
        lwr = quantile(rfreq, 0.025),
        upr = quantile(rfreq, 0.975),
        rfreq = mean(rfreq)
    )

png("TrimResistanceModvsDataRibbon.png", width = 15, height = 8, units = "cm", res = 300)
p <-ggplot(average_modresults, aes(time, rfreq)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr), col = "lightgrey", alpha = 0.3) +
    geom_line(col = "#220f0f", alpha = 0.5) + 
    geom_point(data = trimressel, aes(x = time, y = prop, group = NA, col = NA), col = "black")+
    theme_bw() +
    labs(x = "Time since 2015 in days", 
        y = "Proportion of urine samples of E. coli\nresistant to trimethoprim")
print(p)
dev.off()
