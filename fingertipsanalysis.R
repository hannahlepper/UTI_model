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

modeldf <- read_csv("modeloutput.csv")

png("TrimResistanceModvsData.png", width = 15, height = 8, units = "cm", res = 300)
p <-ggplot(modeldf, aes(timeindays, rfreq)) + 
    geom_line(col = "darkred") + 
    geom_point(data = trimressel, aes(x = timeindays, y = prop))+
    theme_bw() +
    labs(x = "Time since 2015 in days", 
        y = "Proportion of urine samples of E. coli\nresistant to trimethoprim")
print(p)
dev.off()
