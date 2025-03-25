library(tidyverse)

pres <- read_csv("data/foi_waterlow.csv")

trim <- subset(pres, drug_name == "Trimethoprim") %>%
    group_by(YEAR, MONTH, GENDER) %>%
    summarise(.,
        average_usage = mean(ITEMS/population),
        av_usage_items = mean(ITEMS)
    ) %>%
    mutate(., 
        month_year = paste0(YEAR, "/", MONTH, "/1"),
        monthdate = as.Date(month_year, format = "%Y/%m/%d"),
        ab = "trim"
    )

nitro <- subset(pres, drug_name == "Nitrofurantoin") %>%
    group_by(YEAR, MONTH, GENDER) %>%
    summarise(.,
        average_usage = mean(ITEMS/population)
    ) %>%
    mutate(., 
        month_year = paste0(YEAR, "/", MONTH, "/1"),
        monthdate = as.Date(month_year, format = "%Y/%m/%d"),
        ab = "nitro"
    )

trimnitro <- bind_rows(trim, nitro) 

png("TrimVsNitroPrescriptions.png", width = 25, height = 8, units = "cm", res = 300)
p <- ggplot(trimnitro, aes(monthdate, average_usage, col = ab)) + 
    geom_point() + 
    theme_bw() + 
    labs(x = "Month", y = "Usage per capita", col = "antibiotic") + 
    facet_wrap(~GENDER)
print(p)
dev.off()

#convert for working with the model
#need a left join not a bindrow
nitrosub <- nitro %>%
    ungroup %>%
    subset(GENDER == "Female") %>%
    rename(nitro = average_usage) %>%
    select(monthdate, nitro)

trimnitro2 <- trim %>%
    ungroup() %>%
    subset(GENDER == "Female") %>%
    rename(trim = average_usage) %>%
    select(monthdate, trim) %>%
    left_join(., nitrosub, by = "monthdate") 

#needs to be in days from 0
trimnitro2$timeindays <- as.numeric(abs(trimnitro2$monthdate[1] - trimnitro2$monthdate))
#needs to have relative usage rates
trimnitro2$relpres <- trimnitro2$trim/(trimnitro2$trim+trimnitro2$nitro)

trimsub <- trimnitro2 %>%
    select(timeindays, relpres)

write_csv(trimsub, "relativeprescriptionchanges.csv")

