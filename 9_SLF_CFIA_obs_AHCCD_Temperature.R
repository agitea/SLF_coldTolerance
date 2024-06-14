## SLF Cold tolerance and lower lethal limits analysis
## 
## written by Anna J. Turbelin, July 13, 2023 
## aturbelin@gmail.com

library(tidyr)
library(ggplot2)
library(ggpubr)
library(Cairo)
library(readxl)
library(tidyverse)
library(hrbrthemes)
library(RColorBrewer)
library(cshapes)
library(magrittr)
library(plyr)
detach(package:plyr)    
library(dplyr)
library(sp)
library(cowplot)
library(MASS)
install.packages("viridis")
library(viridis)
library(scales)
library(ggthemes)

library(Matrix)
library(lme4)
library(emmeans)

theme_set(theme_bw())

options(stringsAsFactors=FALSE)

rm(list=ls())

## data files ------------------------------------------------------------------
temperature_file <- "data/6_SLF_obs_climateData_AHCCD.csv" ## Climate data from https://climatedata.ca/download/#ahccd-download

## read data -------------------------------------------------------------------
temperature_df <- read.csv(paste0(temperature_file), header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE)

## filter data -----------------------------------------------------------------
temperature_df <- temperature_df %>% filter(year > 1972)

temperature_df <- temperature_df %>% filter(month %in% c(10, 11, 12, 01, 02, 03, 04, 05))

NA_temp <- temperature_df %>% filter(is.na(tas))

temperature_df <- temperature_df %>% filter(!is.na(tasmin))
temperature_df <- temperature_df %>% filter(!is.na(tasmax))

## add winter classification ---------------------------------------------------

temperature_df$year <- as.numeric(temperature_df$year)
temperature_df$month <- as.numeric(temperature_df$month)

for (i in 1:nrow(temperature_df)) {
  
  yr <- temperature_df$year[i]
  yr0 <- yr - 1
  yr2 <- yr + 1
  
  if(temperature_df$month[i] < 6){
    
    temperature_df$winter[i] <- paste0("winter_", yr0 , "_", yr)
  
    } else {

    temperature_df$winter[i] <- paste0("winter_", yr, "_", yr2)
    
      }
  
}



################################################################################
## identify days below certain temperature function ----------------------------

countDaysBelowT <- function(df, temperature_col = NULL, temperature = NULL, group_col = NULL) {

  df <- df %>%
    group_by({{group_col}}) %>%
    mutate(cumul20 = ifelse({{temperature_col}} < temperature, cumsum({{temperature_col}} < temperature), 0),
           cumul20 = ifelse({{temperature_col}} < temperature, 1, cumul20)) %>%
    ungroup()
  names(df)[names(df) == 'cumul20'] <- paste0("days_", abs(temperature))
  
  return(df)
}

################################################################################

## questions to answer:
## 1. What is the maximum number of accumulated days where the minimum or maximum daily temperature is below -15 and -20 in a winter? 
## 2. How frequently do you get more than 10 or 15 cumulated days where the minimum or maximum daily temperature is below -15 or -20 in a winter? (e.g. zero, once, twice per winter)
## 3. How many days does the minimum daily temperature go below -20 or -27.7 in a winter? 
  
## count days below given temperature ------------------------------------------
df <- temperature_df

## -20 C daily min temperature -------------------------------------------------
temperature <- -20 # set temperature
temperature_col <- "tasmin" 

df <- countDaysBelowT(df, temperature_col = tasmin, temperature = temperature, group_col = winter)

names(df)[names(df) == 'days_20'] <- "days_20_min"

df$dayscount <- df$days_20_min # add new column


for (i in 1:nrow(df)) {
  if (i==1) {df$cumul_days = df$dayscount}
  if (i > 1) { 
    if(isTRUE(df$dayscount[i] == 0)==TRUE) {df$cumul_days[i] = df$dayscount[i]} else { df$cumul_days[i] = df$dayscount[i] + df$cumul_days[i-1]}
  }
}


df <- df %>%
  group_by(station_name, winter, cumul_days) %>%
  mutate(days_frequency = n()) %>%
  ungroup()


names(df)[names(df) == 'cumul_days'] <- paste0("cumul_days_",temperature_col, abs(temperature))
names(df)[names(df) == 'days_frequency'] <- paste0("days_frequency",temperature_col,abs(temperature))


summary <- df  %>% group_by(station_name, winter, cumul_days_tasmin20, days_frequencytasmin20) %>% summarise()

streakFrequency20_10 <- summary %>% filter(cumul_days_tasmin20 == 10)
names(streakFrequency20_10)[names(streakFrequency20_10) == 'cumul_days_tasmin20'] <- "cumul_days_tasmin20_10"
names(streakFrequency20_10)[names(streakFrequency20_10) == 'days_frequencytasmin20'] <- "days_frequencytasmin20_10"


streakFrequency20_15 <- summary %>% filter(cumul_days_tasmin20 == 15)
names(streakFrequency20_15)[names(streakFrequency20_15) == 'cumul_days_tasmin20'] <- "cumul_days_tasmin20_15"
names(streakFrequency20_15)[names(streakFrequency20_15) == 'days_frequencytasmin20'] <- "days_frequencytasmin20_15"


maxstreak20 <- df %>% group_by(station_name, winter) %>% summarise(max_cumul20_min = max(cumul_days_tasmin20))

daysBelow20 <- df  %>% group_by(station_name, winter) %>% summarise(daysbelow20_min_count = sum(dayscount)) # count days when the minimum daily temperature goes below -20

cityTempSummary <- left_join(daysBelow20, maxstreak20)
cityTempSummary <- left_join(cityTempSummary, streakFrequency20_10)
cityTempSummary <- left_join(cityTempSummary, streakFrequency20_15)

df_20_min <- df

## -20 C daily max temperature -------------------------------------------------
df <- temperature_df

temperature <- -20 # set temperature
temperature_col <- "tasmax" 

df <- countDaysBelowT(df, temperature_col = tasmax, temperature = temperature, group_col = winter)

names(df)[names(df) == 'days_20'] <- "days_20_max"
df$dayscount <- df$days_20_max # add new column

for (i in 1:nrow(df)) {
  if (i==1) {df$cumul_days = df$dayscount}
  if (i > 1) { 
    if(isTRUE(df$dayscount[i] == 0)==TRUE) {df$cumul_days[i] = df$dayscount[i]} else { df$cumul_days[i] = df$dayscount[i] + df$cumul_days[i-1]}
  }
}


df <- df %>%
  group_by(station_name, winter, cumul_days) %>%
  mutate(days_frequency = n()) %>%
  ungroup()


names(df)[names(df) == 'cumul_days'] <- paste0("cumul_days_",temperature_col, abs(temperature))
names(df)[names(df) == 'days_frequency'] <- paste0("days_frequency",temperature_col,abs(temperature))

head(df)

summary <- df  %>% group_by(station_name, winter, cumul_days_tasmax20, days_frequencytasmax20) %>% summarise()

streakFrequency20_10 <- summary %>% filter(cumul_days_tasmax20 == 10)
names(streakFrequency20_10)[names(streakFrequency20_10) == 'cumul_days_tasmax20'] <- "cumul_days_tasmax20_10"
names(streakFrequency20_10)[names(streakFrequency20_10) == 'days_frequencytasmax20'] <- "days_frequencytasmax20_10"


streakFrequency20_15 <- summary %>% filter(cumul_days_tasmax20 == 15)
names(streakFrequency20_15)[names(streakFrequency20_15) == 'cumul_days_tasmax20'] <- "cumul_days_tasmax20_15"
names(streakFrequency20_15)[names(streakFrequency20_15) == 'days_frequencytasmax20'] <- "days_frequencytasmax20_15"


maxstreak20 <- df %>% group_by(station_name, winter) %>% summarise(max_cumul20_max = max(cumul_days_tasmax20))

daysBelow20 <- df  %>% group_by(station_name, winter) %>% summarise(daysbelow20_count_max = sum(dayscount)) # count days when the minimum daily temperature goes below -20

cityTempSummary <- left_join(cityTempSummary, daysBelow20)
cityTempSummary <- left_join(cityTempSummary, maxstreak20)
cityTempSummary <- left_join(cityTempSummary, streakFrequency20_10)
cityTempSummary <- left_join(cityTempSummary, streakFrequency20_15)


df_20_max <- df

## -15 C daily min temperature -------------------------------------------------
df <- temperature_df

temperature <- -15 # set temperature
temperature_col <- "tasmin" 

df <- countDaysBelowT(df, temperature_col = tasmin, temperature = temperature, group_col = winter)

names(df)[names(df) == 'days_15'] <- "days_15_min"

df$dayscount <- df$days_15_min # add new column

for (i in 1:nrow(df)) {
  if (i==1) {df$cumul_days = df$dayscount}
  if (i > 1) { 
    if(isTRUE(df$dayscount[i] == 0)==TRUE) {df$cumul_days[i] = df$dayscount[i]} else { df$cumul_days[i] = df$dayscount[i] + df$cumul_days[i-1]}
  }
}


df <- df %>%
  group_by(station_name, winter, cumul_days) %>%
  mutate(days_frequency = n()) %>%
  ungroup()


names(df)[names(df) == 'cumul_days'] <- paste0("cumul_days_",temperature_col, abs(temperature))
names(df)[names(df) == 'days_frequency'] <- paste0("days_frequency",temperature_col,abs(temperature))


summary <- df  %>% group_by(station_name, winter, cumul_days_tasmin15, days_frequencytasmin15) %>% summarise()

streakFrequency15_10 <- summary %>% filter(cumul_days_tasmin15 == 10)
names(streakFrequency15_10)[names(streakFrequency15_10) == 'cumul_days_tasmin15'] <- "cumul_days_tasmin15_10"
names(streakFrequency15_10)[names(streakFrequency15_10) == 'days_frequencytasmin15'] <- "days_frequencytasmin15_10"


streakFrequency15_15 <- summary %>% filter(cumul_days_tasmin15 == 15)
names(streakFrequency15_15)[names(streakFrequency15_15) == 'cumul_days_tasmin15'] <- "cumul_days_tasmin15_15"
names(streakFrequency15_15)[names(streakFrequency15_15) == 'days_frequencytasmin15'] <- "days_frequencytasmin15_15"


maxstreak15 <- df %>% group_by(station_name, winter) %>% summarise(max_cumul15_min = max(cumul_days_tasmin15))

daysBelow15 <- df  %>% group_by(station_name, winter) %>% summarise(daysbelow15_min_count = sum(dayscount)) # count days when the minimum daily temperature goes below -20

cityTempSummary <- left_join(cityTempSummary, daysBelow15)
cityTempSummary <- left_join(cityTempSummary, maxstreak15)
cityTempSummary <- left_join(cityTempSummary, streakFrequency15_10)
cityTempSummary <- left_join(cityTempSummary, streakFrequency15_15)

df_15_min <- df

## -15 C daily max temperature -------------------------------------------------
df <- temperature_df

temperature <- -15 # set temperature
temperature_col <- "tasmax" 

df <- countDaysBelowT(df, temperature_col = tasmax, temperature = temperature, group_col = winter)

names(df)[names(df) == 'days_15'] <- "days_15_max"
df$dayscount <- df$days_15_max # add new column

for (i in 1:nrow(df)) {
  if (i==1) {df$cumul_days = df$dayscount}
  if (i > 1) { 
    if(isTRUE(df$dayscount[i] == 0)==TRUE) {df$cumul_days[i] = df$dayscount[i]} else { df$cumul_days[i] = df$dayscount[i] + df$cumul_days[i-1]}
  }
}


df <- df %>%
  group_by(station_name, winter, cumul_days) %>%
  mutate(days_frequency = n()) %>%
  ungroup()


names(df)[names(df) == 'cumul_days'] <- paste0("cumul_days_",temperature_col, abs(temperature))
names(df)[names(df) == 'days_frequency'] <- paste0("days_frequency",temperature_col,abs(temperature))

head(df)

summary <- df  %>% group_by(station_name, winter, cumul_days_tasmax15, days_frequencytasmax15) %>% summarise()

streakFrequency15_10 <- summary %>% filter(cumul_days_tasmax15 == 10)
names(streakFrequency15_10)[names(streakFrequency15_10) == 'cumul_days_tasmax15'] <- "cumul_days_tasmax15_10"
names(streakFrequency15_10)[names(streakFrequency15_10) == 'days_frequencytasmax15'] <- "days_frequencytasmax15_10"


streakFrequency15_15 <- summary %>% filter(cumul_days_tasmax15 == 15)
names(streakFrequency15_15)[names(streakFrequency15_15) == 'cumul_days_tasmax15'] <- "cumul_days_tasmax15_15"
names(streakFrequency15_15)[names(streakFrequency15_15) == 'days_frequencytasmax15'] <- "days_frequencytasmax15_15"


maxstreak15 <- df %>% group_by(station_name, winter) %>% summarise(max_cumul15_max = max(cumul_days_tasmax15))

daysBelow15 <- df  %>% group_by(station_name, winter) %>% summarise(daysbelow15_count_max = sum(dayscount)) # count days when the minimum daily temperature goes below -15

cityTempSummary <- left_join(cityTempSummary, daysBelow15)
cityTempSummary <- left_join(cityTempSummary, maxstreak15)
cityTempSummary <- left_join(cityTempSummary, streakFrequency15_10)
cityTempSummary <- left_join(cityTempSummary, streakFrequency15_15)

df_15_max <- df




## -27.7 C daily min temperature -------------------------------------------------
df <- temperature_df

temperature <- -27.7 # set temperature
temperature_col <- "tasmin" 

df <- countDaysBelowT(df, temperature_col = tasmin, temperature = temperature, group_col = winter)

df$dayscount <- df$days_27.7 # add new column

daysBelow27min <- df  %>% group_by(station_name, winter) %>% summarise(daysbelow27_count_min = sum(days_27.7)) # count days when the minimum daily temperature goes below -27.7

cityTempSummary <- left_join(cityTempSummary, daysBelow27min)


## -27.7 C daily max temperature -------------------------------------------------
df <- temperature_df

temperature <- -27.7 # set temperature
temperature_col <- "tasmax" 

df <- countDaysBelowT(df, temperature_col = tasmax, temperature = temperature, group_col = winter)

df$dayscount <- df$days_27.7 # add new column

daysBelow27max <- df  %>% group_by(station_name, winter) %>% summarise(daysbelow27_count_max = sum(days_27.7)) # count days when the minimum daily temperature goes below -27.7


cityTempSummary <- left_join(cityTempSummary, daysBelow27max)


## -25 C daily min temperature -------------------------------------------------
df <- temperature_df

temperature <- -25 # set temperature
temperature_col <- "tasmin" 

df <- countDaysBelowT(df, temperature_col = tasmin, temperature = temperature, group_col = winter)

df$dayscount <- df$days_25 # add new column

daysBelow25min <- df  %>% group_by(station_name, winter) %>% summarise(daysbelow25_count_min = sum(days_25)) # count days when the minimum daily temperature goes below -25

cityTempSummary <- left_join(cityTempSummary, daysBelow25min)

## -25 C daily max temperature -------------------------------------------------
df <- temperature_df

temperature <- -25 # set temperature
temperature_col <- "tasmax" 

df <- countDaysBelowT(df, temperature_col = tasmax, temperature = temperature, group_col = winter)

df$dayscount <- df$days_25 # add new column

daysBelow25max <- df  %>% group_by(station_name, winter) %>% summarise(daysbelow25_count_max = sum(days_25)) # count days when the minimum daily temperature goes below -25


cityTempSummary <- left_join(cityTempSummary, daysBelow25max)


## join datasets--- ------------------------------------------------------------
df_all <- left_join(df_20_min, df_20_max)
df_all<- left_join(df_all, df_15_min)
df_all <- left_join(df_all, df_15_max)


## summary 50 years ------------------------------------------------------------

colnames(cityTempSummary)

everyX <- function(x) 50/x

## get the average number of days below given temperature and the average maximum number of accumulated days below -20 and -15
daysBelowCount <- cityTempSummary %>% dplyr::select(station_name, daysbelow20_min_count,
                                                    daysbelow25_count_min,  daysbelow27_count_min,
                                                    max_cumul15_max, max_cumul20_max)


daysBelowCount_mean50 <- daysBelowCount %>% group_by(station_name) %>% summarise(across(everything(), list(mean)))

daysBelowCount <- cityTempSummary %>% dplyr::select(station_name, daysbelow20_min_count,
                                                    daysbelow25_count_min,  daysbelow27_count_min)

daysBelowCount_sum50 <- daysBelowCount %>% group_by(station_name) %>% summarise(across(everything(), list(sum)))

daysBelowCount_sum50$daysbtw20_25min <- daysBelowCount_sum50$daysbelow20_min_count_1 - daysBelowCount_sum50$daysbelow25_count_min_1
# daysBelowCount_sum50$daysbtw20_25max <- daysBelowCount_sum50$daysbelow20_count_max_1 - daysBelowCount_sum50$daysbelow25_count_max_1

daysBelowCount_sum50$everyxYear27 <-  everyX(daysBelowCount_sum50$daysbelow27_count_min_1)
daysBelowCount_sum50$everyxYear20_25 <-  everyX(daysBelowCount_sum50$daysbtw20_25min)

## In the last 50 years what was the maximum number of accumulated days @ -15 and -20

daysConsecutiveMax50maxtemp20 <- cityTempSummary %>% group_by(station_name) %>% summarise(maxCumulDays20_maxtemp = max(max_cumul20_max))
daysConsecutiveMax50mintemp20 <- cityTempSummary %>% group_by(station_name) %>% summarise(maxCumulDays20_mintemp = max(max_cumul20_min))

daysConsecutiveMax50maxtemp15 <- cityTempSummary %>% group_by(station_name) %>% summarise(maxCumulDays15_maxtemp = max(max_cumul15_max))
daysConsecutiveMax50mintemp15 <- cityTempSummary %>% group_by(station_name) %>% summarise(maxCumulDays15_mintemp = max(max_cumul15_min))


daysConsecutiveMax50 <- left_join(daysConsecutiveMax50maxtemp20 , daysConsecutiveMax50mintemp20)
daysConsecutiveMax50 <- left_join(daysConsecutiveMax50, daysConsecutiveMax50maxtemp15)
daysConsecutiveMax50 <- left_join(daysConsecutiveMax50, daysConsecutiveMax50mintemp15)

## get the total number of times when you get 10 or more and 15 or more consecutive days at -15 and -20
daysConsecutive50 <- cityTempSummary %>% dplyr::select(station_name, days_frequencytasmin15_10, days_frequencytasmax15_10,
                                                       days_frequencytasmin15_15, days_frequencytasmax15_15,
                                                       days_frequencytasmin20_10, days_frequencytasmax20_10,
                                                       days_frequencytasmin20_15, days_frequencytasmax20_15)


daysConsecutive50[is.na(daysConsecutive50)] <- 0

daysConsecutive50_sum <- daysConsecutive50 %>% group_by(station_name) %>% summarise(across(everything(), list(sum)))


daysConsecutive50_mean50 <- daysConsecutive50 %>% group_by(station_name) %>% summarise(across(everything(), list(mean)))


everyX <- function(x) 50/x

daysConsecutive50_everyxyear <- data.frame(daysConsecutive50_sum[1], lapply(daysConsecutive50_sum[2:9], everyX) )

## reshape dataset

cityTemp50 <- daysBelowCount_sum50 %>% dplyr::select(station_name, everyxYear27, everyxYear20_25)

cityTemp50a <- daysConsecutive50_everyxyear %>% dplyr::select(station_name, days_frequencytasmax15_10_1,days_frequencytasmax20_15_1)

cityTemp50 <- left_join(cityTemp50, cityTemp50a)



cityTemp50 <- tidyr::gather(cityTemp50, variable, freq_value, 2:5)

## add scale for figure
cityTemp50$Rank <- cityTemp50$freq_value

cityTemp50$Rank[cityTemp50$freq_value >= 10] <- 4
cityTemp50$Rank[cityTemp50$freq_value < 10 & cityTemp50$freq_value >= 5] <- 3
cityTemp50$Rank[cityTemp50$freq_value < 5 & cityTemp50$freq_value >= 2] <- 2
cityTemp50$Rank[cityTemp50$freq_value < 2 & cityTemp50$freq_value >= 1] <- 1

cityTemp50$Rank[cityTemp50$freq_value == Inf] <- 0

cityTemp50$Rank[cityTemp50$freq_value < 1 & cityTemp50$freq_value >= 0.4] <- -1
cityTemp50$Rank[cityTemp50$freq_value < 0.4 & cityTemp50$freq_value >= 0.2] <- -2
cityTemp50$Rank[cityTemp50$freq_value < 0.2 & cityTemp50$freq_value >= 0.1] <- -3
cityTemp50$Rank[cityTemp50$freq_value < 0.1 ] <- -4


## define categorical scale

cityTemp50$freqCat <- cityTemp50$freq_value

cityTemp50$freqCat[cityTemp50$freq_value >= 10] <- "2_Every 10 years +"
cityTemp50$freqCat[cityTemp50$freq_value < 10 & cityTemp50$freq_value >= 5] <- "3_Every 5 - 9 years"
cityTemp50$freqCat[cityTemp50$freq_value < 5 & cityTemp50$freq_value >= 1.5] <- "4_Every 2 - 4 years"
cityTemp50$freqCat[cityTemp50$freq_value < 1.5 & cityTemp50$freq_value >= 0.67] <- "5_Every year"

cityTemp50$freqCat[cityTemp50$freq_value == Inf] <- "1_0 in 50 years"

cityTemp50$freqCat[cityTemp50$freq_value < 0.67 & cityTemp50$freq_value >= 0.23] <- "6_2 to 4 times per year"
cityTemp50$freqCat[cityTemp50$freq_value < 0.23 & cityTemp50$freq_value >= 0.1] <- "7_5 to 10 times per year"
# cityTemp50$freqCat[cityTemp50$freq_value < 0.2 & cityTemp50$freq_value >= 0.1] <- "2 to 4 times per year"
cityTemp50$freqCat[cityTemp50$freq_value < 0.1 ] <- "8_>10 times per year"


##

cityTemp50$station_name_r[cityTemp50$station_name == "WINNIPEG"] <- "1_WINNIPEG"
cityTemp50$station_name_r[cityTemp50$station_name == "WELLAND"] <- "2_WELLAND"
cityTemp50$station_name_r[cityTemp50$station_name == "OTTAWA"] <- "3_OTTAWA"
cityTemp50$station_name_r[cityTemp50$station_name == "ST HUBERT MONT"] <- "4_ST HUBERT MONT"
cityTemp50$station_name_r[cityTemp50$station_name == "QUEBEC"] <- "5_QUEBEC"
cityTemp50$station_name_r[cityTemp50$station_name == "MONCTON"] <- "6_MONCTON"
cityTemp50$station_name_r[cityTemp50$station_name == "HALIFAX"] <- "7_HALIFAX"
cityTemp50$station_name_r[cityTemp50$station_name == "SYDNEY"] <- "8_SYDNEY"


##

cityTemp50$pch_class <- cityTemp50$Rank

cityTemp50$pch_class[cityTemp50$freq_value >= 10] <- "EveryXYR"
cityTemp50$pch_class[cityTemp50$freq_value < 10 & cityTemp50$freq_value >= 5] <- "EveryXYR"
cityTemp50$pch_class[cityTemp50$freq_value < 5 & cityTemp50$freq_value >= 1.5] <- "EveryXYR"
cityTemp50$pch_class[cityTemp50$freq_value < 1.5 & cityTemp50$freq_value >= 0.67] <- "PerYear"

cityTemp50$pch_class[cityTemp50$freq_value == Inf] <- "EveryXYR"

cityTemp50$pch_class[cityTemp50$freq_value < 0.67 & cityTemp50$freq_value >= 0.23] <- "PerYear"
cityTemp50$pch_class[cityTemp50$freq_value < 0.23 & cityTemp50$freq_value >= 0.1] <- "PerYear"
# cityTemp50$pch_class[cityTemp50$freq_value < 0.2 & cityTemp50$freq_value >= 0.1] <- "2 to 4 times per year"
cityTemp50$pch_class[cityTemp50$freq_value < 0.1 ] <- "PerYear"

cityTemp50$variable_r[cityTemp50$variable == "everyxYear27"] <- "2_everyxYear27"
cityTemp50$variable_r[cityTemp50$variable == "everyxYear20_25"] <- "3_everyxYear20_25"
cityTemp50$variable_r[cityTemp50$variable == "days_frequencytasmax15_10_1"] <- "1_days_frequencytasmax15_10_1"
cityTemp50$variable_r[cityTemp50$variable == "days_frequencytasmax20_15_1"] <- "4_days_frequencytasmax20_15_1"
cityTemp50_sub <- dplyr::filter(cityTemp50, variable %in% c("everyxYear27", "days_frequencytasmax15_10_1"))


cityTemp50$Rank2 <- 1/cityTemp50$freq_value



cairo_pdf(file = "output/SLF_cityTemp_summary_square_size2.pdf", width = 10, height = 8)

cityTemp50 %>%
  ggplot(aes(x=station_name_r, y=variable_r, color = freqCat, fill = freqCat, pch = pch_class, size = Rank2)) +
  geom_point(alpha=1) +
  scale_shape_manual(values = c(22,21))+
  scale_size(range = c(4, 20), name="Event occurance", breaks = c(0,0.05, 0.2,0.4, 1, 2, 4,12))+
  scale_color_manual(values = c("#dfdfdf","#f9b97a","#e69f65","#ad8f79" , "#7e8288","#567796", "#275b8b","#001d39"))+ 
  scale_fill_manual(values = c("#F8F8F8","#fecd90","#f3a360","#ba9374" , "#8b8684","#677c90", "#3f719e","#00366c"))
NULL

dev.off()




## make map of observations from the CFIA

library("rnaturalearth")
library("rnaturalearthdata")
library("sf")
library("maps")
library("ggspatial")

station_obs_location <- "data/SLF_obs_station_latlon.csv" ## station lat lon aquired from https://climatedata.ca/download/#ahccd-download
## read data -------------------------------------------------------------------
station_obs_location_data  <- read.csv(paste0(station_obs_location), header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE)
latlon_obs <- filter(station_obs_location_data, Type_point == "Observation")
latlon_station <- filter(station_obs_location_data, Type_point == "Station")

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

lakes  <- read_sf("data/ne_50m_lakes/ne_50m_lakes.shp") # shapefile from Natural earth available online: https://www.naturalearthdata.com/downloads/50m-physical-vectors/


## Map CFIA SLF observations ---------------------------------------------------  
  
  cairo_pdf(file = "output/SLF_Map_obs_purple.pdf", width = 7, height = 10)
  
  ggplot(data = world) +
    geom_sf(fill = "#f4f3f1", color = "#252525") +
    geom_sf(data = lakes,fill = "#d1e5f0", color = "#d1e5f0", alpha = 0.6) +
    geom_point(data = latlon_station, aes(x = Lon, y = Lat), fill = "#7c71a1", color = "#443a67", pch = 23, size = 3, alpha = 0.7) +
    geom_point(data = latlon_obs, aes(x = Lon, y = Lat), fill = "#d6604d", color = "#67001f", pch = 25, size = 1.5) +
    coord_sf(xlim = c(-150, -50), ylim = c(40, 85), expand = FALSE)+
    annotation_scale()+
    annotate("rect", xmin = -100, xmax = -55, ymin = 40, ymax = 52,
             alpha = .08,fill = "blue", color = "#252525")
  
  NULL
  dev.off()
  
  cairo_pdf(file = "output/SLF_Map_obs_focus_purple.pdf", width = 7, height = 5)  
  ggplot(data = world) +
    geom_sf(fill = "#f4f3f1", color = "#252525") +
    geom_sf(data = lakes,fill = "#d1e5f0", color = "#d1e5f0", alpha = 0.6) +
    geom_point(data = latlon_station, aes(x = Lon, y = Lat), fill = "#7c71a1", color = "#443a67", pch = 23, size = 4, alpha = 0.7) +
    geom_point(data = latlon_obs, aes(x = Lon, y = Lat), fill = "#d6604d", color = "#67001f", pch = 25, size = 2) +
    coord_sf(xlim = c(-100, -55), ylim = c(40, 52), expand = FALSE)+
    annotation_scale()
    # theme_void()
    
    NULL
  
  dev.off()
  
  

