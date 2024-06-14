## SLF Cold tolerance and lower lethal limits analysis
## Figure 1 and 
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
field_temperature_file <- "data/1_SLF_PA_GLFC_temp_data.csv" ## Fluctuating temperatures experienced by egg masses in Pennsylvania, USA, prior to collection, followed by laboratory and treatment conditions experienced in the Insect Production and Quarantine Laboratories from September 2022 until May 2023 (YR2).
field_temperature_sup_file <-"data/2_SLF_PA_temp_data_sup.csv" ##  Fluctuating temperatures experienced by egg masses in Pennsylvania, USA https://www.weather.gov/wrh/Climate?wfo=phi

################################################################################
# Field acclimation data - Fig 1
################################################################################

## read data -------------------------------------------------------------------
temperature_df <- read.csv(paste0(field_temperature_file), header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE)


## reshape data ----------------------------------------------------------------
SLF_temp <- gather(temperature_df, "treatment", "value", 3:48)
temp_cat <- as.data.frame(str_split_i(SLF_temp$treatment, "_", 3)) 
colnames(temp_cat) <- paste0("temp_cat")

SLF_temp <- bind_cols(SLF_temp,temp_cat)

Temp_control <- filter(SLF_temp, treatment %in% c("C_Temp_max_C","C_Temp_min_C"))

min <- Temp_control %>% dplyr::filter(temp_cat == "min") %>% dplyr::select(Day, value)
max <- Temp_control %>% dplyr::filter(temp_cat == "max") %>% dplyr::select(Day, value)
colnames(min) <- c(paste0("Day"),paste0("lower"))
colnames(max) <- c(paste0("Day"),paste0("upper"))

Temp_control <- left_join(Temp_control, min)
Temp_control <- left_join(Temp_control, max)


## Plot Fig. 1 -----------------------------------------------------------------
## Temperature exposure of Lycorma delicatula egg masses used in Y2 experiments - bars for days area days labels

Temp_control %>%
  ggplot( aes(x=Day, y=value, group = treatment , color = temp_cat)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha = 0.3, fill = "#3670CE")+
  annotate("rect", xmin = 113, xmax = 140, ymin = -30, ymax = 30,
           alpha = .5,fill = "#fee6ce")+
  annotate("rect", xmin = 140, xmax = 210, ymin = -30, ymax = 30,
           alpha = .5,fill = "#fff7ec")+
  geom_line() +
  scale_color_manual(values = c("#0753A5", "#0753A5")) +
  scale_y_continuous(limits = c(-30, 30)) +
  geom_vline(xintercept = 31, color = "#66c2a5", linetype="dotted") + # Oct
  geom_vline(xintercept = 62, color = "#fc8d62", linetype="dotted") + # Nov
  geom_vline(xintercept = 92, color = "#8da0cb", linetype="dotted") + # Dec
  geom_vline(xintercept = 123, color = "#e78ac3", linetype="dotted") + # Jan
  geom_vline(xintercept = 154, color = "#a6d854", linetype="dotted") + # Feb
  geom_vline(xintercept = 182, color = "#ffd92f", linetype="dotted") + # Mar
  geom_vline(xintercept = 213, color = "#e5c494", linetype="dotted") + # Apr
  geom_vline(xintercept = 243, color = "#b3b3b3", linetype="dotted") + # May
  geom_vline(xintercept = 274, color = "#a65628", linetype="dotted") + # Jun
  
  geom_vline(xintercept = 113, color = "#f16913", linetype="dotted") + # day collected from the field
  geom_vline(xintercept = 140, color = "#f16913", linetype="dotted") + # day received in SSM
  geom_vline(xintercept = 146, color = "#f16913", linetype="dotted") + # day of 1h experiments
  geom_vline(xintercept = 147, color = "#f16999", linetype="dotted") + # day of 12 hr experiments
  geom_vline(xintercept = 148, color = "#91003f", linetype="dotted") + # day of other experiments
  geom_vline(xintercept = 159, color = "#ce1256", linetype="dotted") + # day of 12 hr experiments
  geom_vline(xintercept = 166, color = "#e7298a", linetype="dotted") + # day of 12 hr experiments
  geom_vline(xintercept = 212, color = "#df65b0", linetype="dotted") + # day of 12 hr experiments
  geom_vline(xintercept = 218, color = "#c994c7", linetype="dotted") + # day of 12 hr experiments
  
  geom_line(y = -25, color = "#045a8d", linetype="dotted") +
  geom_line(y = -22.5, color = "#045a8d", linetype="dotted") +
  geom_line(y = -20, color = "#045a8d", linetype="dotted") +
  geom_line(y = -15, color = "#045a8d", linetype="dotted") +
  geom_line(y = -10, color = "#045a8d", linetype="dotted") +
  geom_line(y = -5, color = "#045a8d", linetype="dotted") +
  geom_line(y = 0, color = "#045a8d", linetype="dotted") +
  
  annotate("rect", xmin = 146, xmax = 147, ymin = -25, ymax = 6,
           alpha = .5,fill = "#045a8d")+
  
  annotate("rect", xmin = 149, xmax = 163, ymin = -20, ymax = 6,
           alpha = .5,fill = "#045a8d")+
  
  annotate("rect", xmin = 159, xmax = 169, ymin = -20, ymax = 6,
           alpha = .5,fill = "#045a8d")+
  annotate("rect", xmin = 166, xmax = 177, ymin = -20, ymax = 6,
           alpha = .5,fill = "#045a8d")+
  theme_ipsum() +
  # theme_void() +
  # theme(legend.position="none",axis.line=element_blank()) +
  theme(legend.position="none") +
  ylab(expression(paste(
    "Minimum and maximum temperature (" *degree* "C)"
  ))) +
  xlab("Day")+
  NULL






################################################################################
# Field temperature - Fig S1
################################################################################

## Clear environment except from files needed for analysis ---------------------
rm(list = ls()[!ls() %in% c("field_temperature_file","field_temperature_sup_file")])


## Read data -------------------------------------------------------------------
temperature_df <- read.csv(paste0(field_temperature_sup_file), header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE)

## Reshape data ----------------------------------------------------------------
d21 <- temperature_df %>% dplyr::select(Day,Temp.max.2021,Temp.min.2021)
d21 <- gather(d21, "Category", "value", 2:3)

d22 <- temperature_df %>% dplyr::select(Day,Temp.max.2022,Temp.min.2022)
d22 <- gather(d22, "Category", "value", 2:3)

d2122 <- bind_rows(d21,d22)

temp_2122 <- left_join(temperature_df, d2122)


## Plot Fig. S1 ----------------------------------------------------------------

cairo_pdf(file = "output/FigS1_TemperatureProfile2.pdf", width = 8, height = 6)

temp_2122 %>%
  ggplot( aes(x=Day, y=value, group = Category, color = Category)) +
  geom_ribbon(aes(ymin=Temp.min.2021, ymax=Temp.max.2021), alpha = 0.4, fill = "#5fa2ce")+
  geom_ribbon(aes(ymin=Temp.min.2022, ymax=Temp.max.2022), alpha = 0.2, fill = "#ffbc79")+
  # annotate("rect", xmin = 113, xmax = 140, ymin = -30, ymax = 30,
  #          alpha = .5,fill = "#fee6ce")+
  # annotate("rect", xmin = 140, xmax = 210, ymin = -30, ymax = 30,
  #          alpha = .5,fill = "#fff7ec")+
  geom_line() +
  scale_color_manual(values = c("#1170aa",  "#fc7d0b","#1170aa", "#fc7d0b")) +
  scale_y_continuous(limits = c(-10, 35)) +
  geom_vline(xintercept = 31, color = "#66c2a5", linetype="dotted") + # Oct
  geom_vline(xintercept = 62, color = "#fc8d62", linetype="dotted") + # Nov
  
  geom_line(y = -5, color = "#045a8d", linetype="dotted") +
  geom_line(y = 0, color = "#045a8d", linetype="dotted") +
  
  theme_ipsum() +
  
  theme(legend.position="none") +
  ylab(expression(paste(
    "Minimum and maximum temperature (" *degree* "C)"
  ))) +
  xlab("Day")+
  NULL

dev.off()


## Analysis: differences in temp between months year 1 and 2 -------------------
df <- gather(temperature_df, "Category", "value", 4:7)
df$Category2 = paste0(df$Month, "_", df$Category)

k_test <- kruskal.test(value ~ Category2, data = df) 

wilco <- pairwise.wilcox.test(df$value, df$Category2, p.adjust.method = "BH")

wilco_results <- wilco$p.value

write.csv(wilco_results, "output/wilco_results_temp2021_22.csv")  




df$catMinMax = sub("Temp.", "", df$Category)
df$catMinMax = sub(".2021", "", df$catMinMax)
df$catMinMax = sub(".2022", "", df$catMinMax)

df$Year = sub("Temp.min.", "", df$Category)
df$Year = sub("Temp.max.", "", df$Year)

sept <- df %>% filter(Month == "Sep")
oct <- df %>% filter(Month == "Oct")
nov <- df %>% filter(Month == "Nov")

sept_min <- sept %>% filter(catMinMax == "min")
sept_max <- sept %>% filter(catMinMax == "max")

oct_min <- oct %>% filter(catMinMax == "min")
oct_max <- oct %>% filter(catMinMax == "max")

nov_min <- nov %>% filter(catMinMax == "min")
nov_max <- nov %>% filter(catMinMax == "max")


k_test_s_min <- kruskal.test(value ~ Year, data = sept_min) 
k_test_s_max <- kruskal.test(value ~ Year, data = sept_max)

k_test_o_min <- kruskal.test(value ~ Year, data = oct_min) 
k_test_o_max <- kruskal.test(value ~ Year, data = oct_max) 

k_test_n_min <- kruskal.test(value ~ Year, data = nov_min) 
k_test_n_max <- kruskal.test(value ~ Year, data = nov_max) 

temperature_df %>%
  ggplot( aes(x=Temp.max.2021, y=Temp.max.2022, group = Month, color = Month)) +
  geom_point()+
  geom_abline(intercept = 1)




## reshape data ----------------------------------------------------------------
dmin <- temperature_df %>% dplyr::select(Day,Temp.min.2021, Temp.min.2022)
dmin <- gather(dmin, "Category", "value", 2:3)

dmin$minMax <- "MinimumTemperature"

dmax <- temperature_df %>% dplyr::select(Day,Temp.max.2021,Temp.max.2022)
dmax <- gather(dmax, "Category", "value", 2:3)

dmax$minMax <- "MaximumTemperature"

dminmax <- bind_rows(dmin,dmin)

temp_2122 <- left_join(temperature_df, d2122)


data <- data.frame(var1 = c(0, 0, 1, 1), var2 = c(0, 1, 0, 1))
model <- glm(var2 ~ var1, data = data, family = binomial)
summary(model)








