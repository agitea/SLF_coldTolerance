## SLF Cold tolerance and lower lethal limits analysis
## Temperature / time treatment hatch rate analysis - Fig.3
## written by Anna J. Turbelin, Mar 14, 2024 
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
library(jtools)
library(effects)
library(ggpmisc)
library(ggeffects)

theme_set(theme_bw())

options(stringsAsFactors=FALSE)



## Data file -------------------------------------------------------------------
YR2_hatch_file <- "data/SLF_YR1_YR2_hatch2.csv" ## cold tolerance strategy data 


## read data -------------------------------------------------------------------
df <- read.csv(paste0(YR2_hatch_file), header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE)


## reshape data ----------------------------------------------------------------
df <- df %>% dplyr::filter(TotalEggs != "#N/A")
df <- df %>% dplyr::filter(PercentageHatch != "#DIV/0!")
df$TreatmentTemp_c <- as.character(df$TreatmentTemp)
df$TreatmentTemp <- as.numeric(df$TreatmentTemp)
df$PercentageHatch <- as.numeric(df$PercentageHatch)

df_yr2 <- df %>% dplyr::filter(YR == "2")

df.sub <- df_yr2%>%
  group_by(TreatmentTime, TreatmentTemp) %>% 
  mutate(MaxHatch = max(PercentageHatch))

df.sub <- df.sub%>%
  group_by(TreatmentTime, TreatmentTemp) %>% 
  mutate(MeanHatch = mean(PercentageHatch))

################################################################################
## Mean percentage (a) and maximum (b) hatch of Lycorma delicatula egg masses exposed to 
## different treatments (0°C, -5°C, -10°C, -15°C, -20°C, -22.5°C, -25°C) and durations (1 hour, 12 hours, 10 days, 15 days) 
## during January and February 2023.
## Fig.3 a ----------------------------------------------------------------------
ggplot(df.sub, aes(reorder(TreatmentTemp_c, df.sub$TreatmentTemp), TreatmentTime, fill= MeanHatch)) + 
  scale_fill_viridis(option = "G" ,discrete=FALSE, direction = -1, begin = 0.4) +
  geom_tile()+
  theme(legend.position = "bottom")


## Fig.3 b
ggplot(df.sub, aes(reorder(TreatmentTemp_c, df.sub$TreatmentTemp), TreatmentTime, fill= MaxHatch)) + 
  scale_fill_viridis(option = "G" ,discrete=FALSE, direction = -1) +
  geom_tile()+
  theme(legend.position = "bottom")


## basic stats - Table S1 ------------------------------------------------------
std <- function(x) sd(x)/sqrt(length(x))

df_new <- df %>% dplyr::group_by(YR,TreatmentTemp, TreatmentTime) %>% mutate(meanhatchrate = mean(PercentageHatch)) 
df_new2 <- df %>% dplyr::group_by(YR,TreatmentTemp, TreatmentTime) %>% mutate(medianhatchrate = median(PercentageHatch))
df_new3 <- df %>% dplyr::group_by(YR,TreatmentTemp, TreatmentTime) %>% mutate(std_error = std(PercentageHatch))

df_all <- left_join(df_new, df_new2)
df_all <- left_join(df_all, df_new3)
df_all$YR <- as.numeric(df_all$YR)

df_sub <- dplyr::select(df_all, YR,TreatmentTemp, TreatmentTime, meanhatchrate, medianhatchrate, std_error)
df_sub <- df_all %>% dplyr::group_by(YR,TreatmentTemp, TreatmentTime, meanhatchrate, medianhatchrate, std_error) %>% summarise(count = n())

df_eggcount <- df %>% dplyr::group_by(YR,TreatmentTemp, TreatmentTime) %>% count(TotalHatch == 0)
names(df_eggcount)[names(df_eggcount) == 'TotalHatch == 0'] <- 'HatchCount'
df_eggcount$HatchCount <- as.character(df_eggcount$HatchCount)

df_eggcount <- df_eggcount %>% tidyr::spread(key = HatchCount, value= n)
names(df_eggcount)[names(df_eggcount) == 'FALSE'] <- 'HatchCount'
names(df_eggcount)[names(df_eggcount) == 'TRUE'] <- 'NoHatchCount'

TableS1 <- left_join(df_eggcount, df_sub)

TableS1$PercentageHatch <- TableS1$HatchCount / TableS1$count

TableS1 <- TableS1[,c(1:3,9,4,10,6:8)]

# write.csv(df_sub, "output/SLF_summary_hatchRate_stats_TableS1.csv" ) # Supplementary material Data 



## Kruskal-Wallis rank sum test comparing treatments to control in YR2 ---------

df_yr2$Treatment <- paste(df_yr2$TreatmentTemp, df_yr2$TreatmentTime)

list <- unique(df_yr2$Treatment)

output <- data.frame(matrix(nrow = 0, ncol= 6), stringsAsFactors = FALSE) ## create blank dataframe
colnames(output) <- c("Treatment","TreatmentTemp", "TreatmentTime", "kruskal_pvalue", "chiSquared", "KW_df") ## name column of blank df

for (i in 2:length(list)){
  t <- list[i]
  control <- dplyr::filter(df_yr2, Treatment == "6 control")
  treat_data <- dplyr::filter(df_yr2, Treatment == t)
  data <- bind_rows(control,treat_data)
  k_test <- kruskal.test(PercentageHatch ~ Treatment, data = data)
  temp <- unique(treat_data$TreatmentTemp)
  time <- unique(treat_data$TreatmentTime)
  output[i,which(colnames(output)== "Treatment")] <- t
  output[i,which(colnames(output)== "TreatmentTemp")] <- temp
  output[i,which(colnames(output)== "TreatmentTime")] <- time
  output[i,which(colnames(output)== "kruskal_pvalue")] <- k_test$p.value
  output[i,which(colnames(output)== "chiSquared")] <- k_test$statistic
  output[i,which(colnames(output)== "KW_df")] <- k_test$parameter
  
}

# write.csv(output, "output/SLF_stats_summary_treatment_control_TableS3.csv")  ## Supplementary material Table S3

################################################################################
## compare survival of temperature/time treatments to no exposure survival - contingency table

rm(list = ls()[!ls() %in% c("YR2_hatch_file", "bi.LLT.data", "bi.LLT.data.theo")])

## read data -------------------------------------------------------------------
hatch_df <- read.csv(paste0(YR2_hatch_file), header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE)

## reshape data ----------------------------------------------------------------
hatch_df<- hatch_df %>% dplyr::filter(TotalEggs != "#N/A")
hatch_df$TreatmentTime[hatch_df$TreatmentTime == 'control'] <- (66*24)+3000
hatch_df$TreatmentTemp[hatch_df$TreatmentTemp == 'control'] <- 6
hatch_df <- hatch_df %>% dplyr::filter(PercentageHatch != "#DIV/0!")
hatch_df$TotalHatch <- as.numeric(hatch_df$TotalHatch)
hatch_df$TotalEggs <- as.numeric(hatch_df$TotalEggs)
hatch_df<- hatch_df %>% dplyr::filter(!is.na(TotalEggs))

db.sub <- hatch_df
db.sub$total_alive <- as.numeric(db.sub$TotalHatch)
db.sub$TotalEggs <- as.numeric(db.sub$TotalEggs)
db.sub$total_dead <- db.sub$TotalEggs - db.sub$total_alive
db.sub <- dplyr::filter(db.sub, !is.na(db.sub$total_dead))
db.sub$YR <- as.character(db.sub$YR)
bi_db <- db.sub
bi_db <- bi_db %>% dplyr::filter(YR == "2")

bi_db$TreatmentTime <- as.character(bi_db$TreatmentTime)

t.time <- unique(bi_db$TreatmentTime)

for (j in 1:length(t.time)){
  
  control_ti <- "4584"
  ti = t.time[j]
  sub_time <- filter(bi_db,TreatmentTime%in% c(control_ti, ti))
  list <- unique(sub_time$TreatmentTemp)
  
  for (i in 1:length(list)){
    
    control_te <- "6"
    te <- list[i]
    
    sub_temp <- filter(sub_time,TreatmentTemp %in% c(control_te,te))
    
    sub_data <-  sub_temp %>% group_by(TreatmentTime,TreatmentTemp)%>% summarise(total_alive = sum(total_alive)) 
    sub_data2 <-  sub_temp %>% group_by(TreatmentTime,TreatmentTemp)%>% summarise(total_dead = sum(total_dead)) 
    
    sub_data <- left_join(sub_data, sub_data2)
    data <- sub_data[,3:4]
    
    df <- data.frame(matrix(nrow = 1, ncol= 4), stringsAsFactors = FALSE) ## create blank dataframe
    colnames(df) <- c("TreatmentTime", "TreatmentTemp","pvalue", "Xsquared") ## name column
    
    Xsq <- chisq.test(data)
    
    df$TreatmentTime <- ti
    df$TreatmentTemp <- te
    df$pvalue <- chisq.test(data)$p.value 
    df$Xsquared <- chisq.test(data)$statistic
    
    
    # bind fitness datasets of all temperatures at time j
    if (exists("df_results") == TRUE){
      df_results <- bind_rows(df_results, df)
    } else {df_results <- df}
    
    
    
  }
  
}

# write.csv(df_results, "output/SLF_ChiTest_results_contolComparison.csv")  ## Supplementary material Table S2
