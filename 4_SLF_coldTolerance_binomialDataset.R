## SLF Cold tolerance and lower lethal limits analysis
## Create Binomial dataset
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
library(viridis)
library(scales)
library(ggthemes)
library(Matrix)
library(lme4)
library(emmeans)
library(jtools)
library(effects)


theme_set(theme_bw())

options(stringsAsFactors=FALSE)

## Raw data dataset ------------------------------------------------------------

## Data file -------------------------------------------------------------------
YR2_hatch_file <- "data/SLF_YR1_YR2_hatch2.csv" ## cold tolerance strategy data 

rm(list = ls()[!ls() %in% c("YR2_hatch_file")])

################################################################################
## Read data -------------------------------------------------------------------
hatch_df <- read.csv(paste0(YR2_hatch_file), header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE)

## Format data -----------------------------------------------------------------

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


## Create binomial dataset -----------------------------------------------------

list <- unique(db.sub$EggMass2)

for (i in 1:length(list)){
  E <- list[i]
  nr = as.numeric(db.sub[which(db.sub$EggMass2 == E, arr.ind=TRUE), which(colnames(db.sub)== "total_alive")])
  df0 <- data.frame(matrix(nrow = nr, ncol= 7), stringsAsFactors = FALSE) ## create blank dataframe
  colnames(df0) <- c("EggMass2","EggMass","EggMassClass","YR","TreatmentTime","TreatmentTemp","Fitness") ## name column
  df0$EggMass2 = as.character(df0$EggMass2)
  df0$EggMass = as.character(df0$EggMass)
  df0$EggMassClass = as.character(df0$EggMassClass)
  df0$TreatmentTime = as.numeric(df0$TreatmentTime)
  df0$TreatmentTemp = as.numeric(df0$TreatmentTemp)
  df0$YR = as.character(df0$YR)
  if(nr > 0) {
    df0$EggMass2 <- E
    df0$EggMass <- db.sub[which(db.sub$EggMass2 == E, arr.ind=TRUE),which(colnames(db.sub)== "EggMass")]
    df0$EggMassClass <- db.sub[which(db.sub$EggMass2 == E, arr.ind=TRUE),which(colnames(db.sub)== "EggMassClass")]
    df0$YR <- db.sub[which(db.sub$EggMass2 == E, arr.ind=TRUE),which(colnames(db.sub)== "YR")]
    df0$TreatmentTime <- as.numeric(db.sub[which(db.sub$EggMass2 == E, arr.ind=TRUE),which(colnames(db.sub)== "TreatmentTime")])
    df0$TreatmentTemp <- as.numeric(db.sub[which(db.sub$EggMass2 == E, arr.ind=TRUE),which(colnames(db.sub)== "TreatmentTemp")])
    df0$Fitness <- 1
  } 
  
  nrd = as.numeric(db.sub[which(db.sub$EggMass2 == E, arr.ind=TRUE),which(colnames(db.sub)=="total_dead" )])
  df1 <- data.frame(matrix(nrow = nrd, ncol= 7), stringsAsFactors = FALSE) ## create blank dataframe
  colnames(df1) <- c("EggMass2","EggMass","EggMassClass","YR","TreatmentTime","TreatmentTemp","Fitness") ## name column
  df1$EggMass2 = as.character(df1$EggMass2)
  df1$EggMass = as.character(df1$EggMass)
  df1$EggMassClass = as.character(df1$EggMassClass)
  df1$TreatmentTime = as.numeric(df1$TreatmentTime)
  df1$TreatmentTemp = as.numeric(df1$TreatmentTemp)
  df1$YR = as.character(df1$YR)
  if(nrd > 0) {
    df1$EggMass2 <- E
    df1$EggMass <- db.sub[which(db.sub$EggMass2 == E, arr.ind=TRUE),which(colnames(db.sub)== "EggMass")]
    df1$EggMassClass <- db.sub[which(db.sub$EggMass2 == E, arr.ind=TRUE),which(colnames(db.sub)== "EggMassClass")]
    df1$YR <- db.sub[which(db.sub$EggMass2 == E, arr.ind=TRUE),which(colnames(db.sub)== "YR")]
    df1$TreatmentTime <- as.numeric(db.sub[which(db.sub$EggMass2 == E, arr.ind=TRUE),which(colnames(db.sub)== "TreatmentTime")])
    df1$TreatmentTemp <- as.numeric(db.sub[which(db.sub$EggMass2 == E, arr.ind=TRUE),which(colnames(db.sub)== "TreatmentTemp")])
    df1$Fitness <- 0
  } 
  
  # bind fitness datasets for temperature i at time j
  df_temp <- bind_rows(df0, df1) 
  
  # bind fitness datasets of all temperatures at time j
  if (exists("df_temp_all") == TRUE){
    df_temp_all <- bind_rows(df_temp_all, df_temp)
  } else {df_temp_all = df_temp}
  
}

bi.LLT.data <- df_temp_all


check2 <- df_temp_all %>% group_by(YR,TreatmentTime,TreatmentTemp, EggMass,Fitness)%>% summarise(n = n()) 


###############################################################################
## Theoretical LT100 and LT0 data dataset --------------------------------------

rm(list = ls()[!ls() %in% c("YR2_hatch_file", "bi.LLT.data")])


## read data -------------------------------------------------------------------
hatch_df <- read.csv(paste0(YR2_hatch_file), header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE)


## prepare data ----------------------------------------------------------------

hatch_df<- hatch_df %>% dplyr::filter(TotalEggs != "#N/A")
hatch_df$TreatmentTime[hatch_df$TreatmentTime == 'control'] <- (66*24)+3000
hatch_df$TreatmentTemp[hatch_df$TreatmentTemp == 'control'] <- 6
hatch_df <- hatch_df %>% dplyr::filter(PercentageHatch != "#DIV/0!")
hatch_df$TotalHatch <- as.numeric(hatch_df$TotalHatch)
hatch_df$TotalEggs <- as.numeric(hatch_df$TotalEggs)
hatch_df <- hatch_df %>% dplyr::filter(!is.na(TotalEggs))

hatch_df_sub <- hatch_df[ , -which(names(hatch_df) %in% c("Source","TreatmentDate", "EggMass"))]
hatch_df_sub$total_dead <- hatch_df_sub$TotalEggs - hatch_df_sub$TotalHatch


## create survival proportion data frame with theoretical bounds ---------------
## LLT100 is 27.7 or the lowest supercooling point measurement
## LLT0 is 10.0 - tis is the required temp threshold for xxx (Keena et al. xxxx)

db.sub <- hatch_df_sub %>%
  group_by(YR, TreatmentTime, TreatmentTemp) %>%
  summarise(total = sum(TotalEggs), 
            total_alive = sum(TotalHatch),
            Prop_Survival = total_alive/total)

db.sub$total_dead <- db.sub$total - db.sub$total_alive


## LLT100 data
LLT100 <- hatch_df %>%
  group_by(YR, TreatmentTime) %>%
  summarise(meaneggs = mean(TotalEggs))

LLT100$TreatmentTemp = -27.7

LLT100$meaneggs = round(LLT100$meaneggs, digits = 0)

LLT100 <- LLT100 %>%
  mutate(total = ifelse(YR == 1, meaneggs * 7, meaneggs * 10))

LLT100 <- LLT100 %>% dplyr::filter(TreatmentTime != "4584")

head(LLT100)

LLT100$total_alive = 0
LLT100$Prop_Survival = LLT100$total_alive/LLT100$total
LLT100$total_dead <- LLT100$total - LLT100$total_alive
LLT100$EggMassClass <- "single"
LLT100$EggMass2 <- c(300,301,302,303,304)

## LLT0 data
LLT0 <- hatch_df %>%
  group_by(YR, TreatmentTime) %>%
  summarise(meaneggs = mean(TotalEggs))

LLT0$TreatmentTemp = 10
LLT0$meaneggs = round(LLT0$meaneggs, digits = 0)

LLT0 <- LLT0 %>%
  mutate(total = ifelse(YR == 1, meaneggs * 7, meaneggs * 10))

LLT0 <- LLT0 %>% dplyr::filter(TreatmentTime != "4584")

head(LLT0)

LLT0$total_alive = round(LLT0$total, digits = 0)
LLT0$Prop_Survival = LLT0$total_alive/LLT0$total
LLT0$total_dead <- LLT0$total - LLT0$total_alive
LLT0$EggMassClass <- "single"
LLT0$EggMass2 <- c(400,401,402,403,404)

LLTnew <- bind_rows(LLT100, LLT0)
LLTnew = LLTnew[,c(1,4,2,9,10,5:8)]


db. <- hatch_df_sub %>% dplyr::filter(TreatmentTime != "4584")

names(db.)[names(db.) == 'TotalEggs'] <- 'total'
names(db.)[names(db.) == 'TotalHatch'] <- 'total_alive'
names(db.)[names(db.) == 'PercentageHatch'] <- 'Prop_Survival'

db.$TreatmentTemp = as.character(db.$TreatmentTemp)
LLTnew$TreatmentTemp = as.character(LLTnew$TreatmentTemp)

db.$EggMass2 = as.character(db.$EggMass2)
LLTnew$EggMass2 = as.character(LLTnew$EggMass2)

db.$Prop_Survival = as.numeric(db.$Prop_Survival)
LLTnew$Prop_Survival = as.numeric(LLTnew$Prop_Survival)


db.merge <- bind_rows(db.,LLTnew)

db.sub <- db.merge 


y.yr <- as.character(unique(db.sub$YR))

list <- unique(db.sub$EggMass2)


for (i in 1:length(list)){
  E <- list[i]
  nr = as.numeric(db.sub[which(db.sub$EggMass2 == E, arr.ind=TRUE), which(colnames(db.sub)== "total_alive")])
  df0 <- data.frame(matrix(nrow = nr, ncol= 6), stringsAsFactors = FALSE) ## create blank dataframe
  colnames(df0) <- c("EggMass2", "EggMassClass","YR","TreatmentTime","TreatmentTemp","Fitness") ## name column
  df0$EggMass2 = as.character(df0$EggMass2)
  df0$EggMassClass = as.character(df0$EggMassClass)
  df0$TreatmentTime = as.numeric(df0$TreatmentTime)
  df0$TreatmentTemp = as.numeric(df0$TreatmentTemp)
  df0$YR = as.character(df0$YR)
  if(nr > 0) {
    df0$EggMass2 <- E
    df0$EggMassClass <- db.sub[which(db.sub$EggMass2 == E, arr.ind=TRUE),which(colnames(db.sub)== "EggMassClass")]
    df0$YR <- as.character(db.sub[which(db.sub$EggMass2 == E, arr.ind=TRUE),which(colnames(db.sub)== "YR")])
    df0$TreatmentTime <- as.numeric(db.sub[which(db.sub$EggMass2 == E, arr.ind=TRUE),which(colnames(db.sub)== "TreatmentTime")])
    df0$TreatmentTemp <- as.numeric(db.sub[which(db.sub$EggMass2 == E, arr.ind=TRUE),which(colnames(db.sub)== "TreatmentTemp")])
    df0$Fitness <- 1
  } 
  
  nrd = as.numeric(db.sub[which(db.sub$EggMass2 == E, arr.ind=TRUE),which(colnames(db.sub)=="total_dead" )])
  df1 <- data.frame(matrix(nrow = nrd, ncol= 6), stringsAsFactors = FALSE) ## create blank dataframe
  colnames(df1) <- c("EggMass2","EggMassClass","YR","TreatmentTime","TreatmentTemp","Fitness") ## name column
  df1$EggMass2 = as.character(df1$EggMass2)
  df1$EggMassClass = as.character(df1$EggMassClass)
  df1$TreatmentTime = as.numeric(df1$TreatmentTime)
  df1$TreatmentTemp = as.numeric(df1$TreatmentTemp)
  df1$YR = as.character(df1$YR)
  if(nrd > 0) {
    df1$EggMass2 <- E
    df1$EggMassClass <- db.sub[which(db.sub$EggMass2 == E, arr.ind=TRUE),which(colnames(db.sub)== "EggMassClass")]
    df1$YR <- as.character(db.sub[which(db.sub$EggMass2 == E, arr.ind=TRUE),which(colnames(db.sub)== "YR")])
    df1$TreatmentTime <- as.numeric(db.sub[which(db.sub$EggMass2 == E, arr.ind=TRUE),which(colnames(db.sub)== "TreatmentTime")])
    df1$TreatmentTemp <- as.numeric(db.sub[which(db.sub$EggMass2 == E, arr.ind=TRUE),which(colnames(db.sub)== "TreatmentTemp")])
    df1$Fitness <- 0
  } 
  
  # bind fitness datasets for temperature i at time j
  df_temp <- bind_rows(df0, df1) 
  
  # bind fitness datasets of all temperatures at time j
  if (exists("df_temp_all") == TRUE){
    df_temp_all <- bind_rows(df_temp_all, df_temp)
  } else {df_temp_all = df_temp}
  
}

bi.LLT.data.theo <- df_temp_all
check <- df_temp_all %>% group_by(YR,TreatmentTime,TreatmentTemp, Fitness)%>% summarise(n = n())



rm(list = ls()[!ls() %in% c("YR2_hatch_file", "bi.LLT.data", "bi.LLT.data.theo")])

