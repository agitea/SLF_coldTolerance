## SLF Overwintering Supercooling point graph
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
# install.packages("viridis")
library(viridis)
library(ggthemes)

theme_set(theme_bw())

options(stringsAsFactors=FALSE)

rm(list=ls())

egg_disturbance_file <- "data/3_SLF_2022_Y1_EggManipulation.csv" ## impact of egg disturbance on hatch rate data



#####################################################
## read data
egg_dist_df <- read.csv(paste0(egg_disturbance_file), header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE)

#####################################################
## reshape data
ControlEggs_df <- dplyr::select(egg_dist_df, EggMass, ControlEggs, ControlHatch)
RemovedEggs_df <- dplyr::select(egg_dist_df, EggMass,  RemovedEggs, RemovedHatch)

names(ControlEggs_df)[names(ControlEggs_df) == 'ControlHatch'] <- 'total_alive'
names(RemovedEggs_df)[names(RemovedEggs_df) == 'RemovedHatch'] <- 'total_alive'
names(ControlEggs_df)[names(ControlEggs_df) == 'ControlEggs'] <- 'total'
names(RemovedEggs_df)[names(RemovedEggs_df) == 'RemovedEggs'] <- 'total'
ControlEggs_df$total_dead <- ControlEggs_df$total - ControlEggs_df$total_alive
RemovedEggs_df$total_dead <- RemovedEggs_df$total - RemovedEggs_df$total_alive

ControlEggs_df$Type = "Control"
RemovedEggs_df$Type = "Removed"



# create fitness data frame
db.sub = ControlEggs_df


list <- unique(db.sub$EggMass)

for (i in 1:length(list)){
  t <- list[i]
  nr = as.numeric(db.sub[which(db.sub$EggMass == t, arr.ind=TRUE),which(colnames(db.sub)== "total_alive")])
  df0 <- data.frame(matrix(nrow = nr, ncol= 3), stringsAsFactors = FALSE) ## create blank dataframe
  colnames(df0) <- c("EggMass","Type","Fitness") ## name column
  df0$Type <- as.character(df0$Type)
  if(nr > 0) {
    df0$EggMass <- t
    df0$Type <- as.character(db.sub[which(db.sub$EggMass == t, arr.ind=TRUE),which(colnames(db.sub)== "Type")])
    df0$Fitness <- 1
  } 
  
  nrd = as.numeric(db.sub[which(db.sub$EggMass == t, arr.ind=TRUE),which(colnames(db.sub)=="total_dead" )])
  df1 <- data.frame(matrix(nrow = nrd, ncol= 3), stringsAsFactors = FALSE) ## create blank dataframe
  colnames(df1) <- c("EggMass","Type","Fitness") ## name column
  df1$Type <- as.character(df1$Type)
  if(nrd > 0) {
    df1$EggMass <- t
    df1$Type <- as.character(db.sub[which(db.sub$EggMass == t, arr.ind=TRUE),which(colnames(db.sub)== "Type")])
    df1$Fitness <- 0
  } 
  
  # bind fitness datasets for temperature i at time j
  df_temp <- bind_rows(df0, df1) 
  
  # bind fitness datasets of all temperatures at time j
  if (exists("df_temp_all") == TRUE){
    df_temp_all <- bind_rows(df_temp_all, df_temp)
  } else {df_temp_all = df_temp}
  
}



df_binary_control <- df_temp_all


# create fitness data frame
db.sub = RemovedEggs_df

list <- unique(db.sub$EggMass)

for (i in 1:length(list)){
  t <- list[i]
  nr = as.numeric(db.sub[which(db.sub$EggMass == t, arr.ind=TRUE),which(colnames(db.sub)== "total_alive")])
  df0 <- data.frame(matrix(nrow = nr, ncol= 3), stringsAsFactors = FALSE) ## create blank dataframe
  colnames(df0) <- c("EggMass","Type","Fitness") ## name column
  df0$Type <- as.character(df0$Type)
  if(nr > 0) {
    df0$EggMass <- t
    df0$Type <- as.character(db.sub[which(db.sub$EggMass == t, arr.ind=TRUE),which(colnames(db.sub)== "Type")])
    df0$Fitness <- 1
  } 
  
  nrd = as.numeric(db.sub[which(db.sub$EggMass == t, arr.ind=TRUE),which(colnames(db.sub)=="total_dead" )])
  df1 <- data.frame(matrix(nrow = nrd, ncol= 3), stringsAsFactors = FALSE) ## create blank dataframe
  colnames(df1) <- c("EggMass","Type","Fitness") ## name column
  df1$Type <- as.character(df1$Type)
  if(nrd > 0) {
    df1$EggMass <- t
    df1$Type <- as.character(db.sub[which(db.sub$EggMass == t, arr.ind=TRUE),which(colnames(db.sub)== "Type")])
    df1$Fitness <- 0
  } 
  
  # bind fitness datasets for temperature i at time j
  df_temp <- bind_rows(df0, df1) 
  
  # bind fitness datasets of all temperatures at time j
  if (exists("df_temp_all") == TRUE){
    df_temp_all <- bind_rows(df_temp_all, df_temp)
  } else {df_temp_all = df_temp}
  
}


df_binary_removed <- df_temp_all

# rm(list = ls()[!ls() %in% c("field_temperature_file", "egg_disturbance_file", "SCP_file", "YR2_hatch_file", "df_binary_control")])

df_binary <- bind_rows(df_binary_control, df_binary_removed)

df_binary$EggMass = as.character(df_binary$EggMass)


# glm function for logistic regression with logit link and binomial error distribution
GLM_dist_model <-
  glm(Fitness ~ Type,
      family = binomial(link = 'logit'),
      data = df_binary)

summary(GLM_dist_model)
anova(GLM_dist_model, test = "Chi")

pi <- exp(-0.31425)/ (1+ exp(-0.31425))
pir <- exp(-0.31425-1.10960)/ (1+ exp(-0.31425-1.10960))



####################################################################################################
# Egg disturbance - Fig S2 and analysis
####################################################################################################
rm(list = ls()[!ls() %in% c( "egg_disturbance_file")])

#####################################################
## read data
egg_dist_df <- read.csv(paste0(egg_disturbance_file), header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE)

#####################################################
## reshape data
egg_dist_df_sub <- dplyr::select(egg_dist_df, EggMass, HatchSuccessControl, HatchSuccessRemoved)
egg_dist_df2 <- gather(egg_dist_df_sub, "treatment", "value", 2:3)

#####################################################
## Fig.2 Percentage of hatch for intact egg masses of L. delicatula and eggs removed from 25 egg masses.
egg_dist_df2 %>%
  ggplot( aes(x=treatment, y=value, fill=treatment), alpha = 0.5) +
  geom_boxplot(notch = FALSE, notchwidth = 0.6) +
  # scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  # scale_fill_manual(values = c("#7b848f", "#a3cce9")) +
  scale_fill_manual(values = c("#7b848f", "#c8d0d9")) +
  geom_jitter(color="black", size=2, alpha=0.9) +
  # scale_y_log10()+
  theme_ipsum() +
  ylab("% hatch per egg mass")+
  NULL

#####################################################
## Paired t-test
# Save the data in two different vector
HatchSuccessControl <- egg_dist_df2 %>% filter(treatment == "HatchSuccessControl") 
HatchSuccessRemoved <- egg_dist_df2 %>% filter(treatment == "HatchSuccessRemoved") 
# Compute t-test
res <- t.test(HatchSuccessControl$value, HatchSuccessRemoved$value, paired = TRUE)
res





