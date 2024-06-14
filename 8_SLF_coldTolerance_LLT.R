## SLF Cold tolerance and lower lethal limits analysis
## Effect of time and temperature on survival and GLMM model to estimate LLT
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

## clear environment except from files needed for analysis ---------------------
rm(list = ls())

source("4_SLF_coldTolerance_binomialDataset.R")
source("adjusted_doseP_fun.R")


rm(list = ls()[!ls() %in% c("YR2_hatch_file", "bi.LLT.data", "bi.LLT.data.theo", "calculate_se_single", "dose.adj")])

################################################################################
## GLM raw data ----------------------------------------------------------------
## Effect of hurdled egg masses on hatch rate: 
## Is survival from egg masses with hurdled egg masses different from single egg masses in each treatment?

## Filter data -----------------------------------------------------------------
db. <- bi.LLT.data %>% dplyr::filter(YR == "2")
db. <- db. %>% dplyr::filter(TreatmentTime != "4584")

db.$TreatmentTime <- as.character(db.$TreatmentTime)
db.$TreatmentTime <- factor(db.$TreatmentTime, levels = c("1", "12", "240", "360",ordered = TRUE))
db.$Treatment <- paste0(db.$TreatmentTime, "_", db.$TreatmentTemp)

# Adding EggMassClass as a fixed effect because we only have two categories

LLT_model <-
  glm(Fitness ~ TreatmentTemp * TreatmentTime + EggMassClass,
      family = binomial(link = 'probit'),
      data = db.)

summary(LLT_model)

## 1 ---------------------------------------------------------------------------
obj <- LLT_model
cf <- c(1,2)
p <- c(0.1, 0.5, 0.9)
tp <- NULL
# terms <- c("(Intercept)", "TreatmentTemp")

Time1 <- dose.adj(obj, cf = cf, p = p, tp = NULL, terms = NULL)
Time1$TreatmentTime <- "1"


## 12 --------------------------------------------------------------------------
obj <- LLT_model
cf <- c(1,2)
p <- c(0.1, 0.5, 0.9)
tp <- NULL
terms <- c("(Intercept)", "TreatmentTemp", "TreatmentTime12", "TreatmentTemp:TreatmentTime12")

Time12 <- dose.adj(obj, cf = cf, p = p, tp = NULL, terms = terms)
Time12$TreatmentTime <- "12"


## 240 -------------------------------------------------------------------------
obj <- LLT_model
cf <- c(1,2)
p <- c(0.1, 0.5, 0.9)
tp <- NULL
terms <- c("(Intercept)", "TreatmentTemp", "TreatmentTime240", "TreatmentTemp:TreatmentTime240")

Time240 <- dose.adj(obj, cf = cf, p = p, tp = NULL, terms = terms)
Time240$TreatmentTime <- "240"


## 360 -------------------------------------------------------------------------
obj <- LLT_model
cf <- c(1,2)
p <- c(0.1, 0.5, 0.9)
tp <- NULL
terms <- c("(Intercept)", "TreatmentTemp", "TreatmentTime360", "TreatmentTemp:TreatmentTime360")

Time360 <- dose.adj(obj, cf = cf, p = p, tp = NULL, terms = terms)
Time360$TreatmentTime <- "360"


LLT_raw_ouput <- bind_rows(Time1,Time12,Time240, Time360)


LLT_raw_ouput$EggMassClass <- "multiple"

LLT_raw_ouput$dataType <- "raw"


################################################################################
## GLM Theoretical data ----------------------------------------------------------------
## Effect of hurdled egg masses on hatch rate: 
## Is survival from egg masses with hurdled egg masses different from single egg masses in each treatment?

## Filter data -----------------------------------------------------------------
db. <- bi.LLT.data.theo %>% dplyr::filter(YR == "2")
db. <- db. %>% dplyr::filter(TreatmentTime != "4584")

db.$TreatmentTime <- as.character(db.$TreatmentTime)
db.$TreatmentTime <- factor(db.$TreatmentTime, levels = c("1", "12", "240", "360",ordered = TRUE))
db.$Treatment <- paste0(db.$TreatmentTime, "_", db.$TreatmentTemp)

# Adding EggMassClass as a fixed effect because we only have two categories

LLT_model <-
  glm(Fitness ~ TreatmentTemp * TreatmentTime + EggMassClass,
      family = binomial(link = 'probit'),
      data = db.)

summary(LLT_model)

## 1 ---------------------------------------------------------------------------
obj <- LLT_model
cf <- c(1,2)
p <- c(0.1, 0.5, 0.9)
tp <- NULL
# terms <- c("(Intercept)", "TreatmentTemp")

Time1 <- dose.adj(obj, cf = cf, p = p, tp = NULL, terms = NULL)
Time1$TreatmentTime <- "1"


## 12 --------------------------------------------------------------------------
obj <- LLT_model
cf <- c(1,2)
p <- c(0.1, 0.5, 0.9)
tp <- NULL
terms <- c("(Intercept)", "TreatmentTemp", "TreatmentTime12", "TreatmentTemp:TreatmentTime12")
# terms <- c("(Intercept)", "TreatmentTemp", "TreatmentTime12")

Time12 <- dose.adj(obj, cf = cf, p = p, tp = NULL, terms = terms)
Time12$TreatmentTime <- "12"


## 240 -------------------------------------------------------------------------
obj <- LLT_model
cf <- c(1,2)
p <- c(0.1, 0.5, 0.9)
tp <- NULL
terms <- c("(Intercept)", "TreatmentTemp", "TreatmentTime240", "TreatmentTemp:TreatmentTime240")
# terms <- c("(Intercept)", "TreatmentTemp", "TreatmentTime240")

Time240 <- dose.adj(obj, cf = cf, p = p, tp = NULL, terms = terms)
Time240$TreatmentTime <- "240"


## 360 -------------------------------------------------------------------------
obj <- LLT_model
cf <- c(1,2)
p <- c(0.1, 0.5, 0.9)
tp <- NULL
terms <- c("(Intercept)", "TreatmentTemp", "TreatmentTime360", "TreatmentTemp:TreatmentTime360")

Time360 <- dose.adj(obj, cf = cf, p = p, tp = NULL, terms = terms)
Time360$TreatmentTime <- "360"


LLT_theo_ouput <- bind_rows(Time1,Time12,Time240, Time360)


LLT_theo_ouput$EggMassClass <- "multiple"


LLT_theo_ouput$dataType <- "theoretical"




LLT_multi <- bind_rows(LLT_raw_ouput, LLT_theo_ouput)


## clear environment except from files needed for analysis ---------------------
rm(list = ls()[!ls() %in% c("YR2_hatch_file", "bi.LLT.data", "bi.LLT.data.theo", "LLT_multi", "calculate_se_single", "dose.adj")])


################################################################################
## SINGLE egg mass
## GLM raw data ----------------------------------------------------------------
## Effect of hurdled egg masses on hatch rate: 
## Is survival from egg masses with hurdled egg masses different from single egg masses in each treatment?

## Filter data -----------------------------------------------------------------
db. <- bi.LLT.data %>% dplyr::filter(YR == "2")
db. <- db. %>% dplyr::filter(TreatmentTime != "4584")

db.$TreatmentTime <- as.character(db.$TreatmentTime)
db.$TreatmentTime <- factor(db.$TreatmentTime, levels = c("1", "12", "240", "360",ordered = TRUE))
db.$Treatment <- paste0(db.$TreatmentTime, "_", db.$TreatmentTemp)

# Adding EggMassClass as a fixed effect because we only have two categories

LLT_model <-
  glm(Fitness ~ TreatmentTemp * TreatmentTime + EggMassClass,
      family = binomial(link = 'probit'),
      data = db.)

summary(LLT_model)


## 1 ---------------------------------------------------------------------------
obj <- LLT_model
cf <- c(1,2)
p <- c(0.1, 0.5, 0.9)
tp <- NULL
terms <- c("(Intercept)", "TreatmentTemp", "EggMassClasssingle")

Time1 <- dose.adj(obj, cf = cf, p = p, tp = NULL, terms = terms)
Time1$TreatmentTime <- "1"


## 12 --------------------------------------------------------------------------
obj <- LLT_model
cf <- c(1,2)
p <- c(0.1, 0.5, 0.9)
tp <- NULL
terms <- c("(Intercept)", "TreatmentTemp", "TreatmentTime12", "TreatmentTemp:TreatmentTime12","EggMassClasssingle")

Time12 <- dose.adj(obj, cf = cf, p = p, tp = NULL, terms = terms)
Time12$TreatmentTime <- "12"


## 240 -------------------------------------------------------------------------
obj <- LLT_model
cf <- c(1,2)
p <- c(0.1, 0.5, 0.9)
tp <- NULL
terms <- c("(Intercept)", "TreatmentTemp", "TreatmentTime240", "TreatmentTemp:TreatmentTime240","EggMassClasssingle")

Time240 <- dose.adj(obj, cf = cf, p = p, tp = NULL, terms = terms)
Time240$TreatmentTime <- "240"


## 360 -------------------------------------------------------------------------
obj <- LLT_model
cf <- c(1,2)
p <- c(0.1, 0.5, 0.9)
tp <- NULL
terms <- c("(Intercept)", "TreatmentTemp", "TreatmentTime360", "TreatmentTemp:TreatmentTime360","EggMassClasssingle")

Time360 <- dose.adj(obj, cf = cf, p = p, tp = NULL, terms = terms)
Time360$TreatmentTime <- "360"


LLT_raw_ouput <- bind_rows(Time1,Time12,Time240, Time360)


LLT_raw_ouput$EggMassClass <- "single"

LLT_raw_ouput$dataType <- "raw"


################################################################################
## GLM Theoretical data ----------------------------------------------------------------
## Effect of hurdled egg masses on hatch rate: 
## Is survival from egg masses with hurdled egg masses different from single egg masses in each treatment?

## Filter data -----------------------------------------------------------------
db. <- bi.LLT.data.theo %>% dplyr::filter(YR == "2")
db. <- db. %>% dplyr::filter(TreatmentTime != "4584")

db.$TreatmentTime <- as.character(db.$TreatmentTime)
db.$TreatmentTime <- factor(db.$TreatmentTime, levels = c("1", "12", "240", "360",ordered = TRUE))
db.$Treatment <- paste0(db.$TreatmentTime, "_", db.$TreatmentTemp)

# Adding EggMassClass as a fixed effect because we only have two categories

LLT_model <-
  glm(Fitness ~ TreatmentTemp * TreatmentTime + EggMassClass,
      family = binomial(link = 'probit'),
      data = db.)


## 1 ---------------------------------------------------------------------------
obj <- LLT_model
cf <- c(1,2)
p <- c(0.1, 0.5, 0.9)
tp <- NULL
terms <- c("(Intercept)", "TreatmentTemp", "EggMassClasssingle")

Time1 <- dose.adj(obj, cf = cf, p = p, tp = NULL, terms = NULL)
Time1$TreatmentTime <- "1"


## 12 --------------------------------------------------------------------------
obj <- LLT_model
cf <- c(1,2)
p <- c(0.1, 0.5, 0.9)
tp <- NULL
terms <- c("(Intercept)", "TreatmentTemp", "TreatmentTime12", "TreatmentTemp:TreatmentTime12","EggMassClasssingle")

Time12 <- dose.adj(obj, cf = cf, p = p, tp = NULL, terms = terms)
Time12$TreatmentTime <- "12"


## 240 -------------------------------------------------------------------------
obj <- LLT_model
cf <- c(1,2)
p <- c(0.1, 0.5, 0.9)
tp <- NULL
terms <- c("(Intercept)", "TreatmentTemp", "TreatmentTime240", "TreatmentTemp:TreatmentTime240","EggMassClasssingle")

Time240 <- dose.adj(obj, cf = cf, p = p, tp = NULL, terms = terms)
Time240$TreatmentTime <- "240"


## 360 -------------------------------------------------------------------------
obj <- LLT_model
cf <- c(1,2)
p <- c(0.1, 0.5, 0.9)
tp <- NULL
terms <- c("(Intercept)", "TreatmentTemp", "TreatmentTime360", "TreatmentTemp:TreatmentTime360", "EggMassClasssingle")

Time360 <- dose.adj(obj, cf = cf, p = p, tp = NULL, terms = terms)
Time360$TreatmentTime <- "360"

## bind output -----------------------------------------------------------------
LLT_theo_ouput <- bind_rows(Time1,Time12,Time240, Time360)


LLT_theo_ouput$EggMassClass <- "single"


LLT_theo_ouput$dataType <- "theoretical"




LLT_single <- bind_rows(LLT_raw_ouput, LLT_theo_ouput)



## bind all --------------------------------------------------------------------
LLT_output <- bind_rows(LLT_multi, LLT_single)



LLT_output_mean <- data.frame(matrix(nrow = nrow(LLT_single), ncol= ncol(LLT_single)), stringsAsFactors = FALSE) ## create blank dataframe
colnames(LLT_output_mean) <- colnames(LLT_single) ## name column

LLT_output_mean$Fitness <- LLT_single$Fitness
LLT_output_mean$TreatmentTime <- LLT_single$TreatmentTime
LLT_output_mean$dataType <- LLT_single$dataType
LLT_output_mean$TreatmentTemp <- (LLT_single$TreatmentTemp + LLT_multi$TreatmentTemp) /2
LLT_output_mean$SE <- (LLT_single$SE + LLT_multi$SE) /2
LLT_output_mean$SEup <- (LLT_single$SEup + LLT_multi$SEup) /2
LLT_output_mean$SElo <- (LLT_single$SElo + LLT_multi$SElo) /2
LLT_output_mean$se_pi <- (LLT_single$se_pi + LLT_multi$se_pi) /2
LLT_output_mean$LowerBoundCI <- (LLT_single$LowerBoundCI + LLT_multi$LowerBoundCI) /2
LLT_output_mean$UpperBoundCI <- (LLT_single$UpperBoundCI + LLT_multi$UpperBoundCI) /2


LLT_output_mean_sub <- LLT_output_mean %>%
  dplyr::select(dataType, Fitness, TreatmentTime, TreatmentTemp, SE)

LLT_output_mean_sub$TempSE <- paste0(round(LLT_output_mean_sub$TreatmentTemp,1), "(+- ", round(LLT_output_mean_sub$SE,1) ,")")

LLT_output_mean_sub <- LLT_output_mean_sub %>%
  dplyr::select(dataType, Fitness, TreatmentTime, TempSE)

LLT_output_mean_sub <- LLT_output_mean_sub %>% tidyr::spread(Fitness,TempSE)




LLT_multi_sub <- LLT_multi %>%
  dplyr::select(dataType, Fitness, TreatmentTime, TreatmentTemp, SE)

LLT_multi_sub$TempSE <- paste0(round(LLT_multi_sub$TreatmentTemp, 1), "(+- ", round(LLT_multi_sub$SE,1) ,")")

LLT_multi_sub <- LLT_multi_sub %>%
  dplyr::select(dataType, Fitness, TreatmentTime, TempSE)

LLT_multi_sub <- LLT_multi_sub %>% tidyr::spread(Fitness,TempSE)


LLT_single_sub <- LLT_single %>%
  dplyr::select(dataType, Fitness, TreatmentTime, TreatmentTemp, SE)

LLT_single_sub$TempSE <- paste0(round(LLT_single_sub$TreatmentTemp,1), "(+- ", round(LLT_single_sub$SE, 1) ,")")

LLT_single_sub <- LLT_single_sub %>%
  dplyr::select(dataType, Fitness, TreatmentTime, TempSE)

LLT_single_sub <- LLT_single_sub %>% tidyr::spread(Fitness,TempSE)





####################################################################################################
# Lower Lethal limits analysis - LLT using probit link generalized model Figure
####################################################################################################

## clear environment except from files needed for analysis ---------------------
rm(list = ls()[!ls() %in% c("YR2_hatch_file", "bi.LLT.data", "bi.LLT.data.theo", "calculate_se_single", "dose.adj")])

## plot LLT model data (raw) ---------------------------------------------------
db_mod <- bi.LLT.data %>% dplyr::filter(TreatmentTime != "4584") # filter data to remove control

y.yr <- as.character(unique(db_mod$YR))

for(y in 1:length(y.yr)){
  
  yr = y.yr[y] # get year 
  db.1 <- db_mod %>% dplyr::filter(YR == yr) # filter data based on year
  t.time <- unique(db.1$TreatmentTime) # get list of treatment times in year data
  p = seq(0, 1, length = 1001) # generate survival probability sequence - change to p = c(0.1, 0.5, 0.9) to get LLT 90, 50 and 10.
  
  if( yr == 1){
    
    for (j in 1:length(t.time)){
      
      time = t.time[j] # get time
      db.sub <- db.1 %>% dplyr::filter(TreatmentTime == time) # filter data based on time
      
      # glm function for logistic regression with probit link and binomial error distribution
      LLTemp_model <-
        glm(Fitness ~ TreatmentTemp,
            family = binomial(link = "probit"),
            data = db.sub)
      
      # summary(LLTemp_model) # get summary of model - not necessary in this function
      
      # calculate LLT
      m.output <- dose.adj(LLTemp_model, p = p) # use adjusted dose function to get LLT, SE and 95% CI
      
      m.output$YR = yr  # add year to the dataframe
      m.output$TreatmentTime = time # add time to the dataframe
      
      if (exists("m.output.all") == TRUE){
        m.output.all <- bind_rows(m.output.all, m.output)
      } else {m.output.all = m.output}
      
    }
  } else {
    
    for (j in 1:length(t.time)){
      
      time = t.time[j]
      db.sub <- db.1 %>% dplyr::filter(TreatmentTime == time)
      
      # glm function for logistic regression with probit link and binomial error distribution
      LLTemp_model <-
        glm(Fitness ~ TreatmentTemp + EggMassClass,
            family = binomial(link = "probit"),
            data = db.sub)
      
      # summary(LLTemp_model) # get summary of model - not necessary in this function
      
      # calculate LLT
      m.output <- dose.adj(LLTemp_model, p = p)
      
      m.output$YR = yr
      m.output$TreatmentTime = time
      
      
      if (exists("m.output.all") == TRUE){
        m.output.all <- bind_rows(m.output.all, m.output)
      } else {m.output.all = m.output}
      
    }
    
  }
  
}


m.output.all.raw <- m.output.all


## plot LLT model data (theoretical) -------------------------------------------

rm(list = ls()[!ls() %in% c("YR2_hatch_file","m.output.all.raw", "bi.LLT.data", "bi.LLT.data.theo",  "calculate_se_single", "dose.adj")])

# db. <- bi.LLT.data %>% dplyr::filter(YR == "2")
db_mod <- bi.LLT.data.theo %>% dplyr::filter(TreatmentTime != "4584")

y.yr <- as.character(unique(db_mod$YR))

for(y in 1:length(y.yr)){
  
  yr = y.yr[y]
  db.1 <- db_mod %>% dplyr::filter(YR == yr)
  t.time <- unique(db.1$TreatmentTime)
  p = seq(0, 1, length = 1001)
  
  if( yr == 1){
    
    for (j in 1:length(t.time)){
      
      time = t.time[j]
      db.sub <- db.1 %>% dplyr::filter(TreatmentTime == time)
      
      # glm function for logistic regression with logit link and binomial error distribution
      LLTemp_model <-
        glm(Fitness ~ TreatmentTemp,
            family = binomial(link = "probit"),
            data = db.sub)
      
      summary(LLTemp_model)
      
      # calculation of LT50 and LT90
      m.output <- dose.adj(LLTemp_model, p = p)
      
      m.output$YR = yr
      m.output$TreatmentTime = time
      
      if (exists("m.output.all") == TRUE){
        m.output.all <- bind_rows(m.output.all, m.output)
      } else {m.output.all = m.output}
      
    }
  } else {
    
    for (j in 1:length(t.time)){
      
      time = t.time[j]
      db.sub <- db.1 %>% dplyr::filter(TreatmentTime == time)
      
      # glm function for logistic regression with logit link and binomial error distribution
      LLTemp_model <-
        glm(Fitness ~ TreatmentTemp + EggMassClass,
            family = binomial(link = "probit"),
            data = db.sub)
      
      summary(LLTemp_model)
      
      # calculation of LT50 and LT90
      m.output <- dose.adj(LLTemp_model, p = p)
      
      m.output$YR = yr
      m.output$TreatmentTime = time
      
      
      if (exists("m.output.all") == TRUE){
        m.output.all <- bind_rows(m.output.all, m.output)
      } else {m.output.all = m.output}
      
    }
    
  }
  
}

m.output.all.theo <- m.output.all


## Merge datasets --------------------------------------------------------------
m.output.all.theo$LTtheory = "theoretical"
m.output.all.raw$LTtheory = "raw"
model.outputs = bind_rows(m.output.all.raw, m.output.all.theo)


## Figure S6 - make plot -------------------------------------------------------

db. <-  bi.LLT.data.theo %>% dplyr::filter(TreatmentTime != "4584")
db.$TreatmentTime <- as.character(db.$TreatmentTime)
db.$TreatmentTimeYR <- paste0(db.$YR,"_",db.$TreatmentTime)
model.outputs$TreatmentTimeYR <- paste0(model.outputs$YR,"_",model.outputs$TreatmentTime)
m.output.all.raw$TreatmentTimeYR <- paste0(m.output.all.raw$YR,"_",m.output.all.raw$TreatmentTime)
m.output.all.theo$TreatmentTimeYR <- paste0(m.output.all.theo$YR,"_",m.output.all.theo$TreatmentTime)

cairo_pdf(file = "output/FigS6_SFlTempTime_GLM_raw_theo_facet_v2.pdf", width = 8.5, height = 10)

model.outputs %>%
  ggplot(aes(x = TreatmentTemp, y = Fitness, color = LTtheory)) +
  geom_ribbon(aes(ymin = LowerBoundCI, ymax = UpperBoundCI, fill = LTtheory, alpha = .5)) +
  geom_line() +
  geom_point(data= db., aes(x = TreatmentTemp, y = Fitness),color = "#1170aa", alpha = .2, position = position_jitter(width = 0.9, height = 0.001, seed = 2))+
  theme(legend.position = "none") +
  scale_color_manual(values= c("#5fa2ce","#57606c"))+
  scale_fill_manual(values= c("#a3cce9","#57606c")) +
  scale_x_continuous(limits = c(-50,50), breaks = c(-40,-30,-27.7,-25,-22.5, -20,-15,-10,-5,0,5,10,20,30,40)) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.1,0.5,0.9,1))+
  facet_grid(rows = vars(TreatmentTimeYR)) +
  ggtitle("LLT ") +
  NULL

dev.off()







