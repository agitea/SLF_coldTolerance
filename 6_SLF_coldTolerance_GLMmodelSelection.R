## SLF Cold tolerance GLM model testing
## written by Anna J. Turbelin, Mar 14, 2024 
## aturbelin@gmail.com


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
library(dplyr)
library(tidyr)
# library(ggpmisc)
# library(ggeffects)

theme_set(theme_bw())

options(stringsAsFactors=FALSE)


source("4_SLF_coldTolerance_binomialDataset.R")

## Data file -------------------------------------------------------------------
YR2_hatch_file <- "data/SLF_YR1_YR2_hatch2.csv" ## cold tolerance strategy data 

rm(list = ls()[!ls() %in% c("YR2_hatch_file", "bi.LLT.data", "bi.LLT.data.theo")])


## raw data checking -----------------------------------------------------------

## boxplot of fitness against temperature in each treatment time
db <- bi.LLT.data
db$TreatmentTimeYR <- paste0(db$YR , "_", db$TreatmentTime)
db$Fitness <- as.character(db$Fitness)

db  %>%
  ggplot( aes(x=Fitness, y= TreatmentTemp, fill=Fitness)) +
  geom_boxplot()+
  geom_point(color = "#737373", alpha = .2, 
        position = position_jitter(width = 0.2, height = 0.4, seed = 2)) +
  facet_grid(rows = vars(TreatmentTimeYR)) 

################################################################################
## GLM raw data ----------------------------------------------------------------
## Effect of hurdled egg masses on hatch rate: 
## Is survival from egg masses with hurdled egg masses different from single egg masses in each treatment?

## YR2 All ---------------------------------------------------------------------
## Filter data 
db. <- bi.LLT.data %>% dplyr::filter(YR == "2")
db. <- db. %>% dplyr::filter(TreatmentTime != "4584")

db.$TreatmentTime <- as.character(db.$TreatmentTime)
db.$TreatmentTime <- factor(db.$TreatmentTime, levels = c("1", "12", "240", "360",ordered = TRUE))
db.$Treatment <- paste0(db.$TreatmentTime, "_", db.$TreatmentTemp)

# Adding EggMassClass as a fixed effect because we only have two categories

## probit link - time and temperature interactive effect 

LLT_model <-
  glm(Fitness ~ TreatmentTemp * TreatmentTime + EggMassClass,
      family = binomial(link = 'probit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")


E1 <- resid(LLT_model, type = "pearson")
sum(E1^2) / (LLT_model$df.res) # check for dispersion


resid(LLT_model) #List of residuals
plot(density(resid(LLT_model))) #A density plot
qqnorm(resid(LLT_model)) # A quantile normal plot - good for checking normality
qqline(resid(LLT_model))

drop1(LLT_model, test = "Chi")

plot( predict(LLT_model),resid(LLT_model) ) 

LLT_model %>%
  ggplot( aes(x=TreatmentTime, y=resid(LLT_model), fill=TreatmentTime)) +
  geom_boxplot()

LLT_model %>%
  ggplot( aes(x=TreatmentTime, y= Fitness, fill=TreatmentTime)) +
  geom_boxplot()+
  geom_point(aes(x = TreatmentTime, y = Fitness),
             color = "#4d004b", alpha = .2, 
             position = position_jitter(width = 0.2, height = 0.3, seed = 2))

LLT_model %>%
  ggplot( aes(x=EggMassClass, y=resid(LLT_model), fill=EggMassClass)) +
  geom_boxplot()

LLT_model %>%
  ggplot(aes(x=EggMassClass, y=Fitness, fill=EggMassClass)) +
  geom_boxplot() +
  geom_point(aes(x = EggMassClass, y = Fitness),
             color = "#4d004b", alpha = .2, 
             position = position_jitter(width = 0.2, height = 0.3, seed = 2))

## logit link - time and temperature interactive effect
LLT_model <-
  glm(Fitness ~ TreatmentTemp * TreatmentTime + EggMassClass,
      family = binomial(link = 'logit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")


## probit link - no time interaction
LLT_model <-
  glm(Fitness ~ TreatmentTemp + TreatmentTime + EggMassClass,
      family = binomial(link = 'probit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")


## logit link - no time interaction
LLT_model <-
  glm(Fitness ~ TreatmentTemp + TreatmentTime + EggMassClass,
      family = binomial(link = 'logit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")

## probit link - no egg class
LLT_model <-
  glm(Fitness ~ TreatmentTemp * TreatmentTime,
      family = binomial(link = 'probit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")


## logit link - no egg class
LLT_model <-
  glm(Fitness ~ TreatmentTemp * TreatmentTime,
      family = binomial(link = 'logit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")




## YR1 1H ----------------------------------------------------------------------
## Filter data -----------------------------------------------------------------
db. <- bi.LLT.data %>% dplyr::filter(YR == "1")


## probit link 
LLT_model <-
  glm(Fitness ~ TreatmentTemp,
      family = binomial(link = 'probit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")

## logit link 
LLT_model <-
  glm(Fitness ~ TreatmentTemp,
      family = binomial(link = 'logit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")





## YR2 1H ----------------------------------------------------------------------
## Filter data -----------------------------------------------------------------
db. <- bi.LLT.data %>% dplyr::filter(YR == "2")
db. <- db. %>% dplyr::filter(TreatmentTime == "1")

# Adding EggMassClass as a fixed effect because we only have two categories
## probit link + EggMassClass
LLT_model <-
  glm(Fitness ~ TreatmentTemp + EggMassClass,
      family = binomial(link = 'probit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")

## logit link + EggMassClass
LLT_model <-
  glm(Fitness ~ TreatmentTemp + EggMassClass,
      family = binomial(link = 'logit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")

## probit link 
LLT_model <-
  glm(Fitness ~ TreatmentTemp,
      family = binomial(link = 'probit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")

## logit link 
LLT_model <-
  glm(Fitness ~ TreatmentTemp,
      family = binomial(link = 'logit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")


## YR2 12H ----------------------------------------------------------------------
## Filter data -----------------------------------------------------------------
db. <- bi.LLT.data %>% dplyr::filter(YR == "2")
db. <- db. %>% dplyr::filter(TreatmentTime == "12")

# Adding EggMassClass as a fixed effect because we only have two categories
## probit link + EggMassClass
LLT_model <-
  glm(Fitness ~ TreatmentTemp + EggMassClass,
      family = binomial(link = 'probit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")

## logit link + EggMassClass
LLT_model <-
  glm(Fitness ~ TreatmentTemp + EggMassClass,
      family = binomial(link = 'logit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")

## probit link 
LLT_model <-
  glm(Fitness ~ TreatmentTemp,
      family = binomial(link = 'probit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")

## logit link 
LLT_model <-
  glm(Fitness ~ TreatmentTemp,
      family = binomial(link = 'logit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")



## YR2 240H ----------------------------------------------------------------------
## Filter data -----------------------------------------------------------------
db. <- bi.LLT.data %>% dplyr::filter(YR == "2")
db. <- db. %>% dplyr::filter(TreatmentTime == "240")

# Adding EggMassClass as a fixed effect because we only have two categories
## probit link + EggMassClass
LLT_model <-
  glm(Fitness ~ TreatmentTemp + EggMassClass,
      family = binomial(link = 'probit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")

## logit link + EggMassClass
LLT_model <-
  glm(Fitness ~ TreatmentTemp + EggMassClass,
      family = binomial(link = 'logit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")

## probit link 
LLT_model <-
  glm(Fitness ~ TreatmentTemp,
      family = binomial(link = 'probit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")

## logit link 
LLT_model <-
  glm(Fitness ~ TreatmentTemp,
      family = binomial(link = 'logit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")


## YR2 360H --------------------------------------------------------------------
## Filter data -----------------------------------------------------------------
db. <- bi.LLT.data %>% dplyr::filter(YR == "2")
db. <- db. %>% dplyr::filter(TreatmentTime == "360")

# Adding EggMassClass as a fixed effect because we only have two categories
## probit link + EggMassClass
LLT_model <-
  glm(Fitness ~ TreatmentTemp + EggMassClass,
      family = binomial(link = 'probit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")

## logit link + EggMassClass
LLT_model <-
  glm(Fitness ~ TreatmentTemp + EggMassClass,
      family = binomial(link = 'logit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")

## probit link 
LLT_model <-
  glm(Fitness ~ TreatmentTemp,
      family = binomial(link = 'probit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")

## logit link 
LLT_model <-
  glm(Fitness ~ TreatmentTemp,
      family = binomial(link = 'logit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")




################################################################################
## GLM theoretical data --------------------------------------------------------
## YR2 All ---------------------------------------------------------------------
## Filter data -----------------------------------------------------------------
db. <- bi.LLT.data.theo %>% dplyr::filter(YR == "2")
db. <- db. %>% dplyr::filter(TreatmentTime != "4584")

db.$TreatmentTime <- as.character(db.$TreatmentTime)
db.$TreatmentTime <- factor(db.$TreatmentTime, levels = c("1", "12", "240", "360",ordered = TRUE))
db.$Treatment <- paste0(db.$TreatmentTime, "_", db.$TreatmentTemp)

# Adding EggMassClass as a fixed effect because we only have two categories

## probit link - time and temperature interactive effect

LLT_model <-
  glm(Fitness ~ TreatmentTemp * TreatmentTime + EggMassClass,
      family = binomial(link = 'probit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")



## YR2 1H ----------------------------------------------------------------------
## Filter data -----------------------------------------------------------------
db. <- bi.LLT.data.theo %>% dplyr::filter(YR == "2")
db. <- db. %>% dplyr::filter(TreatmentTime == "1")

# Adding EggMassClass as a fixed effect because we only have two categories
## probit link + EggMassClass
LLT_model <-
  glm(Fitness ~ TreatmentTemp + EggMassClass,
      family = binomial(link = 'probit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")




## YR2 12H ---------------------------------------------------------------------
## Filter data -----------------------------------------------------------------
db. <- bi.LLT.data.theo %>% dplyr::filter(YR == "2")
db. <- db. %>% dplyr::filter(TreatmentTime == "12")

# Adding EggMassClass as a fixed effect because we only have two categories
## probit link + EggMassClass
LLT_model <-
  glm(Fitness ~ TreatmentTemp + EggMassClass,
      family = binomial(link = 'probit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")



## YR2 240H --------------------------------------------------------------------
## Filter data -----------------------------------------------------------------
db. <- bi.LLT.data.theo %>% dplyr::filter(YR == "2")
db. <- db. %>% dplyr::filter(TreatmentTime == "240")

# Adding EggMassClass as a fixed effect because we only have two categories
## probit link + EggMassClass
LLT_model <-
  glm(Fitness ~ TreatmentTemp + EggMassClass,
      family = binomial(link = 'probit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")



## YR2 360H --------------------------------------------------------------------
## Filter data -----------------------------------------------------------------
db. <- bi.LLT.data.theo %>% dplyr::filter(YR == "2")
db. <- db. %>% dplyr::filter(TreatmentTime == "360")

# Adding EggMassClass as a fixed effect because we only have two categories
## probit link + EggMassClass
LLT_model <-
  glm(Fitness ~ TreatmentTemp + EggMassClass,
      family = binomial(link = 'probit'),
      data = db.)

summary(LLT_model)
summ(LLT_model)
anova(LLT_model, test = "Chi")










