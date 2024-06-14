## SLF Cold tolerance 
## Effect of clustering on survival depending on treatment 
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
library(scales)
library(ggthemes)
library(Matrix)
library(lme4)
library(emmeans)
library(jtools)
library(effects)
# library(ggpmisc)
# library(ggeffects)

theme_set(theme_bw())

options(stringsAsFactors=FALSE)

rm(list = ls()[!ls() %in% c("YR2_hatch_file", "bi.LLT.data", "bi.LLT.data.theo")])

source("4_SLF_coldTolerance_binomialDataset.R")


## compare difference between groups  ---------------------------------------------------
rm(list = ls()[!ls() %in% c("YR2_hatch_file", "bi.LLT.data", "bi.LLT.data.theo")])

db_mod <- bi.LLT.data 

db_mod$Treatment <- paste0(db_mod$YR, "_", db_mod$TreatmentTime, "_", db_mod$TreatmentTemp)

y.yr <- as.character(unique(db_mod$Treatment))

for(y in 1:length(y.yr)){
  
  trt = y.yr[y] # get treatment
  db.1 <- db_mod %>% dplyr::filter(Treatment == trt) # filter data based on treatment

  class_list <- unique(db.1$EggMassClass)
  
  if ("multiple" %in% class_list) {
      
      # glm function for logistic regression with probit link and binomial error distribution
      LLT_model <-
        glm(Fitness ~ EggMassClass,
            family = binomial(link = "logit"),
            data = db.1)
      
      summary(LLT_model) # get summary of model - not necessary in this function

      # get coefficients
      b0 <- coef(LLT_model)[1] # get intercept
      b1 <- coef(LLT_model)[2] # get coef for single adjusted 
      
      pi <- exp(b0)/(1+exp(b0))
      pi_single <- exp(b0 + b1) / (1+exp(b0 + b1))
      
      atest <- anova(LLT_model, test = 'Chi')
      
      output <- data.frame(matrix(nrow = 1, ncol= 6), stringsAsFactors = FALSE) ## create blank dataframe
      colnames(output) <- c("TreatmentTime","TreatmentTemp","Treatment", "Survival_multi", "Suvival_single", "Chi_pvalue") ## name column
      output$TreatmentTime <- unique(db.1$TreatmentTime)
      output$TreatmentTemp <- unique(db.1$TreatmentTemp)
      output$Treatment <- trt
      output$Survival_multi <- pi 
      output$Suvival_single <- pi_single
      output$Chi_pvalue <- atest$`Pr(>Chi)`[2]
      
      
      
      if (exists("output.all") == TRUE){
        output.all <- bind_rows(output.all, output)
      } else {output.all = output}
  } 
  
    }
 


# Table S8







