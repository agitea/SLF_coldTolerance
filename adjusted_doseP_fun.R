## dose.p() adjusted function
## This function is used to calculate the lower lethal temperatures for spotted lanternfly
## Function would need to be adjusted if using terms not present in the column names of our dataset
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
# library(plyr)
# detach(package:plyr)    
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

## clear environment except from files needed for analysis ---------------------
rm(list = ls()[!ls() %in% c("YR2_hatch_file", "bi.LLT.data", "bi.LLT.data.theo")])

##------------------------------------------------------------------------------
## Adjusted dose.p function
##------------------------------------------------------------------------------

## Function to calculate standard error of response variable -------------------
calculate_se_single <- function(x.p, obj, cf) {
  vcov_matrix <- vcov(obj)
  se <- sqrt((x.p^2 * vcov_matrix[cf, cf]))
  return(se)
}




## Modified dose.p function to extract and calculate LLT, SE and 95% CI --------
dose.adj <- function (obj, cf = 1:2, p = NULL, tp = NULL, terms = NULL)
{if(is.null(tp)){
  if(is.null(terms)){
    eta <- family(obj)$linkfun(p)
    b <- coef(obj)[cf]
    x.p <- (eta - b[1L])/b[2L]
    pd <- -cbind(1, x.p)/b[2L]
    SE <- sqrt(((pd %*% vcov(obj)[cf, cf]) * pd) %*% c(1, 1))
    z <- qnorm(0.975)
    SEup <- pnorm(eta - z * SE)
    SElo <- pnorm(eta - z * SE)
    
    # Calculate the standard error
    se_values <- sapply(x.p, calculate_se_single, obj = obj, cf = cf)
    se_pi <- se_values[4,]
    
    
    LowerBoundCI <- pnorm(eta - z * se_pi)
    UpperBoundCI <- pnorm(eta + z * se_pi)
    
    se_pi <- matrix(unlist(se_pi), nrow = length(se_pi), byrow = TRUE)
    LowerBoundCI <- matrix(unlist(LowerBoundCI), nrow = length(LowerBoundCI), byrow = TRUE)
    UpperBoundCI <- matrix(unlist(UpperBoundCI), nrow = length(UpperBoundCI), byrow = TRUE)
    
    
    re <- bind_cols(p, eta, x.p, SE, SEup, SElo,  se_pi, LowerBoundCI,UpperBoundCI)
    colnames(re) <- c("Fitness","eta","TreatmentTemp","SE", "SEup", "SElo","se_pi" ,"LowerBoundCI", "UpperBoundCI") ## name column
    re <- as.data.frame(re)
    return(re)
    
  } else {
    if ("EggMassClasssingle" %in% terms) {
      if(length(terms) < 4){
        b <- coef(LLT_model)[terms]
        b1 <- coef(LLT_model)[terms[1]] # intercept
        b2 <- coef(LLT_model)[terms[2]] # temperature
        if(length(terms) ==3 ){
          b3 <- coef(LLT_model)[terms[3]] # egg mass type
        } else { b3 <- 0 }
        
        
        x.p <- (qnorm(p) - b1 -  b3) / b2
        eta <- b1 + (x.p * b2) + b3
        pd <- -cbind(1, x.p)/b[2L]
        cf2 <- c(1,2)
        SE <- sqrt(((pd %*% vcov(obj)[cf2, cf2]) * pd) %*% c(1, 1))
        z <- qnorm(0.975)
        SEup <- pnorm(eta - z * SE)
        SElo <- pnorm(eta - z * SE)
        
        # Calculate the standard error
        se_values <- sapply(x.p, calculate_se_single, obj = obj, cf = cf)
        se_pi <- se_values[4,]
        
        
        LowerBoundCI <- pnorm(eta - z * se_pi)
        UpperBoundCI <- pnorm(eta + z * se_pi)
        
        se_pi <- matrix(unlist(se_pi), nrow = length(se_pi), byrow = TRUE)
        LowerBoundCI <- matrix(unlist(LowerBoundCI), nrow = length(LowerBoundCI), byrow = TRUE)
        UpperBoundCI <- matrix(unlist(UpperBoundCI), nrow = length(UpperBoundCI), byrow = TRUE)
        
        
        re <- bind_cols(p, eta, x.p, SE, SEup, SElo,  se_pi, LowerBoundCI,UpperBoundCI)
        colnames(re) <- c("Fitness","eta","TreatmentTemp","SE", "SEup", "SElo","se_pi" ,"LowerBoundCI", "UpperBoundCI") ## name column
        re <- as.data.frame(re)
        return(re)
        
      }
      if(length(terms) == 4){
        b <- coef(LLT_model)[terms]
        b1 <- coef(LLT_model)[terms[1]] # intercept
        b2 <- coef(LLT_model)[terms[2]] # temp coef
        b3 <- coef(LLT_model)[terms[3]] # time coef
        b4 <- coef(LLT_model)[terms[4]] # egg mass type

        x.p <- (qnorm(p) - b1 - b3 - b4) / b2 
        
        eta <- b1 + (x.p * b2 ) + b3 + b4
        
        pd <- -cbind(1, x.p)/b[2L]
        cf2 <- c(1,2)
        SE <- sqrt(((pd %*% vcov(obj)[cf2, cf2]) * pd) %*% c(1, 1))
        z <- qnorm(0.975)
        SEup <- pnorm(eta - z * SE)
        SElo <- pnorm(eta - z * SE)
        
        # Calculate the standard error
        se_values <- sapply(x.p, calculate_se_single, obj = obj, cf = cf)
        se_pi <- se_values[4,]
        
        
        LowerBoundCI <- pnorm(eta - z * se_pi)
        UpperBoundCI <- pnorm(eta + z * se_pi)
        
        se_pi <- matrix(unlist(se_pi), nrow = length(se_pi), byrow = TRUE)
        LowerBoundCI <- matrix(unlist(LowerBoundCI), nrow = length(LowerBoundCI), byrow = TRUE)
        UpperBoundCI <- matrix(unlist(UpperBoundCI), nrow = length(UpperBoundCI), byrow = TRUE)
        
        
        re <- bind_cols(p, eta, x.p, SE, SEup, SElo,  se_pi, LowerBoundCI,UpperBoundCI)
        colnames(re) <- c("Fitness","eta","TreatmentTemp","SE", "SEup", "SElo","se_pi" ,"LowerBoundCI", "UpperBoundCI") ## name column
        re <- as.data.frame(re)
        return(re)
        
      }
      if(length(terms) > 4){
        b <- coef(LLT_model)[terms]
        b1 <- coef(LLT_model)[terms[1]] # intercept
        b2 <- coef(LLT_model)[terms[2]] # temp coef
        b3 <- coef(LLT_model)[terms[3]] # time coef
        b4 <- coef(LLT_model)[terms[4]] # temp/time intercept
        b5 <- coef(LLT_model)[terms[5]] # egg mass type
        
        x.p <- (qnorm(p) - b1 - b3 - b5) / (b2 + b4)
        
        eta <- b1 + (x.p * (b2 + b4)) + b3 + b5
        
        pd <- -cbind(1, x.p)/b[2L]
        cf2 <- c(1,2)
        SE <- sqrt(((pd %*% vcov(obj)[cf2, cf2]) * pd) %*% c(1, 1))
        z <- qnorm(0.975)
        SEup <- pnorm(eta - z * SE)
        SElo <- pnorm(eta - z * SE)
        
        # Calculate the standard error
        se_values <- sapply(x.p, calculate_se_single, obj = obj, cf = cf)
        se_pi <- se_values[4,]
        
        
        LowerBoundCI <- pnorm(eta - z * se_pi)
        UpperBoundCI <- pnorm(eta + z * se_pi)
        
        se_pi <- matrix(unlist(se_pi), nrow = length(se_pi), byrow = TRUE)
        LowerBoundCI <- matrix(unlist(LowerBoundCI), nrow = length(LowerBoundCI), byrow = TRUE)
        UpperBoundCI <- matrix(unlist(UpperBoundCI), nrow = length(UpperBoundCI), byrow = TRUE)
        
        
        re <- bind_cols(p, eta, x.p, SE, SEup, SElo,  se_pi, LowerBoundCI,UpperBoundCI)
        colnames(re) <- c("Fitness","eta","TreatmentTemp","SE", "SEup", "SElo","se_pi" ,"LowerBoundCI", "UpperBoundCI") ## name column
        re <- as.data.frame(re)
        return(re)
        
      }
      
      
      
    } else {

      b <- coef(LLT_model)[terms]
      b1 <- coef(LLT_model)[terms[1]] # intercept
      b2 <- coef(LLT_model)[terms[2]] # main response variable - here is temperature
      b3 <- coef(LLT_model)[terms[3]] # Time coeficient
      
      if(length(terms) >3){
        b4 <- coef(LLT_model)[terms[4]] # interaction coeficient
        # x.p <- (qnorm(p) - b[1L] -  b[3L]) / ((b[2L] + b[4L]))
        x.p <- (qnorm(p) - b1 -  b3) / ((b2 + b4))
        
        eta <- b1 + (x.p * (b2 + b4)) + b3
        
        pd <- -cbind(1, x.p)/b[2L]
        cf2 <- c(1,2)
        SE <- sqrt(((pd %*% vcov(obj)[cf2, cf2]) * pd) %*% c(1, 1))
        z <- qnorm(0.975)
        SEup <- pnorm(eta - z * SE)
        SElo <- pnorm(eta - z * SE)
        
        # Calculate the standard error
        se_values <- sapply(x.p, calculate_se_single, obj = obj, cf = cf)
        se_pi <- se_values[4,]
        
        
        LowerBoundCI <- pnorm(eta - z * se_pi)
        UpperBoundCI <- pnorm(eta + z * se_pi)
        
        se_pi <- matrix(unlist(se_pi), nrow = length(se_pi), byrow = TRUE)
        LowerBoundCI <- matrix(unlist(LowerBoundCI), nrow = length(LowerBoundCI), byrow = TRUE)
        UpperBoundCI <- matrix(unlist(UpperBoundCI), nrow = length(UpperBoundCI), byrow = TRUE)
        
        
        re <- bind_cols(p, eta, x.p, SE, SEup, SElo,  se_pi, LowerBoundCI,UpperBoundCI)
        colnames(re) <- c("Fitness","eta","TreatmentTemp","SE", "SEup", "SElo","se_pi" ,"LowerBoundCI", "UpperBoundCI") ## name column
        re <- as.data.frame(re)
        return(re)
        
      } else {
        x.p <- (qnorm(p) - b1 -  b3) / b2
        eta <- b1 + (x.p * b2) + b3
        pd <- -cbind(1, x.p)/b[2L]
        cf2 <- c(1,2)
        SE <- sqrt(((pd %*% vcov(obj)[cf2, cf2]) * pd) %*% c(1, 1))
        z <- qnorm(0.975)
        SEup <- pnorm(eta - z * SE)
        SElo <- pnorm(eta - z * SE)
        
        # Calculate the standard error
        se_values <- sapply(x.p, calculate_se_single, obj = obj, cf = cf)
        se_pi <- se_values[4,]
        
        
        LowerBoundCI <- pnorm(eta - z * se_pi)
        UpperBoundCI <- pnorm(eta + z * se_pi)
        
        se_pi <- matrix(unlist(se_pi), nrow = length(se_pi), byrow = TRUE)
        LowerBoundCI <- matrix(unlist(LowerBoundCI), nrow = length(LowerBoundCI), byrow = TRUE)
        UpperBoundCI <- matrix(unlist(UpperBoundCI), nrow = length(UpperBoundCI), byrow = TRUE)
        
        
        re <- bind_cols(p, eta, x.p, SE, SEup, SElo,  se_pi, LowerBoundCI,UpperBoundCI)
        colnames(re) <- c("Fitness","eta","TreatmentTemp","SE", "SEup", "SElo","se_pi" ,"LowerBoundCI", "UpperBoundCI") ## name column
        re <- as.data.frame(re)
        return(re)
        
        
      }
      
      
    }
    
      
      
    }
 
  
}
  if(is.null(p)){
    b <- coef(obj)[cf]
    eta <- b[1L] + tp * b[2L]
    x.p <- tp
    pd <- -cbind(1, x.p)/b[2L]
    SE <- sqrt(((pd %*% vcov(obj)[cf, cf]) * pd) %*% c(1, 1))
    z <- qnorm(0.975)
    SEup <- pnorm(eta - z * SE)
    SElo <- pnorm(eta - z * SE)
    p <- pnorm(eta)
    
    
    # Obtain the variance-covariance matrix of the model coefficients
    vcov_matrix <- vcov(obj)
    
    # Calculate the standard error
    se_values <- sapply(x.p, calculate_se_single, obj = obj, cf = cf)
    se_pi <- se_values[4,]
    
    LowerBoundCI <- pnorm(eta - z * se_pi)
    UpperBoundCI <- pnorm(eta + z * se_pi)
    
    
    se_pi <- matrix(unlist(se_pi), nrow = length(se_pi), byrow = TRUE)
    LowerBoundCI <- matrix(unlist(LowerBoundCI), nrow = length(LowerBoundCI), byrow = TRUE)
    UpperBoundCI <- matrix(unlist(UpperBoundCI), nrow = length(UpperBoundCI), byrow = TRUE)
    
    
    re <- bind_cols(p, eta, tp, SE, SEup, SElo, se_pi, LowerBoundCI,UpperBoundCI)
    colnames(re) <- c("Fitness","eta","TreatmentTemp","SE", "SEup", "SElo","se_pi" ,"LowerBoundCI", "UpperBoundCI") ## name column
    re <- as.data.frame(re)
    
    return(re)
    
    
  }
  
}


