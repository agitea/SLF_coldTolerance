## SLF Cold tolerance and lower lethal limits analysis
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


## data file -------------------------------------------------------------------
SCP_file <- "data/SLF_2022_Y1_SCP.csv" ## Supercooling points of L. delicatula eggs data


################################################################################
# Supercooling points - Fig 2 and analysis
################################################################################

## read data -------------------------------------------------------------------
SCP_df <- read.csv(paste0(SCP_file), header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE)

#
## reshape data ----------------------------------------------------------------
SCP_df$EggMassN <- dplyr::ntile(SCP_df$EggMass, 5)
SCP_df <- SCP_df %>% dplyr::filter(!is.na(SCP_df$EggMass))
SCP_df$colorgp <- as.character(SCP_df$EggMass)
SCP_df$x <- 1


## plot Fig.2 ------------------------------------------------------------------
SCP_df %>% ggplot(aes(x=x, y=SCP), alpha = 0.7) +
  geom_boxplot(notch = FALSE, notchwidth = 0.55,width=0.8, position = position_dodge(width=0.1)) +
  scale_color_tableau(palette = "Color Blind")+
  geom_jitter(aes(color = colorgp), size=2, alpha=0.9) +
  theme_ipsum() +
  NULL


## Kruskal-Wallis test ---------------------------------------------------------

k_test <- kruskal.test(SCP ~ EggMass, data = SCP_df) 
k_test


























