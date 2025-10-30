#### INITIAL CODE ####
rm(list=ls())

#load in packages
library(here)
library(tidyverse)
library(MASS)
library(emmeans)
library(DHARMa)
library(lme4)

#read in data #CHANGE ALL THIS TO HERE() SYNTAX
gam_ooc <- read.csv(here("gametocyst_oocyst.csv"), header = TRUE)
dev <- read.csv(here("development.csv"), header = TRUE)

#### GAMETOCYST/OOCYST - DATA RE-STRUCTURING ####
gam_ooc <- gam_ooc %>%
  filter(!oocysts %in% c("dead", "bad_dissection", "pupated")) %>%
  mutate(across(c(oocysts, gametocysts, par_presence), as.integer)) %>% #turn columns into integers, rather than characters
  mutate(across(c(colony, experimental_repeat, day_removed), as.factor)) %>% #other columns into factors for GLMM
  mutate(oocysts = as.integer(400*oocysts)) %>% #multiply oocyst count by 400 to get number of oocysts in entire sample
  filter(oocysts > 0) #only keep individuals who produced oocysts

#### GAMETOCYST/OOCYST - GLMM ####
gam_mod <- glmer(oocysts ~ gametocysts + colony + experimental_repeat + (1|day_removed) + colony*gametocysts, data = gam_ooc,
                 family = Gamma(link = "log"))
summary(gam_mod)
#AIC 1002

#check fit
gam_output <- simulateResiduals(fittedModel = gam_mod)
plot(gam_output)

#posthoc test for overall difference between colonies
gam_posthoc <- emmeans(gam_mod, specs = pairwise ~ colony, weights = "proportional")
summary(gam_posthoc)
#no significant difference between colonies

#### DEVELOPMENT - DATA RE-STRUCTURING####
dev <- dev %>%
  mutate(across(c(colony, group), as.factor))

#### DEVELOPMENT - GLM ####
dev_mod <- glm(development ~ day + day*group + colony + colony*group + experimental_repeat, data = dev, family = "binomial")
summary(dev_mod)
#AIC 1589 

#check fit
dev_output <- simulateResiduals(fittedModel = dev_mod)
plot(dev_output)

#posthoc test for differences in development between exposed and unexposed larvae
dev_posthoc <- emmeans(dev_mod, specs = pairwise ~ group, weights = "proportional")
summary(dev_posthoc)
#development rate faster in unexposed than exposed
