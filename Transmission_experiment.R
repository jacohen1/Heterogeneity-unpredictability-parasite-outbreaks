#### INITIAL CODE ####
rm(list=ls())

#load in packages
library(here)
library(tidyverse)
library(reshape)
library(MASS)
library(emmeans)
library(DHARMa)
library(lme4)

#read in data
full <- read.csv(here("transmission_experiment.csv"), header = TRUE)
#for the treatment column, the letters correspond to treatments as follows:
#A -- 100% Colony A
#B -- 75% Colony A
#C -- 50% Colony A
#D -- 25% Colony A
#E -- 0% Colony A

#### DATA RE-STRUCTURING ####
#adjust data structure
full <- full %>%
  mutate(treat_rep = paste0(treatment, replicate_no)) %>% #add column combining treatment number and replicate number
  mutate(across(c(treatment, replicate_no, dissection_no, treat_rep), as.factor))

#remove dead, pupated, bad dissections
full <- full %>%
  filter(!trophont %in% c("dead", "bad_dissection", "pupated")) %>%
  mutate(across(c(trophont, gamont, gut_gametocyst, par_presence), as.integer)) %>% #turn columns into integers, rather than characters
  mutate(no_gut_parasites = trophont + (2*gamont) + (2*gut_gametocyst)) %>% #calculate total number of gut parasites
  mutate(no_gut_parasites = as.integer(no_gut_parasites))

#get parasite intensity data (only individuals with > 0 parasites)
int <- full %>%
  filter(no_gut_parasites > 0)

#### PREVALENCE GLMM ####
prev.mod <- glmer(par_presence ~ dissection_no + treatment + (1|treat_rep), data = full, family = "binomial")
summary(prev.mod)
#AIC 2335.3

#check fit
prev.output <- simulateResiduals(fittedModel = prev.mod)
plot(prev.output)

#post hoc test for differences in treatment
prev.posthoc <- emmeans(prev.mod, specs = pairwise ~ treatment, weights = "proportional")
summary(prev.posthoc)
#significant differences between A&D, A&E

#get estimate for mean prevalence per treatment
prev.means <- as.data.frame(summary(prev.posthoc$emmeans, type = "response"))

#### INTENSITY GLMM ####
int.mod <- glmer.nb(no_gut_parasites ~ dissection_no + treatment + (1|treat_rep), data = int)
summary(int.mod)
#AIC 13456.1

#check fit
int.output <- simulateResiduals(fittedModel = int.mod)
plot(int.output)

#post hoc test for differences in treatment
int.posthoc <- emmeans(int.mod, specs = pairwise ~ treatment, weights = "proportional")
summary(int.posthoc)
#significant differences between A and all other treatments, B&D, C&D

#get estimate for mean intensity per treatment
int.means <- as.data.frame(summary(int.posthoc$emmeans, type = "response"))

#### PREVALENCE - GENERATE EXPECTED VALUES ####
#create blank matrix
prev_vals <- matrix(nrow = 10000, ncol = 3)

#get model outputs for prevalence values
a.prev <- prev.means$prob[prev.means$treatment == "A"]
e.prev <- prev.means$prob[prev.means$treatment == "E"]

#get post hoc prevalences with sample sizes
prev.posthoc2 <- emmeans(prev.mod, specs = ~ treatment, weights = "proportional")
#add number of samples for each treatment to dataframe
prev.means$samples <- prev.posthoc2@grid[[".wgt."]]

for(j in 1:10000){
  
  #generate a binomial distribution with a success probability from the estimated per treatment prevalence
  a.dist <- rbinom(prev.means$samples[prev.means$treatment == "A"], size = 1, prob = a.prev)
  e.dist <- rbinom(prev.means$samples[prev.means$treatment == "E"], size = 1, prob = e.prev)
  
  #get mean prev_vals for distributions
  A <- mean(a.dist)
  E <- mean(e.dist)
  
  #generate 25%, 50%, 75% distances between the two points
  B <- (0.25 * (E - A)) + A
  C <- (0.5 * (E - A)) + A
  D <- (0.75 * (E - A)) + A
  
  prev_vals[j,] <- c(B, C, D)
}

#create blank matrix to input summary results of the bootstrapping process
prev_pred <- matrix(nrow = 3, ncol = 5)
#add mean and sd of these prev_vals to the original blank matrix
prev_pred[1,] <- c("B", mean(prev_vals[,1]), sd(prev_vals[,1]), NA, NA)
prev_pred[2,] <- c("C", mean(prev_vals[,2]), sd(prev_vals[,2]), NA, NA)
prev_pred[3,] <- c("D", mean(prev_vals[,3]), sd(prev_vals[,3]), NA, NA)

prev_pred <- as.data.frame(prev_pred)

for(i in 1:nrow(prev_pred)){
  m <- as.numeric(prev_pred[i, 2])
  st.dev <- as.numeric(prev_pred[i, 3])
  high <- m + 1.96*st.dev
  low <- m - 1.96*st.dev
  if(high < 1){
    conf.high <- high
  } else{
    conf.high <- 1
  }
  if(low > 0){
    conf.low <- low
  } else{
    conf.low <- 0
  }
  
  prev_pred[i, 4] <- conf.low
  prev_pred[i, 5] <- conf.high
}

colnames(prev_pred) <- c("treatment", "mean", "sd", "conf.low", "conf.high")

#### INTENSITY - GENERATE EXPECTED VALUES ####
#get dispersion parameter from intensity model
theta <- lme4:::getNBdisp(int.mod)

#create blank matrix
int_vals <- matrix(nrow = 10000, ncol = 3)

#get model outputs for intensity int_vals
a.int <- int.means$response[int.means$treatment == "A"]
e.int <- int.means$response[int.means$treatment == "E"]

#get post hoc intensity int_vals with sample sizes
int.posthoc2 <- emmeans(int.mod, specs = ~ treatment, weights = "proportional")
#add number of samples for each treatment to dataframe
int.means$samples <- int.posthoc2@grid[[".wgt."]]

for(j in 1:10000){
  
  #generate a negative binomial distribution with mu as the model intensity output theta the model dispersion parameter
  a.dist <- rnbinom(int.means$samples[int.means$treatment == "A"], size = theta, mu = a.int)
  e.dist <- rnbinom(int.means$samples[int.means$treatment == "E"], size = theta, mu = e.int)
  
  #get mean int_vals for new distributions
  A <- mean(a.dist)
  E <- mean(e.dist)
  
  #generate 25%, 50%, 75% distances between the two points
  B <- (0.25 * (E - A)) + A
  C <- (0.5 * (E - A)) + A
  D <- (0.75 * (E - A)) + A
  
  int_vals[j,] <- c(B, C, D)
  
}

#create blank matrix
int_preds <- matrix(nrow = 3, ncol = 5)
#add mean and sd of these int_vals to the original blank matrix
int_preds[1,] <- c("B", mean(int_vals[,1]), sd(int_vals[,1]), NA, NA)
int_preds[2,] <- c("C", mean(int_vals[,2]), sd(int_vals[,2]), NA, NA)
int_preds[3,] <- c("D", mean(int_vals[,3]), sd(int_vals[,3]), NA, NA)

int_preds <- as.data.frame(int_preds)

for(i in 1:nrow(int_preds)){
  m <- as.numeric(int_preds[i, 2])
  st.dev <- as.numeric(int_preds[i, 3])
  conf.high <- m + 1.96*st.dev
  low <- m - 1.96*st.dev
  if(low > 0){
    conf.low <- low
  } else{
    conf.low <- 0
  }
  
  int_preds[i, 4] <- conf.low
  int_preds[i, 5] <- conf.high
}

colnames(int_preds) <- c("treatment", "mean", "sd", "conf.low", "conf.high")
