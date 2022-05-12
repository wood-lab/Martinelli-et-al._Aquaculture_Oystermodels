### WORKING DIRECTORY AND PACKAGES
###########################################################
library(tidyverse)
library(tidyr)
library(plyr)
library(corrplot) 
library(zoo)
library(lme4) 
library(car)
library(lubridate)

###########################################################
## ADDING PREVALENCE DATASET
###########################################################
setwd('/Users/jmartine/Desktop/buoys/')
prevalence <- read.csv('master_spreadsheet_all.csv', header=TRUE, row.names = NULL, stringsAsFactors=FALSE,  sep=',') # na.strings=c('NA','\\N')
# remove L valves to avoid duplication of weights etc
prevalence <- subset(prevalence, prevalence$Valve =='R') # only keep R valves
# dim: 4085, 19 

########################
# CHECKING FOR CORRELATION IN oyster metrics
prevalence$Shell_g <- as.numeric(prevalence$Shell_g)     
prevalence$Tissue_g <- as.numeric(prevalence$Tissue_g)
prevalence$Thick <- as.numeric(prevalence$Thick)
prevalence$x <- as.numeric(prevalence$x)
prevalence$y <- as.numeric(prevalence$y)
prevalence$z <- as.numeric(prevalence$z)

prevalence.mat <- prevalence[, c(6,13:18)] # selecting the cols with 
matrix <- cor(prevalence.mat, use="pairwise.complete.obs")
par(xpd = TRUE)
corrplot(matrix,method="number",type="lower",tl.col="black",tl.srt=45,mar=c(0,0,0,0))

#### CREATING SEASON
yq <- as.yearqtr(as.yearmon(prevalence$Date, "%m/%d/%Y") + 1/12)
prevalence$Season <- factor(format(yq, "%q"), levels = 1:4, 
                            labels = c("Winter", "Winter", "Summer", "Summer"))
prevalence <- prevalence %>% mutate(Date1 = Date)
dmy(prevalence$Date)
## FIX HERE!



prevalence <- prevalence %>% separate(Date1, sep="/", into = c('Month', 'Day', 'Year'))

#### SCALING each column individually
###########################################################
prevalence$Lat_sc <- scale(prevalence$Lat, center = TRUE, scale = TRUE)
prevalence$Lng_sc <- scale(prevalence$Lng, center = TRUE, scale = TRUE)
prevalence$Shell_g_sc <- scale(prevalence$Shell_g, center = TRUE, scale = TRUE)
prevalence$Tissue_g_sc <- scale(prevalence$Tissue_g, center = TRUE, scale = TRUE)
prevalence$x_sc <- scale(prevalence$x, center = TRUE, scale = TRUE)
prevalence$y_sc <- scale(prevalence$y, center = TRUE, scale = TRUE)
prevalence$z_sc <- scale(prevalence$z, center = TRUE, scale = TRUE)
prevalence$Thick_sc <- scale(prevalence$Thick, center = TRUE, scale = TRUE)

## MODEL TESTING FOR PREVALENCE - ALL STATES
########################
model1 <- glmer(Infested ~ Tissue_g_sc + y_sc + Season + (1|Date) + Thick_sc + Culture + Ploidy + (1|State/Bay/Farm), family="binomial", data = prevalence)
summary(model1)
anova(model1)
vif(model1)

# try random effect for Date, create that one again (1|Date)
# explore interstate variab between states w/interactions:
# interaction w/ state and season
# interaction w/ state and culture
# plots! ggpredict, ggeffects


# FAILS TO CONVERGE WITH LAT + (1|Lat), its either bay/farm or lat
# The intercept is the predicted value of the dependent variable when the independent variables are 0
# estimate of intercept, -2.363988, is the log odd of infested being 1
# a 1 unit increase in culture on is associate with a 0.44 increase in the log odd of infested being 1

########################
## WASHINGTON # glm as there's 9 single-farm bay, out of total 12 bays
# CHECK PLOIDY WA - where are diploids??!!
wa <- subset(prevalence, prevalence$State =='WA')
modelwa <- glmer(Infested ~ Tissue_g_sc + y_sc + Thick_sc + Culture + Season + (1|Bay/Farm), family="binomial", data = wa)
summary(modelwa) 

## CALIFORNIA
ca <- subset(prevalence, prevalence$State =='CA')
modelca <- glmer(Infested ~ Tissue_g_sc + y_sc + Thick_sc + Culture + Season + Ploidy + (1|Bay/Farm), family="binomial", data = ca)
summary(modelca)

## OREGON (no ploidy data for OR)
or <- subset(prevalence, prevalence$State =='OR')
modelor <- glmer(Infested ~ Tissue_g_sc + y_sc + Thick_sc + Culture + Season + (1|Bay/Farm), family="binomial", data = or)
summary(modelor)

## ALASKA
# glm as there's 4 single-farm bay, out of total 5 bays
ak <- subset(prevalence, prevalence$State =='AK') # only dips in AK
modelak <- glmer(Infested ~ Tissue_g_sc + y_sc + Thick_sc + Culture + Season + (1|Farm), family="binomial", data = ak)
summary(modelak)


