### WORKING DIRECTORY AND PACKAGES
###########################################################
library(plyr)
library(tidyverse)
library(corrplot) 
library(zoo)
library(lme4) 
library(car)
library(lubridate)
library(devtools)
library(ggeffects)
library(statmod)

###########################################################
## ADDING PREVALENCE DATASET
###########################################################
prevalence <- read.table("master_spreadsheet_all.csv", header=T,sep=",")
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
prevalence$Date1 <- mdy(prevalence$Date1)
prevalence <- prevalence %>% separate(Date1, sep="-", into = c('Year','Month', 'Day'))
prevalence$Date <- mdy(prevalence$Date)
prevalence$Year <- as.numeric(prevalence$Year)

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
prevalence$Date_sc <- scale(prevalence$Date, center = TRUE, scale = TRUE)
prevalence$Year_sc <- scale(prevalence$Year, center = TRUE, scale = TRUE)

## MODEL TESTING FOR PREVALENCE - ALL STATES
########################
model1 <- glmer(Infested ~ y_sc + Culture + Ploidy + Season + Thick_sc + (1|Year_sc) + (1|State/Farm), family="binomial", data = prevalence)
summary(model1)
anova(model1)
vif(model1)

model2 <- glmer(Infested ~ y_sc + Ploidy + Season*State + Culture*State + (1|Year_sc) + (1|Farm), family="binomial", data = prevalence)
summary(model2)
anova(model2)
vif(model2)

## trying out a plot
allstates <- ggpredict(model1,c("y_sc","Season"))

allstates_plot<-ggplot(allstates,aes(x,predicted,color=group), color=group)+
        scale_color_manual(values=c("blue","turquoise"))+ #F03B20 #7FCDBB
        geom_point(size=4)+
        geom_errorbar(data=allstates,mapping=aes(x=x,ymin=conf.low,ymax=conf.high),width=0.03)+
        geom_line(aes(group=group))+
        xlab("Shell height")+
        ylab(expression(paste("Predicted infestation")))+
        theme_classic()+
        theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.text.x=element_text(size=14),axis.title.x=element_text(size=14),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey90"),panel.grid.minor=element_line(color=NA))
        #scale_x_discrete(limits=factor(y_reversed))
allstates_plot


#using the attributes to unscale & plot right values on x axis
scaled_y <- as.data.frame(allstates$x)
y_reversed <- (scaled_y * attr(scaled_y, 'scaled:scale')) + attr(scaled_y, 'scaled:center')

## higher infestation in winter, higher pH in winter (pH decreases with higher T)
# modified state/bay/season and got rid of bay
# ggpredict, ggeffects
# The intercept is the predicted value of the dependent variable when the independent variables are 0
# estimate of intercept, -2.363988, is the log odd of infested being 1
# a 1 unit increase in culture on is associate with a 0.44 increase in the log odd of infested being 1


######################## 
## WASHINGTON # similar number bay to farm
wa <- subset(prevalence, prevalence$State =='WA')
modelwa <- glmer(Infested ~ y_sc + Culture + Season + Ploidy + (1|Year_sc) + (1|Bay), family="binomial", data = wa)
summary(modelwa)

## CALIFORNIA
ca <- subset(prevalence, prevalence$State =='CA')
modelca <- glmer(Infested ~ y_sc + Culture + Season + Ploidy + (1|Year_sc) + (1|Bay/Farm), family="binomial", data = ca)
summary(modelca)

## OREGON 
or <- subset(prevalence, prevalence$State =='OR')
modelor <- glmer(Infested ~ y_sc + Culture + Season + Ploidy + (1|Year_sc) + (1|Bay), family="binomial", data = or)
summary(modelor)

## ALASKA
# glm as there's 4 single-farm bay, out of total 5 bays (1|Year)
ak <- subset(prevalence, prevalence$State =='AK') # only dips in AK
modelak <- glmer(Infested ~ y_sc + Culture + Season + (1|Year_sc) + (1|Bay), family="binomial", data = ak)
summary(modelak)








