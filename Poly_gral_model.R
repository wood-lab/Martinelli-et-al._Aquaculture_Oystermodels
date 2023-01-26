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
library(wesanderson)

###########################################################
## ADDING DATASETS
###########################################################
summary_all <- read.table("sampling summary.csv", header=T,sep=",")

prevalence <- read.table("master_spreadsheet_all.csv", header=T,sep=",")
prevalence <- subset(prevalence, prevalence$Valve =='R') # only keep R valves
# dim: 4085, 19 

#########################
#PLOTTING PREV PER STATE

level_order <- c('CA', 'OR', 'WA', 'AK') 

prev_state <- ggplot(summary_all, aes(x= factor(State, level= level_order), y= Prevalence, fill= State, show.legend = FALSE))  +
        stat_boxplot(geom= 'errorbar' , width = 0.3, position = position_dodge(width = 0.75) ) +
        geom_boxplot(alpha=0.7, lwd=1, outlier.shape = NA, show.legend = FALSE) + 
        geom_point(stat = "identity", size= 4, shape = 21, lwd=2, show.legend = FALSE) + 
        ylim(0,1) +
        scale_fill_manual(values=wes_palette("GrandBudapest1", n = 4)) + 
        labs(x = 'State', y = 'Prevalence', size=16) 
prev_state + theme_classic(base_size = 18) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

data_mean <- summary_all %>%
        group_by(State) %>% #
        summarize_at(vars(Prevalence), list(mean = mean), na.rm=TRUE)

data_median <- summary_all %>%
        group_by(State) %>% #
        summarize_at(vars(Prevalence), list(median = median), na.rm=TRUE)

data_state <- summary_all %>%
        group_by(State) %>% #
        summarize_at(vars(Total.oysters), list(sum = sum), na.rm=TRUE)

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
model1 <- glmer(Infested ~ y + Culture + Ploidy + Season + Thick_sc + (1|Year_sc) + (1|State/Farm), family="binomial", data = prevalence)
summary(model1)
anova(model1) 
vif(model1)
car::Anova(model1, type=3) # getting p-values 


model2 <- glmer(Infested ~ y_sc + Season*State + Culture*State + (1|Year_sc) + (1|Farm), family="binomial", data = prevalence)
summary(model2)
anova(model2) 
vif(model2) 
car::Anova(model2, type=3) # getting p-values 

## PLOTTING MODEL 1
########################
allstates <- ggpredict(model1,c("y","Season"))

allstates_plot<-ggplot(allstates,aes(x,predicted,color=group), color=group) +
        scale_color_manual(values=wes_palette("GrandBudapest1", n = 2)) + 
        geom_point(size=4) +
        geom_errorbar(data=allstates, mapping=aes(x=x, ymin=conf.low, ymax=conf.high), width=0.03) +
        geom_line(aes(group=group)) +
        xlab("Shell height (cm)") +
        ylab(expression(paste("Predicted infestation"))) +
        theme_classic() +
        guides(color=guide_legend("Season")) +
        theme(plot.title=element_text(size=14,hjust=0.5,face="plain"), axis.text.y=element_text(size=14), 
        axis.title.y=element_text(size=14), axis.text.x=element_text(size=14), axis.title.x=element_text(size=14),
        panel.grid.minor=element_line(color=NA))
allstates_plot

## PLOTTING MODEL 2
########################
summary(model2)

# allstates_plot2 <- plot(allstates2) +
#         scale_color_manual(values=wes_palette("GrandBudapest1", n = 4)) + 
#         ylim(0,1) +
#         #scale_x_discrete(limits = c("CA","OR","WA", "AK")) +
#         labs(x = 'State', y = 'Predicted infestation', size=16)
# 
# allstates_plot2 +  theme_classic(base_size = 18) + theme(plot.title=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5))

allstates2 <- ggpredict(model2,c("State", "Culture"))
mod2_dd <- as.data.frame(allstates2)
mod2_dd$culture <- paste(mod2_dd$x, mod2_dd$group, sep=" ")
level_order <- c('CA on','CA off','OR on','OR off','WA on','WA off','AK on','AK off') 


Fig5B <- ggplot(mod2_dd, aes(x= factor(culture, level= level_order), y= predicted, fill= x, show.legend = FALSE))  +
        stat_boxplot(geom= 'errorbar', width = 0.3, position = position_dodge(width = 0.75)) +
        geom_errorbar(data=mod2_dd, mapping=aes(x=culture, ymin=conf.low, ymax=conf.high), width=0.03) +
        geom_boxplot(alpha=0.7, lwd=1, outlier.shape = NA, show.legend = FALSE) + 
        geom_point(stat = "identity", size= 4, shape = 21, lwd=2, show.legend = FALSE) + 
        ylim(0,1) +
        scale_fill_manual(values=wes_palette("GrandBudapest1", n = 4)) + 
        labs(x = 'State', y = 'Predicted Infestation', size=16) 
Fig5B + theme_classic(base_size = 18) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


allstates3 <- ggpredict(model2,c("State","Season"))
allstates3_dd <- as.data.frame(allstates3)
allstates3_dd$season <- paste(allstates3_dd$x, allstates3_dd$group, sep=" ")
level_order <- c('CA Summer','CA Winter','OR Summer','OR Winter','WA Summer','WA Winter','AK Summer','AK Winter') 


Fig5A <- ggplot(allstates3_dd, aes(x= factor(season, level= level_order), y= predicted, fill= x, show.legend = FALSE))  +
        stat_boxplot(geom= 'errorbar', width = 0.3, position = position_dodge(width = 0.75)) +
        geom_errorbar(data=allstates3_dd, mapping=aes(x=season, ymin=conf.low, ymax=conf.high), width=0.03) +
        geom_boxplot(alpha=0.7, lwd=1, outlier.shape = NA, show.legend = FALSE) + 
        geom_point(stat = "identity", size= 4, shape = 21, lwd=2, show.legend = FALSE) + 
        ylim(0,1) +
        scale_fill_manual(values=wes_palette("GrandBudapest1", n = 4)) + 
        labs(x = 'State', y = 'Predicted Infestation', size=16) 
Fig5A + theme_classic(base_size = 18) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


allstates_plot3 <- plot(allstates3) +
        scale_color_manual(values=wes_palette("GrandBudapest1", n = 2)) + 
        ylim(0,1) +
        labs(x = 'State', y = 'Predicted infestation', size=16)

allstates_plot3 +  theme_classic(base_size = 18) + theme(plot.title=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5))



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








