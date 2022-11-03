#Laura F. White
#Boston University
#Pre-chemotherapy TB Analysis

##############################################################################
# This program generates the plots used to show all the cure data reported in 
# the supplemental materials
##############################################################################

options(scipen=999)
options(digits = 10)

rm(list = ls())
source("R/utils.R")
reload_source()

# Reading in the study_id correspondence table
studyid <- read.csv("data/study_id.csv")

# Reading in cure data
cureData <- read.csv("data/cure_data_all.csv")

# add in author and date of study
cureData2 <- inner_join(cureData,studyid,by="study_id")

# make clear that positive and negative from Baart de la Faille is smear
cureData2$severity2 <- cureData2$severity
cureData2$severity2[cureData2$severity=="Positive"] <- "Smear Positive"
cureData2$severity2[cureData2$severity=="Negative"] <- "Smear Negative"

ggplot(cureData2[cureData2$severity2!="None",],
       aes(x=interval_r,y=cureRate,group=cohort_id,color=severity2))+
  geom_point(aes(shape=first_author))+geom_line()+xlim(0,10)+ylim(0,1)+
  labs(title="Studies with severity reported",
       x="Year since diagnosis",y="Probability of Natural Recovery",color="Severity",shape="Author")+
  scale_color_brewer(palette="Dark2")

# cure data by sanatorium versus not
ggplot(cureData2,
       aes(x=interval_r,y=cureRate,group=cohort_id,color=sanatorium.x))+
  geom_point()+geom_line()+xlim(0,10)+ylim(0,1)+
  labs(title="Natural recovery for Sanatorium and non-Sanatorium Studies",
       x="Year since diagnosis",y="Probability of Natural Recovery",color="Sanatorium")+
  scale_color_brewer(palette="Dark2")

# cure data by time period
ggplot(cureData2,
       aes(x=interval_r,y=cureRate,group=cohort_id,color=time_period))+
  geom_point()+geom_line()+xlim(0,10)+ylim(0,1)+
  labs(title="Natural recovery by time period",
       x="Year since diagnosis",y="Probability of Natural Recovery",color="Time Period")+
  scale_color_brewer(palette="Dark2")

# cure data by location
ggplot(cureData2,
       aes(x=interval_r,y=cureRate,group=cohort_id,color=location))+
  geom_point()+geom_line()+xlim(0,10)+ylim(0,1)+
  labs(title="Natural recovery by location",
       x="Year since diagnosis",y="Probability of Natural Recovery",color="Location")+
  scale_color_brewer(palette="Dark2")

# combine together in one plot
ggplot(cureData2,
       aes(x=interval_r,y=cureRate,group=cohort_id,color=sanatorium.x))+
  geom_point(aes(shape=location))+
  geom_line(aes(linetype=time_period))+xlim(0,10)+ylim(0,1)+
  labs(title="Natural recovery data",
       x="Year since diagnosis",y="Probability of Natural Recovery",color="Sanatorium",linetype="Time Period",
       shape="Location")+
  scale_color_brewer(palette="Dark2")



### Create a data frame with all the data -------------------------------------------
cureData3 <- cureData2[,c('study_id','severity','interval_r','n','cureRate','first_author','location','time_period','sanatorium.x','year')]

# round interval_r so we are looking at cureRate within a one year interval of the integer time
# also remove interval data that is greater than 10
cureData3$interval_r <- ceiling(cureData3$interval_r)

cureData4 <- subset(cureData3,cureData3$interval_r<=10 & cureData3$interval_r>0)
cureSummary <- spread(cureData4,key="interval_r",value="cureRate")

write.csv(cureSummary,file="data/cureDataSummary.csv")


