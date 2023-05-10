library(here)
library(dplyr)
library(rcarbon)

# Read dates ----
burial  <- read.csv(here('data','raw','All_burialdates_0620.csv'))|> select(SiteID,Country,CaseStudyArea,LabID,Material,CRA,Error,Longitude,Latitude,Disarticulated,Cremation) 
burial$Cremation  <- ifelse(burial$Cremation=='Yes',1,0)
burial$Disarticulated  <- ifelse(burial$Disarticulated=='Yes',1,0)

save(burial,file=here('data','burialdata.RData'))
