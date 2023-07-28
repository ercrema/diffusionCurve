library(here)
library(dplyr)
library(rcarbon)

# Read dates ----
burial  <- read.csv(here('data','raw','burialdates.csv'))|> select(SiteID,Country,CaseStudyArea,LabID,Material,CRA,Error,Longitude,Latitude,Disarticulated,Cremation) 
burial$Cremation  <- ifelse(burial$Cremation=='Yes',1,0)
burial$Disarticulated  <- ifelse(burial$Disarticulated=='Yes',1,0)
cal  <- calibrate(burial$CRA,burial$Error)
burial  <- burial[which.CalDates(cal,BP<5500 & BP>2000,p=0.99),]
nrow(burial)
table(burial$Cremation)

save(burial,file=here('data','burialdata.RData'))
