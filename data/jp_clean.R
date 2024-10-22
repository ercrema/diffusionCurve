# Load Library ----
library(here)
library(rcarbon)
library(dplyr)

# Load Data ----
# Read CSV data from https://www.rekihaku.ac.jp/up-cgi/login.pl?p=param/esrd_en/db_param
c14data  <- read.csv(here('data','raw','c14db_1.1.0.csv'))

# Prioritise unrounded dates ----
c14data$C14Age = as.numeric(c14data$UnroundedCRA)
c14data$C14Error = as.numeric(c14data$UnroundedCRAError)
i = which(is.na(c14data$UnroundedCRA))
j = which(is.na(c14data$UnroundedCRAError))
c14data$C14Age[i] = c14data$CRA[i]
c14data$C14Error[j] = c14data$CRAError[j]

# Window of analyses 4000 to 1000 14C age
c14data  <- subset(c14data,C14Age<=4000 & C14Age >=1000)

# Consider only seeds
seeds  <- subset(c14data,MaterialDetails %in% c('Charred Seed','Seed')) |> select(LabCode,MaterialDetails,Taxa=MaterialTaxa,Prefecture,Region=Region,SiteNameJp,Latitude,Longitude,C14Age,C14Error)


# Read rice data and pool samples TKA-23237 amd TKA-23238
rice  <- read.csv(here('data','raw','R14CDB.csv')) |> subset(UseForAnalyses == 'TRUE')
poolLabCodes  <- c('TKA-23237','TKA-23238')
i  <- which(rice$LabCode%in%poolLabCodes)
pooledDates  <- poolDates(rice$C14Age[i],rice$C14Error[i])[2:3]
rice  <- rbind.data.frame(rice,rice[i[1],])
rice$ID[nrow(rice)] <- max(rice$ID) + 1
rice$LabCode[nrow(rice)] <- 'Combined_TKA23237_TKA23238'
rice$C14Age[nrow(rice)] <- as.numeric(pooledDates[1])
rice$C14Error[nrow(rice)] <- as.numeric(pooledDates[2])
rice$Context[nrow(rice)] <- 'Charred rice grain embedded in sherd fabric (Combined Dates)'
rice  <- rice[-i,]
any(seeds$LabCode == 'TKA-23237')
any(seeds$LabCode == 'TKA-23238')

# Cross reference the two tables
seeds$rice = FALSE #initially assume seeds are all not rice.
seeds$rice[which(seeds$LabCode%in%rice$LabCode)] = TRUE #labcodes present in the rice database
i  <- which(!rice$LabCode%in%seeds$LabCode)
rice  <-  rice[i,]
rice$rice  <- TRUE
rice$MaterialDetails  <- 'Charred Seed'
rice$Taxa  <- 'Oryza sativa'


rice  <- select(rice,LabCode,MaterialDetails,Taxa,Prefecture,Region,SiteNameJp=SiteName_jp,Latitude,Longitude,C14Age,C14Error,rice)
jpdata  <- rbind.data.frame(seeds,rice) |> subset(Region != 'Hokkaido')

taxa_classification  <- read.csv(here('data','raw','Taxa_Edible_Classifications.csv'))
jpdata <- left_join(jpdata,taxa_classification)

jpdata$rice_nuts  <- NA
jpdata$rice_nuts[which(jpdata$rice==TRUE)]  <- 1
jpdata$rice_nuts[which(jpdata$Classification=='nuts (tree)')]  <- 0
jpdata  <- subset(jpdata,!is.na(rice_nuts))

# Subset to chronological window of analyses
xx <- calibrate(jpdata$C14Age,jpdata$C14Error)
jpdata <- jpdata[which.CalDates(xx,BP<4000&BP>1700,p=0.5),]

# Summaries
nrow(jpdata) #551
length(unique(jpdata$SiteNameJp))

# Store Output ----
save(jpdata,file=here('data','jpdata.RData'))
