# Load Library ----
library(here)
library(rcarbon)
library(sf)
library(rnaturalearth)

# Read Data ----
# Download data from Bevan et al 2017
temp <- tempfile()
download.file("https://discovery.ucl.ac.uk/id/eprint/10025178/4/Bevan_gbie14Csub.zip",temp)
gbdata <- read.csv(unz(temp, "gbie14Csub/dates/archdates.csv"))


# Subset dates for Wheat, Barley, and Hazelnuts ----
gbdata$cat  <- NA
toMatchMat <- c("nutshell","grain","fruit","seed")
check1 <- grepl(paste(toMatchMat,collapse="|"), gbdata$Material)
gbdata$cat[which(grepl("Corylus",gbdata$Species) & check1)] = 'hazelnut'

check <- grepl("Triticum",gbdata$Species)
toExclude <- c("Hordeum/Triticum","Triticum/Hordeum vulgare ", "Hordeum vulgare/Triticum spelta", "Hordeum vulgare/Triticum","Avena/Triticum")
gbdata$cat[which(check & !gbdata$Species %in% toExclude)] = 'wheat'

check <- grepl("Hordeum",gbdata$Species)
toExclude <- c("Hordeum/Triticum","Triticum/Hordeum vulgare ", "Hordeum vulgare/Triticum spelta", "Hordeum vulgare/Triticum","Avena/Triticum")
check <- check & !gbdata$Species %in% toExclude
gbdata$cat[which(check & !gbdata$Species %in% toExclude)] = 'barley'


gbdata  <- subset(gbdata,!is.na(cat))
gbdata$cat2  <- "Wheat + Barley"
gbdata$cat2[which(gbdata$cat=='hazelnut')] = "Hazelnut"
gbdata$cat2 = factor(gbdata$cat2,levels=c('Wheat + Barley','Hazelnut'),ordered=T)

# Focus on 8000 to 3000 BP
xx <- rcarbon::calibrate(gbdata$CRA,gbdata$Error)
gbdata <- gbdata[which.CalDates(xx,BP<7000&BP>3000,p=0.5),]

# Exclude Ireland
gbdata  <- subset(gbdata,Region!='Ireland')

# Count Data
nrow(gbdata) #928
length(unique(gbdata$SiteID)) #314
table(gbdata$cat)

# Store Results
save(gbdata,file=here('data','gbdata.RData'))
