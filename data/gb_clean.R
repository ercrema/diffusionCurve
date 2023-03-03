# Load Library ----
library(here)
library(rcarbon)
library(sf)
library(rnaturalearth)

# Read Data ----
temp <- tempfile()
download.file("https://discovery.ucl.ac.uk/id/eprint/10025178/4/Bevan_gbie14Csub.zip",temp)
gbdata <- read.csv(unz(temp, "gbie14Csub/dates/archdates.csv"))

# Subset to Great Britain without scottish Isles ----
uk <- ne_states(country = "united kingdom",returnclass = 'sf')
scotland  <- subset(uk,region%in%c('Highlands and Islands','North Eastern','South Western',"Eastern")) |> st_union() |> st_cast(to='POLYGON')
scotland  <- scotland[order(st_area(scotland),decreasing=TRUE)[1]]
england_and_wales  <- subset(uk,!region%in%c('Highlands and Islands','North Eastern','South Western', 'Eastern','Northern Ireland')) |> st_union() |> st_cast(to='POLYGON')
win = st_union(scotland,england_and_wales) |> st_union()

# Subset only sites within window ----
sites  <- st_as_sf(gbdata,coords=c('Longitude','Latitude'),crs=4326)
i   <- st_contains(win,sites,sparse = F,prepared = F)
gbdata  <- gbdata[i,]


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

save(gbdata,file=here('data','gbdata.RData'))





