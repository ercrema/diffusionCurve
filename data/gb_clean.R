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
gbdata <- gbdata[which.CalDates(xx,BP<8000&BP>3000,p=0.5),]

# Define Geographic Regions ----
# Exclude Ireland
gbdata  <- subset(gbdata,Region!='Ireland')
# Download UK
uk <- ne_states(country = "united kingdom",returnclass = 'sf')
# Project to OSGB
uk  <- st_transform(uk,27700) #project to OSGB
# Identify sites outside
sites  <- st_as_sf(gbdata,coords=c('Longitude','Latitude'),crs=4326) |> st_transform(x=_,27700)
inside <- lengths(st_intersects(sites, uk)) > 0

# Extract Scotland
scotland  <- subset(uk,region%in%c('Highlands and Islands','North Eastern','South Western',"Eastern")) |> st_union() |> st_cast(to='POLYGON')
# Isolate Mainland Scotland
scotland_mainland  <- scotland[order(st_area(scotland),decreasing=TRUE)[1]]
# Create Buffer 1km Mainland Scotland
scotland_mainland_buffer <- st_buffer(scotland_mainland,500)
# Define England and Wales
england_and_wales  <- subset(uk,!region%in%c('Highlands and Islands','North Eastern','South Western', 'Eastern','Northern Ireland')) |> st_union() |> st_cast(to='POLYGON')
# Creat Buffer Englan and Wales
england_and_wales_buffer <- st_buffer(england_and_wales,500)

# Define binary intersections
scotland_mainland.int <- lengths(st_intersects(sites, scotland_mainland)) > 0
scotland_mainland_buffer.int <- lengths(st_intersects(sites, scotland_mainland_buffer)) > 0
england_and_wales.int  <- lengths(st_intersects(sites,england_and_wales)) > 0
england_and_wales_buffer.int  <- lengths(st_intersects(sites,england_and_wales_buffer)) > 0

# Mainland scotland: in scotland and in 1km buffer but not in england and wales
i.scotland.mainland  <- (scotland_mainland.int | scotland_mainland_buffer.int) & !england_and_wales.int
# England and Wales: in E&W and in 1km buffer but not in scotland main
i.england_and_wales  <- (england_and_wales.int | england_and_wales_buffer.int) & !scotland_mainland.int
# Scottish Isles: recorded as scotland in gbdata, but not in buffer or mainland scotland sf
i.scottish_isles  <- (gbdata$Region=='Scotland') & !scotland_mainland.int & !scotland_mainland_buffer.int

# Check if there are any overlaps (should return FALSE):
any((i.scotland.mainland + i.england_and_wales + i.scottish_isles) > 2)

# Sanity Check plot
uk.win  <- st_union(uk) |> st_cast(to='POLYGON') 
plot(uk.win)
plot(sites,pch=19,col='lightgrey',add=TRUE,cex=0.5)
plot(sites[i.scotland.mainland,],pch=19,col='blue',add=TRUE,cex=0.5)
plot(sites[i.scottish_isles,],pch=19,col='orange',add=TRUE,cex=0.5)
plot(sites[i.england_and_wales,],pch=19,col='darkgreen',add=TRUE,cex=0.5)

# Assign Regions
gbdata$Region2  <- NA
gbdata$Region2[i.scotland.mainland] <- 'scotland_main'
gbdata$Region2[i.scottish_isles] <- 'scottish_isles'
gbdata$Region2[i.england_and_wales] <- 'england_wales'
gbdata$Region2[!inside] <- 'sea'
write.csv(gbdata,file=here('data','gb_check.csv'),row.names=FALSE)

# Final Sanity Check

plot(uk.win)
sites2  <- st_as_sf(gbdata,coords=c('Longitude','Latitude'),crs=4326) |> st_transform(x=_,27700)
plot(sites2,pch=19,col='lightgrey',add=TRUE,cex=0.5)
plot(sites2[gbdata$Region2=='scotmain',],pch=19,col='blue',add=TRUE,cex=0.5)
plot(sites2[gbdata$Region2=='scotisles',],pch=19,col='orange',add=TRUE,cex=0.5)
plot(sites2[gbdata$Region2=='engwal',],pch=19,col='darkgreen',add=TRUE,cex=0.5)

gbdata  <- subset(gbdata,!is.na(Region2))

save(gbdata,file=here('data','gbdata.RData'))





