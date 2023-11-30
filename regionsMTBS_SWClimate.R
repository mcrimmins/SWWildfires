# SW Wildfire analysis - using MTBS perimeters dataset AND EcoRegional Analysis, other regions
# MAC 11/28/23
# adapted from burn period research fodAnalysis.R, FOD_SWClimate.R, regionsFOD_SWClimate.R
# using data from https://www.mtbs.gov/direct-download

# library(tidyverse)
# library(ggplot2)
# 
# ##### CREATE MTBS dataframe with ecoregions, vegtypes -----
# # load/filter mtbs perims, ./BurnPeriodResearch/mtbs_perims.R ----
#  load("~/RProjects/BurnPeriodResearch/data/AZNM_mtbs_perims_1984_2022.RData")
# 
# # remove factors -- from ./BurnPeriodResearch/fireProgression_loop.R
#  perimsCrop$BurnBndAc<-as.numeric(as.character(perimsCrop$BurnBndAc))
#  perimsCrop$Incid_Name<-as.character(perimsCrop$Incid_Name)
#  perimsCrop$Ig_Date<-as.Date(perimsCrop$Ig_Date, "%Y-%m-%d")
#  perimsCrop$year<-as.numeric(format(as.Date(perimsCrop$Ig_Date, format="%Y/%m/%d"),"%Y"))
#  perimsCrop<-subset(perimsCrop, Incid_Type!="Prescribed Fire")
# 
# # map layers
# states <- raster::getData('GADM', country='United States', level=1)
# aznm<-subset(states, NAME_1=="Arizona" | NAME_1=="New Mexico")
# ecoreg<-rgdal::readOGR(dsn="~/RProjects/SOMs/monsoonPrecip/shapes", layer="us_eco_l3")
# ecoreg <- sp::spTransform(ecoreg,raster::crs(states))
# # intersect polys
# ecoreg<-raster::intersect(ecoreg,aznm)
# 
# #  add in ecoreg to points
# #  https://gis.stackexchange.com/questions/270537/joining-attributes-by-location-for-more-than-one-located-feature-using-r
# ecoreg <- sp::spTransform(ecoreg,raster::crs(perimsCrop))
# 
# # # us sf to join
# sf::sf_use_s2(FALSE)
#  perimsCrop <- sf::st_join(sf::st_as_sf(perimsCrop), sf::st_as_sf(ecoreg), left=TRUE)
#  perimsCrop<-perimsCrop[!duplicated(perimsCrop$Event_ID),]
#  perimsCrop <- subset(perimsCrop, !is.na(US_L3CODE))
# 
# # drop some columns
# perimsCrop <- perimsCrop[,c(1:4,6,8:11,23,25,40)]
# 
# # read in landfire data, use Terra library
# # #bps<-raster("./data/landfire/LF2020_BPS_220_CONUS/LC20_BPS_220.tif")
# bps<-terra::rast("./data/landfire/LF2020_BPS_220_CONUS/LC20_BPS_220.tif")
#   bps<- terra::project(bps, "EPSG:4326", method = "near")
# 
# # extract BPS name and groupveg
# perimsCrop$groupVeg<-exactextractr::exact_extract(bps,perimsCrop,'mode')
#   bpsMeta<- read.csv('https://landfire.gov/CSV/LF2020/LF20_BPS_220.csv')
#   perimsCrop<-merge(perimsCrop, bpsMeta[,c(1,6)], by.x="groupVeg",by.y="VALUE")
# # need code to map to groupVeg https://landfire.gov/CSV/LF2020/LF20_BPS_220.csv
# # https://landfire.gov/bps.php
# save(perimsCrop,ecoreg, aznm, file="./data/mtbsSW_ecoreg.RData")
#####

##### mtbs diagnostic plots/summaries----
# yearFire<-perimsCrop %>%
#   group_by(year) %>%
#   summarise(totAc=sum(BurnBndAc),
#             nFires=n())
# 
# ggplot(yearFire, aes(year,nFires))+
#   geom_bar(stat = "identity")
# 
# perimsCropSF <- sf::st_as_sf(perimsCrop)
# ggplot()+
#   geom_sf(data=perimsCropSF, aes(fill=as.factor(year)), color='black', linewidth=0.5, alpha=1)+
#   coord_sf()+
#   ggtitle("MTBS Wildfires, 1984-2022")
#####

##### extract ecoRegion climate time series -----
# 
# library(SPEI)
# library(tidyverse)
# 
# load("~/RProjects/SWWildfires/data/fcSW_ecoreg.RData")
# 
# # combine ecoregions
# ecoregSW <- raster::aggregate(ecoreg, 'US_L3NAME')
#  #sp::spplot(ecoregSW, "US_L3NAME")
# 
# ###
# load PRISM data, adapted from ./WinterSummerPrecip/extractClusterTS.R
# dates 1895-2022 PRISM data
# dates=seq(as.Date("1895-01-01"), as.Date("2022-12-31"), by="month")
# dates<-as.data.frame(dates)
# dates$month<-as.numeric(format(dates$dates, "%m"))
# dates$year<-as.numeric(format(dates$dates, "%Y"))
# length(unique(dates$year))
# 
# # update scratch dir with PRISM data 1895-2022 monthly precip/mean temp
# # use ~/RProjects/PRISMDownload/monthyDownloadPRISM.R
# # process to subset using ~/RProjects/WinterSummerPrecip/processPRISM.R
# climFiles<-as.data.frame(list.files(path = "/scratch/crimmins/PRISM/monthly/processed/SW", pattern="*.grd",full.names = TRUE, recursive = TRUE))
# climFiles$fileName<-(list.files(path = "/scratch/crimmins/PRISM/monthly/processed/SW", pattern="*.grd",full.names = FALSE, recursive = TRUE))
# colnames(climFiles)<-c("path","name")
# climFiles<-climFiles %>% 
#   separate(name, sep="_", into=c("dataset", "var", "being","end"))
# climFiles<-subset(climFiles, end=="2022.grd")
# climFiles$path<-as.character(climFiles$path)
# ###
# extract time series from fire perimeter
#climAnom<-list()
# tictoc::tic()
# ecoClim<-list()
# for(j in 1:nrow(ecoregSW)){
#   
#   print(j)
#   
#   climTS<-list()
#   for(i in 1:nrow(climFiles)){
#     # temp raster var
#     #var<-stack(climFiles$path[i])
#     var<-terra::rast(climFiles$path[i])
#     # extract time series
#     ext<-as.data.frame(t(terra::extract(var, terra::vect(ecoregSW[j,]), ID=FALSE, raw=FALSE, fun=mean)))
#     ext<-as.data.frame(ext[2:nrow(ext),])
#     colnames(ext)<-climFiles$var[i]  
#     climTS[[i]]<-ext
#     #print(climFiles$var[i])
#   }
#   # create fire dataframe
#   climTS<-do.call(cbind, climTS)
#   climTS<-cbind.data.frame(ecoregSW@data$US_L3NAME[j],
#                            dates,climTS)
#   colnames(climTS)[1]<-c("eco1")
#   # calculate SPI
#   climTS$spi3<-spi(climTS$prec,3, na.rm = TRUE)$fitted
#   climTS$spi6<-spi(climTS$prec,6, na.rm = TRUE)$fitted
#   climTS$spi12<-spi(climTS$prec,12, na.rm = TRUE)$fitted
#   # calculate SPEI
#   climTS$spei3<-spei(climTS$prec-climTS$hargreaves,3, na.rm=TRUE)$fitted
#   climTS$spei6<-spei(climTS$prec-climTS$hargreaves,6, na.rm=TRUE)$fitted
#   climTS$spei12<-spei(climTS$prec-climTS$hargreaves,12, na.rm=TRUE)$fitted
#   
#   # 3-month moving sums/avgs
#   clim3mo<-cbind.data.frame(zoo::rollapply(climTS[,c("hargreaves","prec")], FUN = sum, width = 3,
#                                            fill=NA,align="right", by.column = TRUE),
#                             zoo::rollapply(climTS[,c("tdmean","tmax","tmean","tmin","vpdmax")], FUN = mean, width = 3,
#                                            fill=NA,align="right", by.column = TRUE))
#   #colnames(clim3mo)<-paste0(colnames(clim3mo),"3mo")
#   
#   # 3-month anoms, swap in rolling avg/sum
#   anoms<-cbind.data.frame(climTS$dates,climTS$month,clim3mo)
#   colnames(anoms)[1:2]<-c("dates","month")  
#   # calculate anomalies
#   #anoms<-climTS[,6:14]
#   meanMO<- anoms %>% group_by(month) %>%
#     summarise(meanHargreaves=mean(hargreaves, na.rm=TRUE),
#               meanPrec=mean(prec, na.rm=TRUE),
#               meanTDmean=mean(tdmean, na.rm=TRUE),
#               meanTmax=mean(tmax, na.rm=TRUE),
#               meanTmean=mean(tmean, na.rm=TRUE),
#               meanTmin=mean(tmin, na.rm=TRUE),
#               meanVPDmax=mean(vpdmax, na.rm=TRUE))
#   anoms<-merge(climTS,meanMO, by=c("month"))
#   # calculate anoms
#   anoms$anomHarg<-anoms$hargreaves-anoms$meanHargreaves
#   anoms$anomPrec<-anoms$prec-anoms$meanPrec
#   anoms$anomTDmean<-anoms$tdmean-anoms$meanTDmean
#   anoms$anomTmax<-anoms$tmax-anoms$meanTmax
#   anoms$anomTmean<-anoms$tmean-anoms$meanTmean
#   anoms$anomTmin<-anoms$tmin-anoms$meanTmin
#   anoms$anomVPDmax<-anoms$vpdmax-anoms$meanVPDmax
#   anoms <- anoms[order(anoms$dates),]
#   
#   # fixed 3-mo seasons
#   anoms$seas<-cut(anoms$month,c(0,3,6,9,12))
#   levels(anoms$seas) = c("JFM","AMJ","JAS","OND")
#   anoms<-subset(anoms, month %in% c(3,6,9,12))
#   
#   # write to list
#   ecoClim[[j]]<-anoms  
#   
#   # save full anoms
#   #climAnom[[j]]<-anoms
#   
# }
# tictoc::toc()
# 
# save(ecoClim, file="./data/AZNM_ecoregion_climate_3moRoll.RData")
######


#####

##### analyze fires and ecoregion climate ----
library(tidyverse)
library(ggplot2)

load("~/RProjects/SWWildfires/data/mtbsSW_ecoreg.RData")
# see regionsMTBS_SWClimate.R
load("./data/AZNM_ecoregion_climate_3moRoll.RData")

# filter out border state fires based on event ID
perimsCrop$event_st<-substr(perimsCrop$Event_ID,1,2)
perimsCrop<-subset(perimsCrop, event_st %in% c("AZ","NM"))

# create DF
mtbsDF<-sf::st_drop_geometry(perimsCrop)
##### summary of fc table ----
# map of ecoreg
ecoregSF <- sf::st_as_sf(ecoreg)
ggplot()+
  geom_sf(data=ecoregSF, aes(), color='grey75', linewidth=0.5, alpha=1)+
  geom_sf(data=perimsCrop, aes(fill=as.factor(year)), color='black', linewidth=0.25, alpha=1)+
  coord_sf()+
  ggtitle("MTBS Fire Locations, 1984-2021")

#plotly::ggplotly(p)

# create summary table
temp<-mtbsDF %>%
  group_by(year) %>%
  summarise(nFires=n(),
            totAc=sum(BurnBndAc))
temp$ratio<-temp$totAc/temp$nFires

ggplot(temp, aes(year,nFires))+
  geom_bar(position="stack", stat="identity")

temp<-mtbsDF %>%
  group_by(year,GROUPVEG) %>%
  summarise(nFires=n(),
            totAc=sum(BurnBndAc))

ggplot(temp, aes(year,totAc, fill=GROUPVEG))+
  geom_bar(position="stack", stat="identity")
#####



