#Author: Elly C Knight
#Date: Oct 26, 2019
#Probability of occurrence models for targeted and general nightjar survey data 

options(scipen = 999)

setwd("/Volumes/ECK001/TrendsAnalysis")

library(tidyverse)
library(lme4) #for mixed effects models
library(usdm) #for VIF
library(MuMIn) #for AIC
library(raster) #for raster data
library(rgdal) #for shapefiles

bbs <- read.csv("BBSWrangledAll.csv")
cns <- read.csv("WRWrangledAll.csv")

#1. WRANGLE####
#1a. Subset to western provinces & territories, after year 2000, center & standardize----
#BBS
bbs.coni <- bbs %>% 
  dplyr::filter(coni.range==1,
                year >=2000) %>% 
  dplyr::filter(region %in% c(4, 11, 43, 79, 93))

bbs.coni.st <- bbs.coni %>% 
  dplyr::select(sun, doy, lat, long) %>% 
  scale() %>% 
  data.frame() %>% 
  rename(doy.st=doy, sun.st=sun, lat.st=lat, long.st=long) %>% 
  cbind(bbs.coni)

#CNS
cns.coni <- cns %>% 
  dplyr::filter(coni.range==1,
                year >=2000) %>% 
  dplyr::filter(region %in% c(4, 11, 43, 79, 93))

cns.coni.st <- cns.coni %>% 
  dplyr::select(sun, doy, lat, long) %>% 
  scale() %>% 
  data.frame() %>% 
  rename(doy.st=doy, sun.st=sun, lat.st=lat, long.st=long) %>% 
  cbind(cns.coni)

#2. VISUALIZE####
#BBS
ggplot(bbs.coni.st, aes(x=sun.st, y=coni.pres)) +
  geom_smooth() #3rd order

ggplot(bbs.coni.st, aes(x=doy.st, y=coni.pres)) +
  geom_smooth() #2nd order

#CNS
ggplot(cns.coni.st, aes(x=sun.st, y=coni.pres)) +
  geom_smooth() #2nd order

ggplot(cns.coni.st, aes(x=doy.st, y=coni.pres)) +
  geom_smooth() #2nd order

#3. 'DETECTION' OCCURRENCE MODELS####
#3a. Find best model----
#BBS
bbs.coni.det.global <- glmer(coni.pres ~ poly(sun.st, degree=3)*poly(doy.st, degree=2) +
                         (1|route),
                       family="binomial", data=bbs.coni.st, na.action="na.fail", 
                       control = glmerControl(optimizer ="Nelder_Mead"))
summary(bbs.coni.det.global)
bbs.coni.dredge <- dredge(bbs.coni.det.global, trace=2)
bbs.coni.dredge
importance(bbs.coni.dredge)

#CNS
cns.coni.det.global <- glmer(coni.pres ~ poly(sun.st, degree=2)*poly(doy.st, degree=3) +
                               (1|route), 
                             family="binomial", data=cns.coni.st, na.action="na.fail", 
                             control = glmerControl(optimizer ="Nelder_Mead"))
summary(cns.coni.det.global)
cns.coni.dredge <- dredge(cns.coni.det.global, trace=2)
cns.coni.dredge
importance(cns.coni.dredge)

#3b. Fit best model----
#BBS
bbs.coni.det <- glmer(coni.pres ~ poly(sun.st, degree=3) + poly(doy.st, degree=2) +
                        (1|route), 
                      family="binomial", data=bbs.coni.st, na.action="na.fail", 
                      control = glmerControl(optimizer ="Nelder_Mead"))
summary(bbs.coni.det)

#CNS
cns.coni.det <- glmer(coni.pres ~ poly(sun.st, degree=2) +
                       (1|route), 
                     family="binomial", data=cns.coni.st, na.action="na.fail", 
                     control = glmerControl(optimizer ="Nelder_Mead"))
summary(cns.coni.det)

#3c. Calculate offsets----
#BBS
off <- predict(bbs.coni.det, bbs.coni.st, type="response", re.form=NA)
bbs.coni.off <- cbind(bbs.coni.st, off)

#CNS
off <- predict(cns.coni.det, cns.coni.st, type="response", re.form=NA)
cns.coni.off <- cbind(cns.coni.st, off)

#4. HABITAT COVARIATES - BCR 11####
#4a. Subset data to BCR 11----
#BBS
bbs.11 <- bbs.coni.off %>% 
  dplyr::filter(BCR==11) %>% 
  rename(pres=coni.pres)

#CNS
cns.11 <- cns.coni.off %>% 
  dplyr::filter(BCR==11) %>% 
  rename(pres=coni.pres)

#4b. Read in spatial data----
#Spatial rasters ("Layers_WGS84V2_BCR11.tif") available at https://drive.google.com/file/d/1iHRGk4Ft2YxMrnGp0c-Z1jp-UW3iHIpu/view?usp=sharing

ag.rd <- raster("Layers_WGS84V2_BCR11.tif", band=1)
names(ag.rd) <- "Agriculture"
water.rd <- raster("Layers_WGS84V2_BCR11.tif", band=2)
names(water.rd) <- "Water"
dev.rd <- raster("Layers_WGS84V2_BCR11.tif", band=3)
names(dev.rd) <- "Development"
grass.rd <- raster("Layers_WGS84V2_BCR11.tif", band=4)
names(grass.rd) <- "Grassland"
shrub.rd <- raster("Layers_WGS84V2_BCR11.tif", band=5)
names(shrub.rd) <- "Shrubland"
min.rd <- raster("Layers_WGS84V2_BCR11.tif", band=6)
names(min.rd) <- "Wetland"
taiga.rd <- raster("Layers_WGS84V2_BCR11.tif", band=7)
names(taiga.rd) <- "Taiga"
wet.rd <- raster("Layers_WGS84V2_BCR11.tif", band=8)
names(wet.rd) <- "Wet"
decid.rd <- raster("Layers_WGS84V2_BCR11.tif", band=9)
names(decid.rd) <- "Deciduous"
conif.rd <- raster("Layers_WGS84V2_BCR11.tif", band=10)
names(conif.rd) <- "Coniferous"
mix.rd <- raster("Layers_WGS84V2_BCR11.tif", band=11)
names(mix.rd) <- "Mixedwood"
tree.rd <- raster("Layers_WGS84V2_BCR11.tif", band=12)
names(tree.rd) <- "Tree"
temp.rd <- raster("Layers_WGS84V2_BCR11.tif", band=13)
names(temp.rd) <- "Temperature"
pre.rd <- raster("Layers_WGS84V2_BCR11.tif", band=14)
names(pre.rd) <- "Precipitation"
dem.rd <- raster("Layers_WGS84V2_BCR11.tif", band=15)
names(dem.rd) <- "Elevation"
height.rd <- raster("Layers_WGS84V2_BCR11.tif", band=16)
names(height.rd) <- "Height"

layers.cut <- stack(ag.rd, 
                water.rd,
                dev.rd,
                grass.rd,
                shrub.rd,
                min.rd,
                taiga.rd,
                wet.rd,
                decid.rd,
                conif.rd,
                mix.rd,
                tree.rd,
                temp.rd,
                pre.rd,
                dem.rd,
                height.rd)

#4c. Extract covariates for survey points & filter out those outside the spatial data extent----
#BBS
covs.all <- data.frame(raster::extract(layers.cut, bbs.11[,c("long", "lat")]))
bbs.covs <- cbind(bbs.11, covs.all)
bbs.use1 <- bbs.covs %>% 
  dplyr::filter(!is.na(Water)) %>% 
  dplyr::mutate(dat="bbs",
                ind=1:n(),
                ID=paste0(dat,"-",ind))

#CNS
covs.all <- data.frame(raster::extract(layers.cut, cns.11[,c("long", "lat")]))
cns.covs <- cbind(cns.11, covs.all)
cns.use1 <- cns.covs %>%  
  dplyr::filter(!is.na(Water)) %>% 
  mutate(dat="cns",
         ind=1:n(),
         ID=paste0(dat,"-",ind))

#Combined
cnsbbs.use1 <- cns.use1 %>% 
  rbind(bbs.use1)

#4d. Plot to check extent----
#BBS
bbs.present <- subset(bbs.use1, pres==1)
bbs.absent <- subset(bbs.use1, pres==0)
plot(layers.cut, 6)
points(bbs.absent$long, bbs.absent$lat, pch = 16, col="black")
points(bbs.present$long, bbs.present$lat, pch=16, col="red")

#CNS
cns.present <- subset(cns.use1, pres==1)
cns.absent <- subset(cns.use1, pres==0)
plot(layers.cut, 6)
points(cns.absent$long, cns.absent$lat, pch = 16, col="black")
points(cns.present$long, cns.present$lat, pch=16, col="red")

#Combined
cnsbbs.present <- subset(cnsbbs.use1, pres==1)
cnsbbs.absent <- subset(cnsbbs.use1, pres==0)
plot(layers.cut, 6)
points(cnsbbs.absent$long, cnsbbs.absent$lat, pch = 16, col="black")
points(cnsbbs.present$long, cnsbbs.present$lat, pch=16, col="red")

#4e. Center & standardize survey data----
#BBS
bbs.use2 <- bbs.use1 %>% 
  dplyr::select(doy,
                sun,
                Agriculture,
                Water,
                Grassland,
                Shrubland,
                Wetland,
                Taiga,
                Wet,
                Deciduous,
                Coniferous,
                Mixedwood,
                Tree,
                Temperature,
                Precipitation,
                Elevation,
                Height) %>% 
  mutate(doy=as.numeric(doy))

bbs.covs.st1 <- scale(bbs.use2)
bbs.covs.st <- data.frame(bbs.covs.st1)
names(bbs.covs.st) <- paste0(names(bbs.covs.st), ".st")
bbs.use3 <- cbind(ID=bbs.use1$ID, pres=bbs.use1$pres, route=bbs.use1$route, lat=bbs.use1$lat, long=bbs.use1$long, off=bbs.use1$off, bbs.use2, bbs.covs.st)
write.csv(bbs.use3, "BBSStandardizedData_BCR11.csv", row.names = FALSE)

#CNS
cns.use2 <- cns.use1 %>% 
  dplyr::select(doy,
                sun,
                Agriculture,
                Water,
                Grassland,
                Shrubland,
                Wetland,
                Taiga,
                Wet,
                Deciduous,
                Coniferous,
                Mixedwood,
                Tree,
                Temperature,
                Precipitation,
                Elevation,
                Height) %>% 
  mutate(doy=as.numeric(doy))

cns.covs.st1 <- scale(cns.use2)
cns.covs.st <- data.frame(cns.covs.st1)
names(cns.covs.st) <- paste0(names(cns.covs.st), ".st")
cns.use3 <- cbind(ID=cns.use1$ID, pres=cns.use1$pres, route=cns.use1$route, lat=cns.use1$lat, long=cns.use1$long, off=cns.use1$off, cns.use2, cns.covs.st)
write.csv(cns.use3, "WRStandardizedData_BCR11.csv", row.names = FALSE)

#Combined
cnsbbs.use2 <- cnsbbs.use1 %>% 
  dplyr::select(doy,
                sun,
                Agriculture,
                Water,
                Grassland,
                Shrubland,
                Wetland,
                Taiga,
                Wet,
                Deciduous,
                Coniferous,
                Mixedwood,
                Tree,
                Temperature,
                Precipitation,
                Elevation,
                Height) %>% 
  mutate(doy=as.numeric(doy))

cnsbbs.covs.st1 <- scale(cnsbbs.use2)
cnsbbs.covs.st <- data.frame(cnsbbs.covs.st1)
names(cnsbbs.covs.st) <- paste0(names(cnsbbs.covs.st), ".st")
cnsbbs.use3 <- cbind(ID=cnsbbs.use1$ID, pres=cnsbbs.use1$pres, route=cnsbbs.use1$route, lat=cnsbbs.use1$lat, long=cnsbbs.use1$long, off=cnsbbs.use1$off, cnsbbs.use2, cnsbbs.covs.st)
write.csv(cnsbbs.use3, "WRBBSStandardizedData_BCR11.csv", row.names = FALSE)

#4f. Center & standardize spatial data (for prediction)----
#BBS
#Extract means and sds
center <- attr(bbs.covs.st1, "scaled:center")
scale <- attr(bbs.covs.st1, "scaled:scale")

#Standardize rasters
ag.rad.st <- (layers.cut[[1]]-center[["Agriculture"]])/scale[["Agriculture"]]
names(ag.rad.st) <- "Agriculture.st"
water.rad.st <- (layers.cut[[2]]-center[["Water"]])/scale[["Water"]]
names(water.rad.st) <- "Water.st"
grass.rad.st <- (layers.cut[[4]]-center[["Grassland"]])/scale[["Grassland"]]
names(grass.rad.st) <- "Grassland.st"
shrub.rad.st <- (layers.cut[[5]]-center[["Shrubland"]])/scale[["Shrubland"]]
names(shrub.rad.st) <- "Shrubland.st"
min.rad.st <- (layers.cut[[6]]-center[["Wetland"]])/scale[["Wetland"]]
names(min.rad.st) <- "Wetland.st"
wet.rad.st <- (layers.cut[[8]]-center[["Wet"]])/scale[["Wet"]]
names(wet.rad.st) <- "Wet.st"
decid.rad.st <- (layers.cut[[9]]-center[["Deciduous"]])/scale[["Deciduous"]]
names(decid.rad.st) <- "Deciduous.st"
conif.rad.st <- (layers.cut[[10]]-center[["Coniferous"]])/scale[["Coniferous"]]
names(conif.rad.st) <- "Coniferous.st"
mix.rad.st <- (layers.cut[[11]]-center[["Mixedwood"]])/scale[["Mixedwood"]]
names(mix.rad.st) <- "Mixedwood.st"
tree.rad.st <- (layers.cut[[12]]-center[["Tree"]])/scale[["Tree"]]
names(tree.rad.st) <- "Tree.st"
temp.cut.st <- (layers.cut[[13]]-center[["Temperature"]])/scale[["Temperature"]]
names(temp.cut.st) <- "Temperature.st"
precip.cut.st <- (layers.cut[[14]]-center[["Precipitation"]])/scale[["Precipitation"]]
names(precip.cut.st) <- "Precipitation.st"
elev.cut.st <- (layers.cut[[15]]-center[["Elevation"]])/scale[["Elevation"]]
names(elev.cut.st) <- "Elevation.st"
height.cut.st <- (layers.cut[[16]]-center[["Height"]])/scale[["Height"]]
names(height.cut.st) <- "Height.st"

#Create new stack
layers.bbs <- stack(water.rad.st, grass.rad.st, shrub.rad.st, min.rad.st, wet.rad.st, decid.rad.st, conif.rad.st, mix.rad.st, tree.rad.st, temp.cut.st, precip.cut.st, elev.cut.st, height.cut.st)
writeRaster(layers.bbs, "BBSStandardizedLayersAll_BCR11.tif", format="GTiff", overwrite=TRUE)

#WR
#Retrieve means & sds
center <- attr(cns.covs.st1, "scaled:center")
scale <- attr(cns.covs.st1, "scaled:scale")

#Standardize rasters
ag.rad.st <- (layers.cut[[1]]-center[["Agriculture"]])/scale[["Agriculture"]]
names(ag.rad.st) <- "Agriculture.st"
water.rad.st <- (layers.cut[[2]]-center[["Water"]])/scale[["Water"]]
names(water.rad.st) <- "Water.st"
grass.rad.st <- (layers.cut[[4]]-center[["Grassland"]])/scale[["Grassland"]]
names(grass.rad.st) <- "Grassland.st"
shrub.rad.st <- (layers.cut[[5]]-center[["Shrubland"]])/scale[["Shrubland"]]
names(shrub.rad.st) <- "Shrubland.st"
min.rad.st <- (layers.cut[[6]]-center[["Wetland"]])/scale[["Wetland"]]
names(min.rad.st) <- "Wetland.st"
wet.rad.st <- (layers.cut[[8]]-center[["Wet"]])/scale[["Wet"]]
names(wet.rad.st) <- "Wet.st"
decid.rad.st <- (layers.cut[[9]]-center[["Deciduous"]])/scale[["Deciduous"]]
names(decid.rad.st) <- "Deciduous.st"
conif.rad.st <- (layers.cut[[10]]-center[["Coniferous"]])/scale[["Coniferous"]]
names(conif.rad.st) <- "Coniferous.st"
mix.rad.st <- (layers.cut[[11]]-center[["Mixedwood"]])/scale[["Mixedwood"]]
names(mix.rad.st) <- "Mixedwood.st"
tree.rad.st <- (layers.cut[[12]]-center[["Tree"]])/scale[["Tree"]]
names(tree.rad.st) <- "Tree.st"
temp.cut.st <- (layers.cut[[13]]-center[["Temperature"]])/scale[["Temperature"]]
names(temp.cut.st) <- "Temperature.st"
precip.cut.st <- (layers.cut[[14]]-center[["Precipitation"]])/scale[["Precipitation"]]
names(precip.cut.st) <- "Precipitation.st"
elev.cut.st <- (layers.cut[[15]]-center[["Elevation"]])/scale[["Elevation"]]
names(elev.cut.st) <- "Elevation.st"
height.cut.st <- (layers.cut[[16]]-center[["Height"]])/scale[["Height"]]
names(height.cut.st) <- "Height.st"

#Create new stack
layers.cns <- stack(water.rad.st, grass.rad.st, shrub.rad.st, min.rad.st, wet.rad.st, decid.rad.st, conif.rad.st, mix.rad.st, tree.rad.st, temp.cut.st, precip.cut.st, elev.cut.st, height.cut.st)
setwd("/Volumes/ECK001/TrendsAnalysis")
writeRaster(layers.cns, "CNSStandardizedLayersAll_BCR11.tif", format="GTiff", overwrite=TRUE)

#Combined
#Retrieve means and sds
center <- attr(cnsbbs.covs.st1, "scaled:center")
scale <- attr(cnsbbs.covs.st1, "scaled:scale")

#Standardize rasters
ag.rad.st <- (layers.cut[[1]]-center[["Agriculture"]])/scale[["Agriculture"]]
names(ag.rad.st) <- "Agriculture.st"
water.rad.st <- (layers.cut[[2]]-center[["Water"]])/scale[["Water"]]
names(water.rad.st) <- "Water.st"
grass.rad.st <- (layers.cut[[4]]-center[["Grassland"]])/scale[["Grassland"]]
names(grass.rad.st) <- "Grassland.st"
shrub.rad.st <- (layers.cut[[5]]-center[["Shrubland"]])/scale[["Shrubland"]]
names(shrub.rad.st) <- "Shrubland.st"
min.rad.st <- (layers.cut[[6]]-center[["Wetland"]])/scale[["Wetland"]]
names(min.rad.st) <- "Wetland.st"
wet.rad.st <- (layers.cut[[8]]-center[["Wet"]])/scale[["Wet"]]
names(wet.rad.st) <- "Wet.st"
decid.rad.st <- (layers.cut[[9]]-center[["Deciduous"]])/scale[["Deciduous"]]
names(decid.rad.st) <- "Deciduous.st"
conif.rad.st <- (layers.cut[[10]]-center[["Coniferous"]])/scale[["Coniferous"]]
names(conif.rad.st) <- "Coniferous.st"
mix.rad.st <- (layers.cut[[11]]-center[["Mixedwood"]])/scale[["Mixedwood"]]
names(mix.rad.st) <- "Mixedwood.st"
tree.rad.st <- (layers.cut[[12]]-center[["Tree"]])/scale[["Tree"]]
names(tree.rad.st) <- "Tree.st"
temp.cut.st <- (layers.cut[[13]]-center[["Temperature"]])/scale[["Temperature"]]
names(temp.cut.st) <- "Temperature.st"
precip.cut.st <- (layers.cut[[14]]-center[["Precipitation"]])/scale[["Precipitation"]]
names(precip.cut.st) <- "Precipitation.st"
elev.cut.st <- (layers.cut[[15]]-center[["Elevation"]])/scale[["Elevation"]]
names(elev.cut.st) <- "Elevation.st"
height.cut.st <- (layers.cut[[16]]-center[["Height"]])/scale[["Height"]]
names(height.cut.st) <- "Height.st"

#Create new stack
layers.cnsbbs <- stack(water.rad.st, grass.rad.st, shrub.rad.st, min.rad.st, wet.rad.st, decid.rad.st, conif.rad.st, mix.rad.st, tree.rad.st, temp.cut.st, precip.cut.st, elev.cut.st, height.cut.st)
setwd("/Volumes/ECK001/TrendsAnalysis")
writeRaster(layers.cnsbbs, "WRBBSStandardizedLayersAll_BCR11.tif", format="GTiff", overwrite=TRUE)

#4g. Check for collinearity----
#BBS
bbs.covs.colin <- bbs.use3 %>% 
  dplyr::select(Water,
                Agriculture,
                Grassland,
                Shrubland,
                Wetland,
                Tree,
                Temperature,
                Precipitation,
                Elevation,
                Height)

vifstep(bbs.covs.colin, th=10) #Take out agriculture

#CNS
cns.covs.colin <- cns.use3 %>% 
  dplyr::select(Water,
                Agriculture,
                Grassland,
                Shrubland,
                Wetland,
                Tree,
                Temperature,
                Precipitation,
                Elevation,
                Height)

vifstep(cns.covs.colin, th=10) #Take out agriculture

#Combined
cnsbbs.covs.colin <- cnsbbs.use3 %>% 
  dplyr::select(Water,
                Agriculture,
                Grassland,
                Shrubland,
                Wetland,
                Tree,
                Temperature,
                Precipitation,
                Elevation,
                Height)

vifstep(cnsbbs.covs.colin, th=10) #Take out agriculture

#5. BOOTSTRAPPED 'HABITAT' OCCURRENCE MODELS & PREDICTIONS####

#5a. Set up for bootstap----
boot <- 100
coeff <- data.frame()
eval <- data.frame()
eval.df <- data.frame()

for (i in 1:boot){
  
  time <- Sys.time()
  
  #6b. Split data into training & testing
  bbs.ind <- sample(1:nrow(bbs.use3), 0.7*nrow(bbs.use3))
  bbs.train <- bbs.use3[bbs.ind,]
  bbs.test <- bbs.use3[-bbs.ind,] %>% 
    dplyr::filter(pres==1)
  
  cns.ind <- sample(1:nrow(cns.use3), 0.7*nrow(cns.use3))
  cns.train <- cns.use3[cns.ind,]
  cns.test <- cns.use3[-cns.ind,]
  
  cnsbbs.train <- cns.train %>% 
    rbind(bbs.train) %>% 
    dplyr::select(ID) %>% 
    left_join(cnsbbs.use3, by="ID")
  
  test <- cns.test %>% 
    rbind(bbs.test)
  
  #5b. Model----
  #BBS
  bbs.glmer <- glmer(pres ~ Water.st +
                       Grassland.st +
                       Shrubland.st +
                       Wetland.st +
                       Tree.st +
                       Temperature.st +
                       Precipitation.st +
                       Elevation.st +
                       Height.st + (1|route),
                     offset = off,
                     family="binomial",
                     data=bbs.train,
                     na.action = "na.fail")
  
  bbs.summary <- summary(bbs.glmer)
  bbs.coeff <- data.frame(bbs.summary$coefficients)
  bbs.coeff$Var <- row.names(bbs.coeff)
  bbs.coeff$boot <- i
  bbs.coeff$mod <- "BBS"
  
  print(paste0("Completed BBS model ", i))
  
  #CNS
  cns.glmer <- glmer(pres ~ Water.st +
                       Grassland.st +
                       Shrubland.st +
                       Wetland.st +
                       Tree.st +
                       Temperature.st +
                       Precipitation.st +
                       Elevation.st +
                       Height.st + (1|route),
                    offset = off,
                    family="binomial",
                    data=cns.train,
                    na.action = "na.fail")
  
  cns.summary <- summary(cns.glmer)
  cns.coeff <- data.frame(cns.summary$coefficients)
  cns.coeff$Var <- row.names(cns.coeff)
  cns.coeff$boot <- i
  cns.coeff$mod <- "cns"
  
  print(paste0("Completed CNS model ", i))
  
  #Combined
  cnsbbs.glmer <- glmer(pres ~ Water.st +
                          Grassland.st +
                          Shrubland.st +
                          Wetland.st +
                          Tree.st +
                          Temperature.st +
                          Precipitation.st +
                          Elevation.st +
                          Height.st + (1|route),
                       family="binomial",
                       data=cnsbbs.train,
                       na.action = "na.fail")
  
  cnsbbs.summary <- summary(cnsbbs.glmer)
  cnsbbs.coeff <- data.frame(cnsbbs.summary$coefficients)
  cnsbbs.coeff$Var <- row.names(cnsbbs.coeff)
  cnsbbs.coeff$boot <- i
  cnsbbs.coeff$mod <- "CNSBBS"
  
  print(paste0("Completed CNSBBS model ", i))
  
  coeff <- rbind(coeff, bbs.coeff, cns.coeff, cnsbbs.coeff)
  
  #5c. Predict----
  bbs.pred <- raster::predict(layers.bbs, bbs.glmer, type="response", re.form=~0)
  writeRaster(bbs.pred, paste0("bbs.pred.11.",i,".tiff"), format="GTiff", overwrite=TRUE)
  cns.pred <- raster::predict(layers.cns, cns.glmer, type="response", re.form=~0)
  writeRaster(cns.pred, paste0("cns.pred.11.",i,".tiff"), format="GTiff", overwrite=TRUE)
  cnsbbs.pred <- raster::predict(layers.cnsbbs, cnsbbs.glmer, type="response", re.form=~0)
  writeRaster(cnsbbs.pred, paste0("cnsbbs.pred.11.",i,".tiff"), format="GTiff", overwrite=TRUE)
  
  #5d. Evaluate----
  bbs.suit <- data.frame(suitability = raster::extract(bbs.pred, test[,c("long", "lat")]))
  bbs.val <- cbind(test, bbs.suit)
  bbs.e <- dismo::evaluate(subset(bbs.val, pres==1)$suitability, subset(bbs.val, pres==0)$suitability)
  bbs.eval <- data.frame(boot=i, mode="BBS")
  bbs.eval$np <- bbs.e@np
  bbs.eval$na <- bbs.e@na
  bbs.eval$auc <- bbs.e@auc
  bbs.eval$cor <- bbs.e@cor
  bbs.eval$cor <- bbs.e@pcor
  bbs.eval$odp <- mean(bbs.e@ODP)
  bbs.eval.df <- data.frame(t = bbs.e@t) %>% 
    mutate(prev = bbs.e@prevalence,
           ccr = bbs.e@CCR,
           tpr = bbs.e@TPR,
           tnr = bbs.e@TNR,
           fpr = bbs.e@FPR,
           fnr = bbs.e@FNR,
           ppp = bbs.e@PPP,
           npp = bbs.e@NPP,
           mcr = bbs.e@MCR,
           or = bbs.e@OR,
           kappa = bbs.e@kappa,
           boot = i,
           mode="BBS")
  
  cns.suit <- data.frame(suitability = raster::extract(cns.pred, test[,c("long", "lat")]))
  cns.val <- cbind(test, cns.suit)
  cns.e <- dismo::evaluate(subset(cns.val, pres==1)$suitability, subset(cns.val, pres==0)$suitability)
  cns.eval <- data.frame(boot=i, mode="CNS")
  cns.eval$np <- cns.e@np
  cns.eval$na <- cns.e@na
  cns.eval$auc <- cns.e@auc
  cns.eval$cor <- cns.e@cor
  cns.eval$cor <- cns.e@pcor
  cns.eval$odp <- mean(cns.e@ODP)
  cns.eval.df <- data.frame(t = cns.e@t) %>% 
    mutate(prev = cns.e@prevalence,
           ccr = cns.e@CCR,
           tpr = cns.e@TPR,
           tnr = cns.e@TNR,
           fpr = cns.e@FPR,
           fnr = cns.e@FNR,
           ppp = cns.e@PPP,
           npp = cns.e@NPP,
           mcr = cns.e@MCR,
           or = cns.e@OR,
           kappa = cns.e@kappa,
           boot = i,
           mode="CNS")
  
  cnsbbs.suit <- data.frame(suitability = raster::extract(cnsbbs.pred, test[,c("long", "lat")]))
  cnsbbs.val <- cbind(test, cnsbbs.suit)
  cnsbbs.e <- dismo::evaluate(subset(cnsbbs.val, pres==1)$suitability, subset(cnsbbs.val, pres==0)$suitability)
  cnsbbs.eval <- data.frame(boot=i, mode="CNSBBS")
  cnsbbs.eval$np <- cnsbbs.e@np
  cnsbbs.eval$na <- cnsbbs.e@na
  cnsbbs.eval$auc <- cnsbbs.e@auc
  cnsbbs.eval$cor <- cnsbbs.e@cor
  cnsbbs.eval$cor <- cnsbbs.e@pcor
  cnsbbs.eval$odp <- mean(cnsbbs.e@ODP)
  cnsbbs.eval.df <- data.frame(t = cnsbbs.e@t) %>% 
    mutate(prev = cnsbbs.e@prevalence,
           ccr = cnsbbs.e@CCR,
           tpr = cnsbbs.e@TPR,
           tnr = cnsbbs.e@TNR,
           fpr = cnsbbs.e@FPR,
           fnr = cnsbbs.e@FNR,
           ppp = cnsbbs.e@PPP,
           npp = cnsbbs.e@NPP,
           mcr = cnsbbs.e@MCR,
           or = cnsbbs.e@OR,
           kappa = cnsbbs.e@kappa,
           boot = i,
           mode="CNSBBS")
  
  eval <- rbind(eval, cns.eval, bbs.eval, cnsbbs.eval)
  eval.df <- rbind(eval.df, cns.eval.df, bbs.eval.df, cnsbbs.eval.df)
  
  elapsed <- Sys.time() - time
  
  print(paste0("Completed bootstrap ", i, " in ", elapsed, " minutes"))
  
  write.csv(eval, "Evaluation11.csv", row.names=FALSE)
  write.csv(coeff, "Coefficients11.csv", row.names=FALSE)
  write.csv(eval.df, "EvaluationDataFrame11.csv", row.names=FALSE)
  
}
