##
## Second exercises on Ground truthing
# June 11th 2024

# Packages
install.packages("devtools")
install_github("akamoske/canopyLazR")
install.packages("rLiDAR")
install.packages("lidR")
install.packages("terra")
install.packages("rayshader")

# Load the library
library("usethis")
library(devtools)
library(fields)
library(spam)
library(plyr)
library(rlas)
library(raster)
library(sp)
library(canopyLazR)
library(rLiDAR) #nict verfügbar
library(forestr)
library(lidR)
library(rayshader)
library(RCSF)

setwd("C:\\Users\\abdul\\Desktop\\Data Claudia")

rm(list=ls())
## Question:

## How estimates of height based on canopy height model obtained by the TLS is comparable with Sentinel-2?

## Data

## 1) Sentinel - 2 height data
## You already downloaded from the moddle 
# (link:https://drive.google.com/drive/folders/1Jn1_xpkWe2nE3kYgdDhV0wgtwlFw-R6D?usp=sharing)

## 2) TLS data
# Please download from here:
#(link: https://drive.google.com/drive/folders/1jN8P54ZMgyU54ybEqrmTtQUiGxS_gwCu?usp=sharing)

### Read TLS data 

tls <- readLAS("Gehöltzbiomass_1.las") 
summary(tls)
las_check(tls)
epsg(tls)


### Ground classification
# algorithm pmf means Ground segmentation algorithm described by Zhang et al. 2003 who implemented ground classification
# based on 3D point cloud and not on raster-based. See help function with pmf (adapted to airborne data). 
# the computation occurs on a sequence of windows sizes and thresholds.
# There are other ground segmentation algorithms as gnd_csf and gnd_mcc

tls_gc <- classify_ground(tls, algorithm = csf()) ## cloth filter

plot(tls_gc, color = "Classification", size = 3, bg = "white")

###############################################################################
# Digital terrain model

# After processing the classification of what are the ground points and above ground points,
# one can describe an "image" of the ground points as mesh. 
# There many different algorithms to calculate DTM, such as triangular irregular network, Invert distance weighting,
# Kriging, etc. It is basically a interpolation of points to describe the ground surface.

###############################################################################
# Calculate a DTM using a Triangular irregular network (TIN) approach
###############################################################################

# It uses the nearest neighbour to complete the missing pixel out of the convex hull of the ground points

dtm_tin_tls <- rasterize_terrain(tls_gc, res= 0.05, algorithm= tin())

plot_dtm3d(dtm_tin_tls, bg= "white")

###############################################################################

library(rayshader)
dtm_tls <- raster::raster(dtm_tin_tls)
dtm_tls <- as(dtm_tls, "Raster")
plot(dtm_tls, main="Rasterized Digital Terrain Model")


###############################################################################
## Height normalization
################################################################################

# Remove the influence of the terrain on above ground measurements. This allows
# a directly comparison of vegetation heights and different analysis across area,
# plots, etc.
# We use the terrain surface to normalize with 

nlas <- tls - dtm_tin_tls # normalization
plot(nlas, size =4, bg="white", main = "Normalized Heights")


#############################################################################
# Calculate the CHM from the DTM
#############################################################################
# Digital surface models has as input non-normalized point cloud with absolute elevations (uses the sea level as references)
# Canopy height model input is a normalized point cloud which means the derived surface represented by the canopy height (vegetation area)

col <- height.colors(25)

chm <- rasterize_canopy(nlas, res =0.05, algorithm =p2r())
# plot(chm, col= col)

# ## CHM with interpolation of empty points
library(geometry)
chm_IEP <- rasterize_canopy(nlas, res = 0.05, p2r(0.2, na.fill = tin()))
plot(chm_IEP, col = col)
print(class(chm_IEP))

library(terra)
chm_IEP <- as(chm_IEP, "Raster")
crs(chm_IEP) <- ("+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs")
plot(chm_IEP,main = "Canopy height model (CHM)")

##############################################################################
# Correct the extent of CHM (TLS)
##############################################################################

# setting an extent to tls raster
aoi <- shapefile("TLS_area.shp")

# check for the extent
aoi
chm_IEP

# setting same extent to chm_IEP
bb <- extent(12.30909, 12.3099, 51.36535, 51.36598) ## this is the extent of the aoi
extent(chm_IEP) <- bb
chm_IEP <- setExtent(chm_IEP, bb, keepres = T, snap = T)

# plot(chm_IEP)
# chm_IEP
# hist(chm_IEP)

# chm_IEP_aoi <- writeRaster(chm_IEP,'chm_IEP.tif',options=c('TFW=YES'), overwrite=TRUE) #if you want to save

##############################################################################
## Calculate standard deviation of CHM (TLS)
##############################################################################
h_tls_sd = aggregate(chm_IEP, fact = 10, fun = "sd")
plot(h_tls_sd,main = "Standard Deviation of the CHM")
hist(h_tls_sd, xlab = "Normalized Canopy Height")

##############################################################################
## Load the Sentinel-2 data
#############################################################################
h_s2_top <- raster("ETH_GlobalCanopyHeight_DE_CHM_clip.tif") 
plot(h_s2_top)
h_s2_top
h_s2_sd <- raster("ETH_GlobalCanopyHeight_DE_STD_clip.tif")
plot(h_s2_sd)

##############################################################################
# Reproject the TLS raster to match the Sentinel product projection
##############################################################################

library(raster)
h_tls_top_repro <- projectRaster(chm_IEP, h_s2_top, alignOnly = F)
h_tls_top_repro
h_tls_sd_repro <- projectRaster(h_tls_sd, h_s2_top, alignOnly = F)
h_tls_sd_repro

# Reproject the TLS raster to match the Sentinel product
# plot(extent(h_s2_top), col="green", main="Comparison of Sentinel 2 and TLS Data Extents", xlab="Longitude", ylab="Latitude")
#plot(extent(chm_IEP), col="red", add=TRUE)

# Crop Sentinel-2 data to TLS extent
#h_s2_top_cropped <-crop(h_s2_top, extent(chm_IEP))
#h_s2_top_cropped

#plot(extent(h_s2_top_clip), col="green", main="Comparison after Cropping", xlab="Longitude", ylab="Latitude")
#plot(extent(chm_IEP), col="red", add=TRUE)

##
library(terra)
h_tls_top_res <- resample(chm_IEP, h_s2_top, method = "bilinear")
h_tls_sd_res <- resample(h_tls_sd, h_s2_top, method = "bilinear")

# Crop to area of interest of the Sentinel-2 layer 
library(rgdal)
h_s2_top_crop<- crop(h_s2_top, aoi)
h_s2_sd_crop<- crop(h_s2_sd, aoi)
h_tls_top_res<- crop(h_tls_top_res, aoi)
h_tls_sd_res<- crop(h_tls_sd_res, aoi)

###############################################################################
# Calculate differences between the TLS and Sentinel products
###############################################################################
# CHM mean
residuals_top <- h_s2_top_crop - h_tls_top_res
residuals_top
# Standard deviation
residuals_sd <- h_s2_sd_crop - h_tls_sd_res
residuals_sd 

###############################################################################
# Visualize the results through plots
###############################################################################

par(mfrow=c(1,3))
plot(h_tls_top_res, main="Canopy Height (TLS)")
plot(residuals_top, main="Canopy Height (m) Residuals")
plot(h_s2_top_crop, main="Canopy Height(m) (Sentinel-2)")

par(mfrow=c(1,3))
plot(h_tls_sd_res, main="SD Canopy Height (TLS)",xlab = "Normalized Canopy Height")
plot(residuals_sd, main="SD Canopy Height (m) Residuals",xlab = "Normalized Canopy Height")
plot(h_s2_sd_crop, main="SD Canopy Height(m) (Sentinel-2)",xlab = "Normalized Canopy Height")

par(mfrow=c(1,3))
hist(h_tls_top_res, main="Canopy Height (TLS)",xlab = "Normalized Canopy Height",
     xlim=c(0,35))
hist(residuals_top, main="Canopy Height (m) Residuals",xlab = "Normalized Canopy Height",
     xlim=c(0,35))
hist(h_s2_top_crop, main="Canopy Height(m) (Sentinel-2)",xlab = "Normalized Canopy Height",
     xlim=c(0,35))

par(mfrow=c(1,3))
hist(h_tls_sd_res, main="SD Canopy Height (TLS)",xlab = "Normalized Canopy Height",
     xlim=c(0,7))
hist(residuals_sd, main="SD Canopy Height (m) Residuals",xlab = "Normalized Canopy Height",
     xlim=c(0,7))
hist(h_s2_sd_crop, main="SD Canopy Height(m) (Sentinel-2)",xlab = "Normalized Canopy Height",
     xlim=c(0,7))

###############################################################################
# Aggregate TLS data and calculate statistical measures (mean, minimum, max-
# imum, standard deviation). Do the same for Sentinel-2 data and Residuals data
###############################################################################

## ## Mean TLS CHM & SD
ras_final <- aggregate(chm_IEP, fact= 10,fun=max )
res_mean <-cellStats(ras_final, 'mean')
res_sd <-cellStats(ras_final, 'sd')

## Mean Sentinel-2 CHM & SD
ras_final2 <- aggregate(h_s2_top, fact= 10,fun=max )
res_mean2 <-cellStats(h_s2_top, 'mean')
res_sd2 <-cellStats(h_s2_top, 'sd')

# Mean residuals CHM & SD
mean_residuals_top <- cellStats(residuals_top, stat='mean', na.rm=TRUE)
mean_residuals_sd <- cellStats(residuals_sd, stat='mean', na.rm=TRUE)

#end
###############################################################################