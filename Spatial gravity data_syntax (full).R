###############################################
## Spatial analysis of bilateral trade flows ##
###############################################

# R version 3.2.3 Patched (2016-02-15 r70179)

# Clear workspace in RStudio
rm(list=ls())


####################################
# Install (patched) "splm" package #
####################################

## Option 1: Install patched "splm" package from local source using Terminal
# (1) Open Terminal
# (2) Navigate to folder where "splm_1.3-7.tar.gz" is saved, e.g., by running: cd /Users/username/Desktop (changes working directory to desktop; use pwd command to print current working directory)
# (3) Run the following command: R CMD INSTALL splm_1.3-7.tar.gz

## Option 2: Install patched "splm" package from local source using R
# Save source package "splm_1.3-7.tar.gz" to, e.g., desktop and then run:
# R> install.packages("~/Desktop/splm_1.3-7.tar.gz", repos = NULL, type="source")

## Note: After installing the patched package from local source, DO NOT run R> install.packages("splm") as this will override the patched version with the original version!


###################################
# Install and load other packages #
###################################

# Install packages
install.packages("countrycode")
install.packages("fBasics")
install.packages("foreign")
install.packages("gdata")
install.packages("ggplot2")
install.packages("maptools")
install.packages("Matrix")
install.packages("RANN")
install.packages("raster")
install.packages("reshape2")
install.packages("rgdal")
install.packages("rgeos")
install.packages("rworldmap")
install.packages("sp")
install.packages("spdep")
install.packages("tripack")

# Load packages
library(countrycode)
library(fBasics)
library(foreign)
library(gdata)
library(ggplot2)
library(maptools)
library(Matrix)
library(RANN)
library(raster)
library(reshape2)
library(rgdal)
library(rgeos)
library(rworldmap)
library(sp)
library(spdep)
library(splm)
library(tripack)

# Create folder and set working directory
dir.create("~/Desktop/Spatial gravity model")
setwd("~/Desktop/Spatial gravity model")


##########################
# Dyadic data management #
##########################

# Import CEPII gravity data (source: http://www.cepii.fr/CEPII/en/bdd_modele/presentation.asp?id=8; last accessed on 2016-03-20)
#gravity <- read.dta("gravdata_cepii.dta") 

# Alternatively: Download, unzip and import gravity data directly from CEPII webpage
temp <- tempfile()
download.file("http://www.cepii.fr/anglaisgraph/bdd/gravity/gravdata_cepii.zip", temp, mode="w")
unzip(temp, "gravdata_cepii.dta")
gravity <- read.dta("gravdata_cepii.dta")
unlink(temp)

# Import CEPII distance data (source: http://www.cepii.fr/CEPII/en/bdd_modele/presentation.asp?id=6; last accessed on 2016-03-20)
#geo_distances <- read.dta("dist_cepii.dta") 

# Alternatively: Download and import distance data directly from CEPII webpage
download.file("http://www.cepii.fr/distance/dist_cepii.dta", "~/Desktop/Spatial gravity model/dist_cepii.dta", mode="w")
distance_geo <- read.dta("dist_cepii.dta")

# Import CoW trade data (source: http://correlatesofwar.org/data-sets/bilateral-trade; last accessed on 2016-03-20)
#trade <- read.csv("dyadic_trade_3.0.csv", sep=",", header=T, fill=T, quote="", row.names=NULL, stringsAsFactors=F)

# Alternatively: Download, unzip and import trade data directly from CoW webpage
temp2 <- tempfile()
download.file("http://correlatesofwar.org/data-sets/bilateral-trade/cow_trade_3.0/at_download/file", temp2, mode="w")
unzip(temp2)
setwd("~/Desktop/Spatial gravity model/COW_Trade_3.0")
trade <- read.csv("dyadic_trade_3.0.csv", sep=",", header=T, fill=T, quote="", row.names=NULL, stringsAsFactors=F)
setwd("~/Desktop/Spatial gravity model")
unlink(temp2)

# Extract distance variables from distance data frame
distance_subset <- distance_geo[, c(1:2, 12)]
colnames(distance_subset) <- c("iso3_o", "iso3_d", "distcap") # rename columns for merging

# Merge geodesic distances with gravity data frame
gravity <- merge(gravity, distance_subset, by=c("iso3_o","iso3_d"), all.x=T, all.y=F) 

# Convert country codes from CoW to ISO3 
ccode1 <- trade$ccode1 
iso3_o <- countrycode(ccode1, "cown", "iso3c")
ccode2 <- trade$ccode2
iso3_d <- countrycode(ccode2, "cown", "iso3c")
trade$iso3_o <- iso3_o
trade$iso3_d <- iso3_d

# Subset trade data frame
trade_subset <- trade[, c(15:16, 3, 6:9)]
trade_subset_t <- subset(trade_subset, year >= 2002 & year <= 2006)

# Set missing trade flows (source==-9) to NA
trade_subset_t$flow1[trade_subset_t$source1==-9] = NA
trade_subset_t$flow2[trade_subset_t$source2==-9] = NA

# Create dyadic trade data frames
trade_od <- trade_subset_t[, c(1:3, 5)]
colnames(trade_od) <- c("iso3_o", "iso3_d", "year", "trade")
trade_do <- trade_subset_t[, c(2, 1, 3:4)]
colnames(trade_do) <- c("iso3_o", "iso3_d", "year", "trade") 

# Subset gravity data frame
gravity_subset <- gravity[, c(1:5, 10:11, 15:16, 32, 39)]
gravity_subset_t <- subset(gravity_subset, year >= 2002 & year <= 2006)

# Merge gravity with trade data frames (note: append non-matching rows of gravity df to the resulting df)
dyad_data <- merge(gravity_subset_t, trade_od, by=c("iso3_o","iso3_d", "year"), all.x=T, all.y=F) 
dyad_data <- merge(dyad_data, trade_do, by=c("iso3_o","iso3_d", "year"), all.x=T, all.y=F) 
dyad_data$trade <- rowMeans(dyad_data[, c("trade.x", "trade.y")], na.rm=T)
dyad_data$trade.x <- NULL
dyad_data$trade.y <- NULL

# Create dyad identifiers
dyad_data$countrypair <- paste(dyad_data$iso3_o, dyad_data$iso3_d, sep = ":") # character
dyad_id_chr <- with(dyad_data, paste(iso3_o, iso3_d))
dyad_data <- within(dyad_data, dyad_id <- match(dyad_id_chr, unique(dyad_id_chr))) # numeric

# Remove observations where origin == destination
dyad_data <- dyad_data[!(dyad_data$iso3_o==dyad_data$iso3_d),]

# Create log-transformed variables
logged_vars <- c("trade", "gdp_o", "gdp_d", "pop_o", "pop_d", "distcap")
dyad_data[logged_vars] <- log(dyad_data[logged_vars])

# Replace -Inf with NA for log-transformed variables
is.na(dyad_data) <- do.call(cbind,lapply(dyad_data, is.infinite))

# Transform unbalanced panel to balanced panel
dyads_incomplete <- unique(dyad_data$countrypair[!complete.cases(dyad_data)])
dyad_data_balanced <- dyad_data[!(dyad_data$countrypair %in% dyads_incomplete),]


###########################
# Spatial data management #
###########################

# Retrieve shapefile from rworldmap package
data(countriesCoarse) 
map <- countriesCoarse

# Extract geographic information from SpatialPolygonsDataFrame
#regions_o <- data.frame(map$ISO_A3, map$GEO3, map$GEO3major, map$continent)
#colnames(regions_o) <- c("iso3_o", "geo3_o", "geo3major_o", "continent_o")
#regions_d <- regions_o # duplicate for subsequent merging
#colnames(regions_d) <- c("iso3_d", "geo3_d", "geo3major_d", "continent_d")

# Merge geographic information with gravity data frames (note: append non-matching rows of gravity df to the resulting df)
#dyad_data <- merge(dyad_data, regions_d, by=c("iso3_d"), all.x=T, all.y=F) 
#dyad_data <- merge(dyad_data, regions_o, by=c("iso3_o"), all.x=T, all.y=F) 

# Important: Reorder dataframe to match connectivity matrices
#dyad_data <- dyad_data[with(dyad_data, order(iso3_o, iso3_d, year)),]

# Build country list based on intersections
countries_o <- unique(dyad_data_balanced$iso3_o) # extract list of countries
countries_d <- unique(dyad_data_balanced$iso3_d)
countries_p <- unique(map$ISO_A3)
countries_od <- intersect(countries_o, countries_d)
country_list <- intersect(countries_od, countries_p)
country_list_full <- country_list # duplicate for plot

# Replicate sample of Baier & Bergstrand (2007)
bb <- c("AGO", "ALB", "ARE", "ARG", "AUS", "AUT", "BEL", "BFA", "BGD", "BGR", "BOL", "BRA",
        "CAN", "CHE", "CHL", "CHN", "CIV", "CMR", "COD", "COG", "COL", "CRI", "CYP", "DEU", "DNK",
        "DOM", "DZA", "ECU", "EGY", "ESP", "ETH", "FIN", "FRA", "GAB", "GBR", "GHA", "GMB", "GRC",
        "GTM", "GUY", "HKG", "HND", "HTI", "HUN", "IDN", "IND", "IRL", "IRN", "ISR", "ITA", "JAM",
        "JPN", "KEN", "KOR", "LKA", "LUX", "MAR", "MDG", "MEX", "MLI", "MOZ", "MRT", "MUS", "MWI",
        "MYS", "NER", "NGA", "NIC", "NLD", "NOR", "NZL", "PAK", "PAN", "PER", "PHL", "POL", "PRT",
        "PRY", "ROU", "SAU", "SDN", "SEN", "SGP", "SLE", "SLV", "SWE", "SYR", "THA", "TTO", "TUN",
        "TUR", "UGA", "URY", "USA", "VEN", "ZMB", "ZWE")

countries_region <- c("AND", "AUT", "BEL", "CHI", "CYP", "CZE", "DNK", "EST",
            "FRO", "FIN", "FRA", "DEU", "GIB", "GRC", "GRL", "HUN", "ISL", "IRL", "IMY",
            "ITA", "LVA", "LIE", "LTU", "LUX", "MLT", "MCO", "NLD", "NOR", "POL", "PRT",
            "SMR", "SVK", "SVN", "ESP", "SWE", "CHE", "GBR", "ALB", "BLR", "BIH", "BGR",
            "HRV", "MKD", "MDA", "MNE", "ROU", "RUS", "SRB", "TUR", "UKR", "BLZ", "CRI",
            "SLV", "GTM", "HND", "NIC", "PAN", "CAN", "MEX", "USA", "ARG", "BOL", "BRA",
            "CHL", "COL", "ECU", "FLK", "GUF", "GUY", "PRY", "PER", "SUR", "URY", "VEN")

# Build country list based on selected GEO3 regions from SpatialPolygonsDataFrame 
#subset_region <- map[which(map$GEO3major=="Europe" | map$GEO3major=="North America" | map$GEO3major=="Latin America and the Caribbean Polar"), ]
#subset_region$ISO_A3 <- factor(subset_region$ISO_A3) # match factor levels 
#countries_region <- unique(subset_region$ISO_A3)

# Build final country list based on intersections
## Note: Sample can easily be changed by adjusting country_list (no further modifications of syntax required)
country_list_subset <- intersect(bb, countries_region)
country_list <- intersect(country_list, country_list_subset)

# Subset spatial polygons df and dyad df to match country list
map_subset <- map[map$ISO_A3 %in% country_list, ] 
map_subset$ISO_A3 <- factor(map_subset$ISO_A3)
dyad_data_subset <- dyad_data_balanced[dyad_data_balanced$iso3_o %in% country_list & dyad_data_balanced$iso3_d %in% country_list, ]

# Plot geographical distribution of countries (full sample)
data(countriesCoarseLessIslands) # retrieve map for subsequent plots
map_plot <- countriesCoarseLessIslands
map_plot <- map_plot[map_plot$ADMIN!="Antarctica", ] # remove Antarctica polygon for aesthetic reasons

map_full <- map[map$ISO_A3 %in% country_list_full, ] 
map_full$ISO_A3 <- factor(map_full$ISO_A3)

ggplot() + 
  geom_polygon(data=map_plot, aes(x=long, y=lat, group=group)) + 
  geom_polygon(data=map_full, aes(x=long, y=lat, group=group), color="white", fill="white", alpha=0.3) +
  labs(x="", y="") +
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +   
  coord_equal()

# Plot geographical distribution of countries (subset with full sample)
ggplot() + 
  geom_polygon(data=map_plot, aes(x=long, y=lat, group=group)) + 
  geom_polygon(data=map_full, aes(x=long, y=lat, group=group), color="white", fill="white", alpha=0.3) +
  geom_polygon(data=map_subset, aes(x=long, y=lat, group=group), color="red", fill="red", alpha=0.4) +
  labs(x="", y="") +
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +   
  coord_equal()

# Create matrix of polygon centroids
coord_map <- coordinates(map_subset) 
ID <- map_subset$ISO_A3 

# Build neighbor list based on sphere of influence
tri_connect <- tri2nb(coord_map, row.names=ID)
connect <- graph2nb(soi.graph(tri_connect, coord_map), row.names=ID)

# Visualize neighbor list
plot(map_plot, border="grey30", col="grey90", xlab="", ylab="", axes=F) 
plot(map_subset, add=T, border="grey30", col="grey60", xlab="", ylab="", axes=F)
plot(connect, coord_map, add=T, pch=19, cex=0.9, col="red")

# Transform neighbor list to n x n connectivity matrix
connect_mat <- nb2mat(connect, style="B", zero.policy=T) #  B = binary weights
colnames(connect_mat) <- rownames(connect_mat)


#############################################
# Create dyadic spatial dependence matrices #
#############################################

## Warning: This might take a while to run!

# Create n x n identity matrix based on connectivity matrix
dim_conn <- dim(connect_mat) # retrieve n
connect_mat_id <- diag(1,dim_conn) 
dimnames(connect_mat_id) <- dimnames(connect_mat)

# Create N x N origin, destination and dyad connectivity matrices 
w_origin <- kronecker(connect_mat, connect_mat_id, make.dimnames = T) # origin dependence spatial weight matrix
w_destination <- kronecker(connect_mat_id, connect_mat, make.dimnames = T) # destination matrix
w_dyad <- kronecker(connect_mat, connect_mat, make.dimnames = T) # origin/destination matrix

# Remove redundant country pairs in N x N connectivity matrices 
countrypair_list <- unique(dyad_data_subset$countrypair) # extract country pairs in dyadic df
w_origin_final <- w_origin[rownames(w_origin) %in% countrypair_list, colnames(w_origin) %in% countrypair_list]
w_destination_final <- w_destination[rownames(w_destination) %in% countrypair_list, colnames(w_destination) %in% countrypair_list]
w_dyad_final <- w_dyad[rownames(w_dyad) %in% countrypair_list, colnames(w_dyad) %in% countrypair_list]

# Create origin- + destination matrix (sum)
w_dyad_sum <- 0.5*(w_origin_final + w_destination_final)

# Row-standardize N x N matrices
get.ZeroPolicyOption()
set.ZeroPolicyOption(TRUE)
w_origin_final_rs <- w_origin_final / apply(w_origin_final, 1, sum, zero.policy=T)
w_destination_final_rs <- w_destination_final / apply(w_destination_final, 1, sum, zero.policy=T)
w_dyad_sum_rs <- w_dyad_sum / apply(w_dyad_sum, 1, sum, zero.policy=T)
w_dyad_final_rs <- w_dyad_final / apply(w_dyad_final, 1, sum, zero.policy=T)

# Create neighborhood lists from row-standardized matrices
w_origin_list <- mat2listw(w_origin_final_rs)
w_destination_list <- mat2listw(w_destination_final_rs)
w_dyad_sum_list <- mat2listw(w_dyad_sum_rs)
w_dyad_list <- mat2listw(w_dyad_final_rs)


##########################################
# Moran's I  for spatial autocorrelation #
##########################################

## Warning: This might take a while to run!

# Clear workspace for subsequent analyses
rm(list= ls()[!(ls() %in% c("dyad_data_subset", "w_origin_list", "w_destination_list", "w_dyad_sum_list", "w_dyad_list"))])

# 2002
data_2002 <- subset(dyad_data_subset, year==2002)
moran.test(data_2002$trade, listw=w_origin_list, zero.policy=T) 
moran.test(data_2002$trade, listw=w_destination_list, zero.policy=T) 
moran.test(data_2002$trade, listw=w_dyad_sum_list, zero.policy=T) 
moran.test(data_2002$trade, listw=w_dyad_list, zero.policy=T) 

# 2003
data_2003 <- subset(dyad_data_subset, year==2003)
moran.test(data_2003$trade, listw=w_origin_list, zero.policy=T) 
moran.test(data_2003$trade, listw=w_destination_list, zero.policy=T) 
moran.test(data_2003$trade, listw=w_dyad_sum_list, zero.policy=T) 
moran.test(data_2003$trade, listw=w_dyad_list, zero.policy=T) 

# 2004
data_2004 <- subset(dyad_data_subset, year==2004)
moran.test(data_2004$trade, listw=w_origin_list, zero.policy=T) 
moran.test(data_2004$trade, listw=w_destination_list, zero.policy=T) 
moran.test(data_2004$trade, listw=w_dyad_sum_list, zero.policy=T)
moran.test(data_2004$trade, listw=w_dyad_list, zero.policy=T) 

# 2005
data_2005 <- subset(dyad_data_subset, year==2005)
moran.test(data_2005$trade, listw=w_origin_list, zero.policy=T) 
moran.test(data_2005$trade, listw=w_destination_list, zero.policy=T) 
moran.test(data_2005$trade, listw=w_dyad_sum_list, zero.policy=T) 
moran.test(data_2005$trade, listw=w_dyad_list, zero.policy=T) 

# 2006
data_2006 <- subset(dyad_data_subset, year==2006)
moran.test(data_2006$trade, listw=w_origin_list, zero.policy=T) 
moran.test(data_2006$trade, listw=w_destination_list, zero.policy=T) 
moran.test(data_2006$trade, listw=w_dyad_sum_list, zero.policy=T) 
moran.test(data_2006$trade, listw=w_dyad_list, zero.policy=T) 


#############################
# Gravity regression models #
#############################

## Warning: This definitely takes a while to run!

# Function to calculate AIC for objects of class splm
aic_spml <- function(object, k=2) {
  s <- summary(object)
  ll <- s$logLik
  npar <- length(coef(s))
  N <- nrow(s$model)
  aic <- -2*ll+k*npar
  names(aic) <- "AIC"
  return(aic)
}

# Function to calculate % effect of binary variables in log-transformed models
effect_binary <- function(object) {
  rta <- object$coefficients["rta"]
  p <- 100*(exp(rta)-1)
  names(p) <- "effect in %"
  return(p)
}

# Basic model specifications 
fm_re <- trade ~ gdp_o + gdp_d + distcap + contig + comlang_off + rta + factor(year) 
fm_fe <- trade ~ gdp_o + gdp_d + rta + factor(year) 


#######################################
# Random effects spatial panel models #
#######################################

# (1) Origin-based dependence
re_sar_origin <- spml(fm_re, dyad_data_subset, index=c("dyad_id", "year"), listw=w_origin_list, model="random", spatial.error="none", lag=T)
summary(re_sar_origin)

# (2) Destination-based dependence
re_sar_destination <- spml(fm_re, dyad_data_subset, index=c("dyad_id", "year"), listw=w_destination_list, model="random", spatial.error="none", lag=T)
summary(re_sar_destination)

# (3) Origin- + destination-based dependence (sum)
re_sar_dyad_sum <- spml(fm_re, dyad_data_subset, index=c("dyad_id", "year"), listw=w_dyad_sum_list, model="random", spatial.error="none", lag=T)
summary(re_sar_dyad_sum)

# (4) Origin-destination-based dependence
re_sar_dyad <- spml(fm_re, dyad_data_subset, index=c("dyad_id", "year"), listw=w_dyad_list, model="random", spatial.error="none", lag=T)
summary(re_sar_dyad)

# Calculate % effects of RTAs
effect_binary(re_sar_origin)
effect_binary(re_sar_destination)
effect_binary(re_sar_dyad_sum)
effect_binary(re_sar_dyad)


######################################
# Fixed effects spatial panel models #
######################################

# (1) Origin-based dependence
fe_sar_origin <- spml(fm_fe, dyad_data_subset, index=c("dyad_id", "year"), listw=w_origin_list, model="within", spatial.error="none", lag=T)
summary(fe_sar_origin) 

# (2) Destination-based dependence
fe_sar_destination <- spml(fm_fe, dyad_data_subset, index=c("dyad_id", "year"), listw=w_destination_list, model="within", spatial.error="none", lag=T)
summary(fe_sar_destination)

# (3) Origin- + destination-based dependence (sum)
fe_sar_dyad_sum <- spml(fm_fe, dyad_data_subset, index=c("dyad_id", "year"), listw=w_dyad_sum_list, model="within", spatial.error="none", lag=T)
summary(fe_sar_dyad_sum)

# (4) Origin-destination-based dependence
fe_sar_dyad <- spml(fm_fe, dyad_data_subset, index=c("dyad_id", "year"), listw=w_dyad_list, model="within", spatial.error="none", lag=T)
summary(fe_sar_dyad)

# Calculate % effects of RTAs
effect_binary(fe_sar_origin)
effect_binary(fe_sar_destination)
effect_binary(fe_sar_dyad_sum)
effect_binary(fe_sar_dyad)


#######################################################
# Fixed effects spatial panel models: Lee & Yu (2010) #
#######################################################

# (1) Origin-based dependence
fe_sar_origin_ly <- spml(fm_fe, dyad_data_subset, index=c("dyad_id", "year"), listw=w_origin_list, model="within", spatial.error="none", lag=T, LeeYu=T, Hess=F)
summary(fe_sar_origin_ly)

# (2) Destination-based dependence
fe_sar_destination_ly <- spml(fm_fe, dyad_data_subset, index=c("dyad_id", "year"), listw=w_destination_list, model="within", spatial.error="none", lag=T, LeeYu=T, Hess=F)
summary(fe_sar_destination_ly)

# (3) Origin- + destination-based dependence (sum)
fe_sar_dyad_sum_ly <- spml(fm_fe, dyad_data_subset, index=c("dyad_id", "year"), listw=w_dyad_sum_list, model="within", spatial.error="none", lag=T, LeeYu=T, Hess=F)
summary(fe_sar_dyad_sum_ly)

# (4) Origin-destination-based dependence
fe_sar_dyad_ly <- spml(fm_fe, dyad_data_subset, index=c("dyad_id", "year"), listw=w_dyad_list, model="within", spatial.error="none", lag=T, LeeYu=T, Hess=F)
summary(fe_sar_dyad_ly)


####################
# Impact estimates #
####################

# Retrieve dimension of T
time <- length(unique(dyad_data_subset$year))

# Modify impacts() function to allow calculation of impacts for style=M of weights list object
## Note: Temporary modification of function is necessary as row-standardization was carried out manually and, when converting to weights list object, style is not retrieved automatically by spdep/mat2listw() but set to "M" instead).
impacts.splm <- function(obj, listw = NULL, time = NULL, ..., tr=NULL, R=200, type="mult", empirical=FALSE, Q=NULL){
  
  if(is.null(listw) && is.null(tr)) stop("either listw or tr should be provided")
  
  if(!is.null(listw) ){	
    if(is.null(time) && is.null(tr)) stop("time periods should be provided")
  }
  if(is.null(tr)){
    sparse.W <- listw2dgCMatrix(listw)
    s.lws <- kronecker(Diagonal(time) , sparse.W)
    tr <- trW(s.lws, type= type)
  }
  if(is.na(match(obj$type, c("fixed effects lag","fixed effects sarar","random effects ML", "fixed effects GM","lag GM","fixed effects GM")))) stop("object type not recognized")
  
  if(obj$type == "fixed effects lag"){
    
    class(obj)<- "gmsar"	
    obj$type <- "SARAR"
    obj$data <- as.vector(obj$model)
    obj$s2 <- obj$sigma2
    obj$secstep_var <- obj$vcov
    imp <- impacts(obj, tr=tr, R=R, ...)
  }
  
  if(obj$type == "fixed effects sarar"){
    
    class(obj)<- "gmsar"	
    obj$type <- "SARAR"
    rho <- obj$coefficients[2]
    obj$coefficients <- obj$coefficients[-2]
    obj$data <- as.vector(obj$model)
    obj$s2 <- obj$sigma2
    obj$secstep_var <- obj$vcov[-2,-2]
    imp <- impacts(obj, tr=tr, R=R,...)		
  }
  if(obj$type == "fixed effects error") stop("Impacts Estimates are not available for Error Model")
  
  if(obj$type == "random effects ML")	{
    
    if(!is.null(obj$arcoef)) {
      class(obj)<- "gmsar"	
      obj$type <- "SARAR"
      
      obj$coefficients <- c(obj$arcoef, obj$coefficients)
      obj$data <- as.vector(obj$model)
      obj$s2 <- obj$sigma2
      obj$secstep_var <- matrix(0,nrow(obj$vcov)+1,nrow(obj$vcov)+1)
      obj$secstep_var[1,1] <- obj$vcov.arcoef
      obj$secstep_var[(2:(nrow(obj$vcov)+1)),(2:(nrow(obj$vcov)+1))] <- obj$vcov
      imp <- impacts(obj, tr=tr, R=R, ...)		
    }
    else stop("Impacts Estimates are not available for Error Model")		
  }
  
  if(obj$type == "fixed effects GM"){
    
    if(is.null(obj$endog)) {
      obj$secstep_var <- vcov(obj)			
      class(obj)<- "gmsar"	
      obj$type <- "SARAR"
      obj$data <- as.vector(obj$model)
      obj$s2 <- obj$sigma2
      
      imp <- impacts(obj, tr=tr, R=R, ...)		
    }
    
    else stop("No impacts estimates when endogenous variables are present in the system")
  }
  
  if(obj$type == "lag GM")			{
    
    if(is.null(obj$endog)) {
      
      class(obj)<- "gmsar"	
      obj$type <- "SARAR"
      obj$secstep_var <- obj$var			
      obj$data <- as.vector(obj$model)
      obj$s2 <- obj$sigma2
      
      imp <- impacts(obj, tr=tr, R=R, ...)		
    }
    else stop("No impacts estimates when endogenous variables are present in the system")
  
  }
 
  if(obj$type == "random effects GM")			{
    
    if(is.null(obj$endog)) {
      
      class(obj)<- "gmsar"	
      obj$type <- "SARAR"
      obj$secstep_var <- obj$vcov			
      obj$data <- as.vector(obj$model)
      obj$s2 <- obj$sigma2
      
      imp <- impacts(obj, tr=tr, R=R, ...)		
    }
    else stop("No impacts estimates when endogenous variables are present in the system")
  }
  return(imp)	
}


## Calculate impacts for random effects models

# (1) Origin-based dependence
set.seed(10000)
imp_re_origin <- impacts(re_sar_origin, listw=w_origin_list, time=time)
summary(imp_re_origin, zstats=T, short=T)

# (2) Destination-based dependence
set.seed(10000)
imp_re_destination <- impacts(re_sar_destination, listw=w_destination_list, time=time)
summary(imp_re_destination, zstats=T, short=T)

# (3) Origin- + destination-based dependence (sum)
set.seed(10000)
imp_re_dyad_sum <- impacts(re_sar_dyad_sum, listw=w_dyad_sum_list, time=time)
summary(imp_re_dyad_sum, zstats=T, short=T)

# (4) Origin-destination-based dependence
set.seed(10000)
imp_re_dyad <- impacts(re_sar_dyad, listw=w_dyad_list, time=time)
summary(imp_re_dyad, zstats=T, short=T)


## Calculate impacts for fixed effects models

# (1) Origin-based dependence
set.seed(10000)
imp_fe_origin <- impacts(fe_sar_origin, listw=w_origin_list, time=time)
summary(imp_fe_origin, zstats=T, short=T)

# (2) Destination-based dependence
set.seed(10000)
imp_fe_destination <- impacts(fe_sar_destination, listw=w_destination_list, time=time)
summary(imp_fe_destination, zstats=T, short=T)

# (3) Origin- + destination-based dependence (sum)
set.seed(10000)
imp_fe_dyad_sum <- impacts(fe_sar_dyad_sum, listw=w_dyad_sum_list, time=time)
summary(imp_fe_dyad_sum, zstats=T, short=T)

# (4) Origin-destination-based dependence
set.seed(10000)
imp_fe_dyad <- impacts(fe_sar_dyad, listw=w_dyad_list, time=time)
summary(imp_fe_dyad, zstats=T, short=T)


####################
# Model validation #
####################

# Lagrange multiplier tests
## Important: Note that alternative hypothesis seems to be one of NO random effects (other than suggested by the summary output)!
bsktest(fm_re, dyad_data_subset,index=c("dyad_id", "year"), listw=w_origin_list, test="LM1", standardize=T) 
## Important: Note that alternative hypothesis seems to be one of NO spatial autocorrelation (other than suggested by the summary output)!
bsktest(fm_re, dyad_data_subset,index=c("dyad_id", "year"), listw=w_origin_list, test="LM2", standardize=T) 

# Function to fix sphtest() for spatial Hausman test (adapted from http://lists.r-forge.r-project.org/pipermail/splm-commits/2015-November/000202.html)
## Note: Treat results from sphtest() with a degree of caution.
sphtest_modified <- function (x, ...)
{
  UseMethod("sphtest")
}

sphtest.formula <- function (x, data, index = NULL, listw,
                             spatial.model = c("lag", "error", "sarar"),
                             method = c("ML", "GM"), errors = c("KKP", "BSK"),...) {
  switch(match.arg(spatial.model),
         lag = {
           lag = TRUE
           spatial.error = FALSE
         },
         error = {
           lag = FALSE
           spatial.error = TRUE
         },
         sarar = {
           lag = TRUE
           spatial.error = TRUE
         })
  
  errors <- match.arg(errors)
  x0 <- update(x, .~.-1)
  
  method <- switch(match.arg(method),
                   ML = {
                     spatial.error <- if(spatial.error) {
                       spatial.error <- if(errors=="BSK") "b" else "kkp"
                     } else {
                       spatial.error <- "none"
                     }
                     femod <- spml(x, data = data, index = index, listw = listw, lag = lag,
                                   spatial.error = spatial.error, model = "within")
                     remod <- spml(x, data = data, index = index, listw = listw, lag = lag,
                                   spatial.error = spatial.error, model = "random")
                   },
                   GM = {
                     femod <- spgm(x, data = data, index = index, listw = listw, lag = lag,
                                   spatial.error = spatial.error, model = "within", moments = "fullweights")
                     remod <- spgm(x, data = data, index = index, listw = listw, lag = lag,
                                   spatial.error = spatial.error, model = "random", moments = "fullweights")
                   },
                   stop("\n Unknown method"))
  return(sphtest(femod, remod, ...))
}

sphtest.splm <- function (x, x2, ...){
  is.gm <- !is.null(x$ef.sph)
  
  if(is.gm) {
    if (!all.equal(x$legacy, x2$legacy)) stop("The models are different")
    if(x$ef.sph == x2$ef.sph) stop("Effects should be different")
    
    ran <- match("random", c(x$ef.sph, x2$ef.sph))
    if(ran == 1){
      xwith <- x2
      xbetw <- x
    }
    if(ran == 2){
      xwith <- x
      xbetw <- x2
    }
    
    tc <- match(names(coef(xwith)), names(coef(xbetw)) )
    coef.wi <- coef(xwith)
    coef.re <- coef(xbetw)[tc]
    vcov.wi <- xwith$vcov
    vcov.re <- xbetw$vcov[tc,tc]
    
  } else {
    if(is.null(dimnames(x$vcov))) {
      xwith <- x
      xbetw <- x2
    } else {
      xwith <- x2
      xbetw <- x
    }
    tc <- intersect(names(coef(xwith)), names(coef(xbetw)))
    wtc <- match(tc, names(coef(xwith)))
    coef.wi <- coef(xwith)[wtc]
    coef.re <- coef(xbetw)[tc]
    vcov.wi <- xwith$vcov[wtc,wtc]
    vcov.re <- xbetw$vcov[tc,tc]
  }
  dbeta <- coef.wi - coef.re
  df <- length(dbeta)
  dvcov <- vcov.re - vcov.wi
  stat <- abs(t(dbeta) %*% solve(dvcov) %*% dbeta)
  pval <- pchisq(stat, df = df, lower.tail = FALSE)
  names(stat) <- "chisq"
  parameter <- df
  names(parameter) <- "df"
  data.name <- paste(deparse(x$call$formula))
  alternative <- "one model is inconsistent"
  res <- list(statistic = stat, p.value = pval, parameter = parameter,
              method = "Spatial Hausman test (patched)",
              data.name = data.name, alternative = alternative)
  class(res) <- "htest"
  return(res)
}

# Spatial Hausman test (RE with spatial lag or error vs. FE with corresponding model specification)
re_sar_dyad_sum <- spml(fm_re, dyad_data_subset, index=c("dyad_id", "year"), listw=w_dyad_sum_list, model="random", spatial.error="none", lag=T)
fe_sar_dyad_sum <- spml(fm_fe, dyad_data_subset, index=c("dyad_id", "year"), listw=w_dyad_sum_list, model="within", spatial.error="none", lag=T)

sphtest_modified(re_sar_dyad_sum, fe_sar_dyad_sum)