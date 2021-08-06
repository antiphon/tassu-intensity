# Using previous processing steps (see below), 
#
# 1. aggregate the Tassu-data to a finite grid for Poisson-regression modelling
# 2. add covariate information.
#
# v1, 3.8.2021: Extracted from "wolf_observations_poisson_regression-v2.6.Rmd" for github friendly approach.
#

library(dplyr)
library(tidyr)
library(sf)
library(raster)

# The file name for the output data-frame ready for glm etc.
agg_out_file <- "processed/celdat-v2.6.rds"


####################################################
# Region of interest:
#
crs1 <- 3067 # use this for all
crs1_c <- paste0("EPSG:", crs1)
finland <- readRDS("../data/finland.rds")
finland <- st_transform(finland, crs1)
# Drop North from geometries. Hard cut, no other preference given.
box1 <- st_bbox(c(xmin=0, xmax=100000000, ymin=0, ymax = max(72e5 + 30000)), crs=st_crs(finland))
finland1 <- st_crop(finland, box1)


#################################
## Prepare Tassu

### Load tassu. Latest.
tassu_file <- dir("../data_not_public/processed/", "tassu", full.names = TRUE) %>% last()
tassu <- readRDS(tassu_file)

# Add info about being inside relevant region.
tloc <- tassu %>% dplyr::select(x=st_x, y=st_y) %>% 
  st_as_sf(coords=c("x", "y"), crs = crs1)
overlap <-  st_contains(finland1, tloc, sparse=FALSE)
tassu <- tassu %>% mutate(inside_fin = overlap[1,], 
                          x = st_x, y = st_y)

###########################################################
# Apply filtering. Important part here!
# 
#
#
tassu1 <- tassu %>% 
  # cut to region
  filter( inside_fin &
            # then filter to only those obs with more than 1 animal
            havainto_id_count > 1 & havainto_id_subidx == 1 &
            # only wolves for this model
            laji == "susi")


##########################################################
# Modifications to variables for modelling. 
# This frame will contain all tassu info for the modeling.
#
tassu_dat1 <- tassu1 %>% transmute( vuosi = as.integer(format(timestamp, format="%Y")),
                                    kk    = as.integer(format(timestamp, format="%m")),
                                    vk    = as.integer(format(timestamp, format="%V")),
                                    julian = as.integer(format(timestamp, format="%j")),
                                    jvk    = 1 + floor( (julian / 367) * 52 ),
                                    timestamp = timestamp %>% format(), # to avoid TZ issues
                                    timestamp_utc = as.POSIXct(timestamp, tz ="UTC"), # in case needed for plots and such
                                    timestamp_y = timestamp_y,
                                    # time is in "wolf years"
                                    wolf_start_year = vuosi - 1 * (kk < 4),
                                    wolf_month      = (kk - 1 - 3) %% 12 + 1,
                                    wolf_year       = sprintf("%s-%s", wolf_start_year, wolf_start_year + 1),
                                    x,y,
                                    laji,
                                    havainto_id, 
                                    havainto_id_count)


# Tassu data ready for aggregation.
saveRDS(tassu_dat1, "processed/tassu_dat1-v1.rds")
#
######################################################################
# Then prepare conditioning set: Territories, GPS and expert knowledge

# 
terr_file <- dir("../data_not_public/processed/", "chull.*prior_polygons", full.names = TRUE) %>% last()
terr0 <- readRDS(terr_file)
st_crs(terr0) <- st_crs(finland1)

# check: 
# terr0 %>% filter(type != "gps") %>% dplyr::select(start.tdate)
# 
############################
# Only relevant to Tassu time-region.

terr <- terr0 %>% filter( start.tdate >= min(tassu1$timestamp) & 
                            end.tdate <= max(tassu1$timestamp) ) %>% 
  mutate(area_m2  = st_area(geometry) %>% units::drop_units(), 
         area_km2 = (area_m2 /1000^2)  )
#
# Territories ready.
# store for plots etc.
saveRDS(terr, "processed/territories-v1.rds")

######################################################################
#
# Then subset tassu to keep reports only on the known wolf territories, a surrogate for "wolf present".
#
# 1. spatial inclusion. 

tloc <- tassu_dat1[,c("x","y")]  %>%  
    st_as_sf(coords=c("x","y"), crs = crs1)
tpot_spat <- st_intersects(st_geometry(terr), tloc)

# 2. Temporal inclusion. Enough to check hit or not for now. 
hit_cid  <-  sapply( 1:nrow(tassu_dat1), 
                     function(i){
                       ty <- tassu_dat1$timestamp[i]
                       pot <- which( ty >= terr$start.tdate & ty <= terr$end.tdate  )
                       if(length(pot) == 0) return(FALSE)
                       pot2 <- sapply( tpot_spat[pot], function(a) i %in% a )
                       if(sum(pot2) == 0) return(FALSE)
                       return(TRUE)
                     }  )
#
nhit <- mean(hit_cid)
cat(sprintf("%3.2f%% tassu hits some territory.\n", nhit*100))

# data that hits, i.e. the "wolf was present" conditioned data:
tassu_dat2 <- tassu_dat1 %>% 
  filter( hit_cid )

# store for plots etc.
saveRDS(tassu_dat2, "processed/tassu_dat2-v1.rds")


########################################################################
#
# Add covariate information
#
# Corine 18 proportions in 1km x 1km grid and digiroad "forest road" (class 12)
#
# Augmented covariates.
corine8 <-  crop(stack("../data/corine/clc2018_fi_reclassed1_1km.grd"), 
                 as_Spatial(finland1))
                 
projection(corine8) <- crs1_c

mask1 <- mask(corine8, as_Spatial(finland1))
#
# Smooth each cell, useful later
corine_km_smooth <- 6000/2
weight <- focalWeight(corine8, d = corine_km_smooth, type="Gauss")
#
cat("Corine smooth radius:", corine_km_smooth/1000, "km\n")
#
#ok do for all layers
corine8sl <- list()
for(l in names(corine8)) {
  corine8sl[[l]] <- focal( corine8[[l]], w=weight, pad = TRUE, na.rm=TRUE)
}
# cut back
corine8s <- mask(stack(corine8sl), as_Spatial(finland1))

# rescale to simplex
corine8s <- calc(corine8s, function(v) v/sum(v, na.rm=TRUE))
names(corine8s) <- names(corine8sl)
#
# Proportions: take log-ratio against largest
#
#
i0 <- which.max(cellStats(corine8s, mean))
# most common is "Sulk metstä"
corine8s_logratio <- calc(corine8s, function(v) log(v/v[i0]) )  %>% 
  dropLayer(i = i0) # drop the constant

##
# digiroad class 12 roads “Ajopolku, maastoajoneuvolla ajettavissa olevat tiet”:
#
droad0 <- crop(stack("../data/digiroad_luokka12_1x1km_2021-04-14.grd"), 
               as_Spatial(finland1))
droad1 <- setNames(droad0, "droad.cl12")
# make sure droad aligned with corine... should be.

###
# Now add the info to the modeling dataframe
#
# locate cells
dat2cell <- cellFromXY(corine8s, tassu_dat2[,c("x","y")] %>% as.matrix())
# combine.
tassu_dat3 <- cbind(tassu_dat2, 
                    c18    = as.data.frame( corine8s         [dat2cell] ),
                    lr.c18 = as.data.frame( corine8s_logratio[dat2cell] ),
                    as.data.frame(droad1 [dat2cell] ) )

message("spatial covariates added.")

# store for plots etc.
saveRDS(corine8s, "processed/corine8s-v1.rds")


# 
###########################################################
# Now we aggregate to cells.

# Getting spatial coverage is easy enough:
sr <- mask(corine8s, as_Spatial( terr ) )
sk <- which( !is.na(sr$Tiet[]) )

# Then each grid cell one is active at certain times, depending on territories.
# make a proper monthly grid. (wolf time not yet relevant)
tgrid <- seq(as.POSIXct("2011-01-01"), as.POSIXct("2021-02-01"), by = "month") %>% format() # add +1 end

# make sure only those inside corine raster
tassu_dat4  <-  tassu_dat3 %>% filter( !is.na( corine8s$Tiet[  cellFromXY( corine8s, cbind(x,y) )  ] )  ) 
nlost <- nrow(tassu_dat3)-nrow(tassu_dat4)
cat(sprintf("Lost %i (%3.1f%%) datapoints due to rasterisation.\n ", nlost, 100 * nlost/nrow(tassu_dat3)))


celdat1 <- tibble() 
# main loop for binning the wolf observations in space and time.
for(i in 1:(length(tgrid)-1) ) { # old school for loop. easier to debug.
  t1 <- tgrid[i]
  t2 <- tgrid[i+1]
  # which territories were active
  on <- terr %>% filter( ! ( start.tdate >= t2 | end.tdate < t1 ) )
  if(nrow(on) == 0) next # no active territories
  # buff up a bit ~ 1pixel
  on <- st_buffer(on, 1000)
  # which spatial cells are active at this time interval
  activer <- mask(corine8s$Tiet, as_Spatial( on ) )
  sidx <-  which( !is.na(activer[]) )  # cell index of actives
  nact <- length(sidx)
  # Gather Covariates
  corine <- cbind(c18   = as.data.frame( corine8s [ sidx] ),
                  lr.c18 = as.data.frame( corine8s_logratio[sidx] ),
                  as.data.frame(droad1[sidx])
  )
  # figure out which territory each pixel is in, needed for determining the offset areas
  gxy <- coordinates(activer)[sidx,] %>% 
    data.frame() %>%
    st_as_sf(coords=c("x", "y"), crs = st_crs(finland1))
  terri <- st_contains(st_geometry(on), gxy) %>% 
    as.data.frame() %>% 
    transmute(ter = row.id, cell = col.id) %>%
    group_by(cell) %>% # for each pixel, store some info
    mutate(area_km2 = on$area_km2[ ter ],
           areasum = sum(area_km2), 
           terr_hit = n()) %>% 
    filter(row_number() == which.max(area_km2)) %>%
    ungroup() %>% arrange(cell) %>%
    mutate(cell = NULL)
  #
  # check...
  #terri %>% mutate(a = st_coordinates(gxy), x=a[,1],y = a[,2], a=NULL) %>% ggplot() + geom_point(aes(x, y, col = areasum), alpha = 0.5)
  #
  # Then count Tassu-wolfs:
  cvecw <- cvec <- rep(0, nact)
  thit <- tassu_dat4 %>% filter( timestamp_utc >= as.POSIXct(t1, tz="UTC") &
                                 timestamp_utc <  as.POSIXct(t2, tz="UTC") )  # temporally possible
  if(nrow(thit)){
    cell_thit <- cellFromXY(corine8s, as.matrix(thit[,c("x","y")]))
    xe <- tibble(k=cell_thit, 
                 cnt =  thit$havainto_id_count)  %>% 
      group_by(k) %>% 
      mutate(cnt = n(), wcnt = sum(cnt)) %>% 
      filter(row_number() == 1) %>% ungroup() %>% arrange(k)
    khit <- match(as.integer(xe$k), sidx) # match the cell numbering
    cvec[khit] <- xe$cnt   # this is sum of "wolf reports"
    cvecw[khit] <- xe$wcnt # this is sum of "reported wolves"
  }
  # compile
  celdat1 <- celdat1 %>% 
    bind_rows(
      tibble(nsusi = cvec, 
             nwsusi = cvecw,
             raster_idx = sidx, 
             start.tdate = t1, 
             end.tdate   = t2, # exclusive
             mean.tdate  = mean( as.POSIXct(c(t1, t2) )) ) %>%
        bind_cols(terri, corine)
    )
}


# Add date etc details
cxy <- coordinates(sr)
celdat1 <- celdat1 %>% mutate(x      = cxy[raster_idx,1], 
                              y      = cxy[raster_idx,2],
                              vuosi  = as.integer(format(mean.tdate, "%Y")),
                              vuosif = factor(vuosi),
                              kk     = as.integer(format(mean.tdate, "%m")),
                              kkf    =  factor(kk),
                              viikko =  as.integer(format(mean.tdate, "%V")),
                              special_may = vuosi >= 2017 & kk %in% 1:3, # for expert territories
                              area_m2 = area_km2 * 1000^2,
                              wolf_start_year = vuosi - 1 * (kk < 4),
                              wolf_month = (kk - 4) %% 12 + 1,
                              wolf_year = sprintf("%s-%s", wolf_start_year, wolf_start_year+1),
                              wolf_monthf = factor(wolf_month, levels = 1:12),
                              wolf_yearf = factor(wolf_year)
                              )
# store
saveRDS(celdat1, agg_out_file)
message("Compiling the modeling data done.")



