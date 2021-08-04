# Convert the GPS-tracking data to guesses of territories, 
# and join with the expert territory information.
#
#
# v6.1: Expert territories active Jan - March, not just March
#
# 3.8.2021: Github friendly, code extracted from "reviirit-v6.1.Rmd"
#
#

library(dplyr)
library(ggplot2)
library(raster)
library(sf)

theme_set(theme_bw())

knitr::opts_chunk$set(echo = TRUE, message=FALSE, fig.width = 10)

finland <- readRDS("../data/finland.rds")
crs <- crs1 <- 3067 # "+init=epsg:3067" # old

# not included.
bad_series <- "pinna4"

# Helpers: (should put this somewhere else)
posix2year <- function(x) {
  y <- as.numeric( format(x, "%Y") )
  i1 <- as.numeric( as.POSIXct(paste0(y, "/01/01 00:00:00"), tz = "UTC") )
  i2 <- as.numeric( as.POSIXct(paste0(y, "/12/31 23:59:59"), tz = "UTC") )
  s <- as.numeric(x)
  pmax(y, y + (s-i1)/(i2-i1))
}

######################################################
## GPS territories: Convex hull on space.
# use latest:
gps_file <- dir("processed/", "pannat.*filtered", full.names = TRUE) %>% last()
panta <- readRDS(gps_file)

# rought filtering, some outliers
panta <- panta %>% filter( x  > 2e5 & y < 74e5 )

# Compute convex hull for each series.
panta <- panta %>% 
  group_by(series) %>% 
  mutate(chull_rank = {
    h <- chull(x,y)
    v <- rep(NA, n())
    v[h] <- 1:length(h)
    v
    } ) %>% 
  mutate(x_centered = x - min(x), 
         y_centered = y - min(y)) %>%  # for plots, centering might be nice.
  ungroup()


# create a new storage for hulls, include time intervals
gpshulls <- lapply( split(panta, panta$series), function(x) {
  cxy <- x %>% filter(!is.na(chull_rank)) %>% 
    dplyr::select(c("x", "y", "chull_rank")) %>% 
    arrange(chull_rank)
  ti <- range(x$timestamp)
  ty <- range(x$timestamp_y)
  bbox <- data.frame(apply(cxy[,1:2], 2, range), tdate=ti, tyear=ty)
  list(xy = as.matrix( cxy[,1:2] ) , bbox = bbox,  type =  "gps")
})

# gps hulls compiled. 
#
###################################################
# Process Expert knowledge territories

st_read <- function( ... ) sf::st_read(..., quiet = TRUE)
pri_reviirit1 <- 
  list(
    st_read("reviirit_raw/Reviirit2017-2020/Reviirit2017/Reviirit2017.dbf") %>% mutate(year = 2017),
    st_read("reviirit_raw/Reviirit2017-2020/2018_reviirit/Lauma_tai_pari.dbf") %>% mutate(year = 2018),
    st_read("reviirit_raw/Reviirit2017-2020/2018_reviirit/Laumat.dbf") %>% mutate(year = 2018),
    st_read("reviirit_raw/Reviirit2017-2020/2018_reviirit/Parit.dbf") %>% mutate(year = 2018),
    st_read("reviirit_raw/Reviirit2017-2020/2018_reviirit/Rajalaumat.dbf") %>% mutate(year = 2018),
    st_read("reviirit_raw/Reviirit2017-2020/Kaikki (Reviirit 2019)/Kaikki.dbf") %>% mutate(year = 2019),
    st_read("reviirit_raw/Reviirit2017-2020/Reviirit2020/Reviirit2020_.dbf") %>% mutate(year = 2020)
  )
# same projection
pri_reviirit <- lapply(pri_reviirit1, st_transform, crs = crs1) %>% 
  lapply(function(v) v[,'year']) %>% 
  bind_rows() %>%
  group_by(year) %>% 
  mutate(id = paste0("pri", year, "_", 1:n())) %>% 
  ungroup()

# 
pri_hulls <- lapply( 2017:2020, function(y){
  d <- pri_reviirit %>% filter(year == y)
  times <- data.frame(tdate = c(sprintf("%d-01-01", y),
                                sprintf("%d-03-31", y))) %>% 
                        mutate(tyear = posix2year(as.POSIXct(tdate) ))
  o <- lapply(1:nrow(d), function(i){
    xy <- st_coordinates(d[i,])[,1:2]
    colnames(xy) <- c("x", "y")
    # keep form, not closed polygon
    xy <- xy[ -nrow(xy), ]
    list(xy = xy, bbox = cbind( apply(xy, 2, range),  times), type = "prior")
  })
  names(o) <- d$id
  o
} ) %>% do.call(what = c)

# Expert hulls compiled.
###########
# Join.
cthulls <- c(gpshulls, pri_hulls)
# ok.


###################################
#  Unify, convert all to sf-class.
topes <- cthulls[! names(cthulls) %in% bad_series ]

topel <- lapply(names(topes), function(s) {
  x <- topes[[s]]
  xy <- rbind(x$xy, x$xy[1,])
  # round of to 10 meters
  xy <- round(xy, -1)
  # make sf_polygon of each
  p <- st_sfc(st_polygon(list(xy)), crs = crs1)
  # 
  b <- x$bbox; colnames(b)[1:2] <- c("x","y")
  z <- data.frame(start=b[1,], end=b[2,])
  st_geometry(z) <- p
  z$series <- s
  z$type <- x$type
  z
})
# collect
topes <- bind_rows(topel)

## Done.
#
fout <- paste0("processed/panta_chull_polygons_and_prior_polygons_", format(Sys.Date(),"%F"), ".rds")

saveRDS(topes, fout)
