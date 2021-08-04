# Compile all gps "panta" information,
# and clean up to coherent serieses.
#
# v. 10.10.2019
#
# 3.8.2021: git-friendly
#
#################################################################
library(dplyr)
library(raster)
# Helpers:
posix2year <- function(x) {
  y <- as.numeric( format(x, "%Y") )
  i1 <- as.numeric( as.POSIXct(paste0(y, "/01/01 00:00:00"), tz = "UTC") )
  i2 <- as.numeric( as.POSIXct(paste0(y, "/12/31 23:59:59"), tz = "UTC") )
  s <- as.numeric(x)
  pmax(y, y + (s-i1)/(i2-i1))
}
#
#################################################################
# Gather up the tables:
#
# From Juha first week
a <- read.csv2("panta_raw/allIntensive.csv", sep  = ",", stringsAsFactors = FALSE)
panta1 <- a %>% mutate(x = as.numeric(x),
                       y = as.numeric(y),
                       ts = paste(year, mon, day, hour, min, sep = "-"),
                       timestamp = as.POSIXct( strptime(ts, format = "%Y-%m-%d-%H-%M", tz = "UTC") )
              ) %>% dplyr::select(name, x, y, timestamp)
# Tuike, just one file.
a <- read.csv2(paste0("panta_raw/panta_tuike.csv"), sep  = ",", stringsAsFactors = FALSE)
panta2 <- a %>% mutate(x = as.numeric(x),
                       y = as.numeric(y),
                       ts = paste(year, mon, day, hour, min, sep = "-"),
                       timestamp = as.POSIXct( strptime(ts, format = "%Y-%m-%d-%H-%M", tz = "UTC") )
) %>% dplyr::select(name, x, y, timestamp)
#
# Load newest 2018 panta data  "GPS_pannat_2018", cleaned up long series.
#
library(readxl)
fils <- dir("panta_raw/GPS_pannat_2018/", ".xlsx", full.names = TRUE)
panta3 <- NULL
for(fil in fils) {
  a <- read_xlsx(fil)
  names(a) <- tolower(names(a))
  # some tables missing "id", which is these files mean the name of the wolf.
  if(!"id"%in%names(a)) names(a)[2] <- "id"  # assume its second
  a$name <- a$id; a$id <- NULL
  # process
  tab <- a %>% 
    # drop missing location values
    filter( !(is.na(x) | is.na(y)) ) %>%
    mutate(
    # add POSIX timestamp
      ts = paste(date, gsub(",",":", time)),
      timestamp = as.POSIXct(ts, format =  "%F %X", tz = "UTC")) %>% dplyr::select(name, x, y, timestamp)
  # Convert from latlon to epsg:3067
  xy <- tab[,c("x","y")]
  coordinates(xy) <- c('x','y')
  proj4string(xy) <- CRS("+proj=longlat")
  tab[,c("x","y")] <- coordinates( spTransform(xy, CRS=CRS("+init=epsg:3067")) )
  # done.
  panta3 <- rbind(panta3, tab)
}
#
# Load as much as we can of the "reviireillä elävät alfasudet"-set
#
dirs <- c(dir("panta_raw/GPS_reviireillä_elävät_alfasudet/", full.names = TRUE),
          dir("panta_raw/GPS_reviireillä_elävät_extrat/", full.names = TRUE))
panta4 <- NULL
for(di in dirs) {
  files <- dir(di, "^[^~].*[.].", full.names = TRUE)
  filess <- dir(di, "^[^~].*[.].", full.names = !TRUE)
  iinfo <- grep("info", files)
  if(length(iinfo) == 0) {message(paste0("info missing for ", dire)); next}
  info <- as.matrix(read.table(files[iinfo]))
  idat <- (1:length(files))[-iinfo]
  nf <- nrow(info)
  out <- NULL
  for(i in 1:nf) {
    fn <- files[idat[i]]
    sheetn <- if(ncol(info) == 7) info[i,7] else 1
    a0 <- read_excel(fn, skip = info[i, 6], col_names = FALSE, sheet = sheetn)
    if(all(is.na(a0[,1]))) a0 <- a0[,-1]
    a <- a0[ , info[i, c(1:4)]]
    names(a) <- c("date", "time", "x", "y")
    # add name
    n <- gsub(".xls|.xlsx", "", filess[idat[i]])
    a$name <- n
    a <- a %>% mutate(x = as.numeric(x),
                      y = as.numeric(y),
                      cr = (x > 180 | x < -180 | y > 90 | y < -90) | (is.na(x)|is.na(y))) %>% 
      filter(!cr)
    # add POSIX timestamp
    date <- if(is(a$date, "POSIXct")) a$date else as.POSIXct(a$date)
    time <- if(is(a$time, "POSIXct")) a$time else as.POSIXct(a$date)
    a$timestamp <- as.POSIXct(paste(format(date, "%F"), format(time, "%X")), "%F %X")
    xy <- a[,c("x","y")]
    coordinates(xy) <- c('x','y')
    proj4string(xy) <- CRS("+proj=longlat")
    a[,c("x","y")] <- coordinates( spTransform(xy, CRS=CRS("+init=epsg:3067")) )
    # id ~ obs number
    tab <- a %>% dplyr::select(name, x, y, timestamp)
    out <- rbind(out, tab)
  }
  panta4 <- rbind(panta4, out )
}
######
# Gather.

panta <- rbind(panta1,
               panta2,
               panta3,
               panta4)
#################################################################
#
# Tweaks.
#
# # Drop locations outside finland?
# finl <- readRDS("finland-5x5-rastermask-and-polygon.rds")
# smask <- finl$mask #crop( finl$mask, finl$extsouth ) # Exclude Lappi
# okpoints <- !is.na(  smask[ cellFromXY(smask, panta[, c("x","y")])  ]   )
# panta$inside_finland <- okpoints
#
# Time in years, numeric
panta$timestamp_y <- posix2year(panta$timestamp) # for filtering, easier
#
# Change names to unique Wolf names. This is needed
# To remove overlaps between the sets.
unames0 <- unique(panta$name)
# Drop YYYY as this is in dates; drop "uros" etc.; drop "-_ " 
unames <- tolower( gsub("[0-9].|\\_|-| ", "", gsub( "alfa|naaras|uros|winter|summer|fin", "", unames0  )) )
nmap <- cbind(unames0, unames)
panta$name <- nmap[ match(panta$name, nmap[,1]) ,2]
#
#################################################################
# Then we process each series.
# 
pantal <- split(panta, panta$name)
# Average locations observed in under 1 minute 
time_th_min <- 1.01 # minutes
pantal_cleaned <- lapply(pantal, function(p) { 
  p <- p[order(p$timestamp_y),] 
  z <- p$timestamp_y * 365 * 24 * 60 # minuutteina 
  delta <- diff(z) 
  less <- delta <= time_th_min
  if(!any(less)) return(p) # early exit
  #
  lessi <- which(less)
  # clusters, in case of many
  cl <- rep(1,length(lessi)) 
  k <- 1 
  for(i in seq_along(lessi)[-1] ) cl[i] <- if(lessi[i]-lessi[i-1] < 2) k else (k <- k +1) 
  cll <- lapply(split(lessi,cl), function(v) c(v, max(v)+1)) # add the missing last index of group
  lessi <- unlist(cll)
  # average
  p_ok <- p[setdiff(1:nrow(p), lessi), ]
  pnew <- do.call(rbind, lapply(cll, function(i)  {
    pa <- p[i,]
    pn <- pa[1,]
    pn$x <- mean(pa$x)
    pn$y <- mean(pn$y)
    pn
  }   ))
  pout <- rbind(p_ok, pnew)
  pout <- pout[order(pout$timestamp_y), ]
  pout
})
panta2 <- do.call(rbind, pantal_cleaned)

############
# Store.
saveRDS(panta2, paste0("processed/pannat-", format(Sys.Date(), "%d%m%y"), ".rds"))
