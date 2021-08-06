# Split and filter the gps serieses to contiguous serieses
library(dplyr)
library(tidyr)

# Use latest by date
fin <- dir("processed", "pannat", full.names = TRUE) %>% last()
panta0 <- readRDS(fin)

bywolf0 <- split(panta0, panta0$name)

print("start with %i serieses\n")

#
# Split long series to pieces with min th months in between
panta_splitter <- function(p, th_m, atleast=0, maxjump = 1000*100){
  p$name2 <- p$name
  p$rubbish <- FALSE
  ### Some holes present. Split over th_m long breaks
  dt <- diff(p$timestamp_y) * 12 
  xy <- p[,c("x","y")]
  n <- nrow(xy)
  dx2 <- rowSums(  ( xy[-1,]  - xy[-n,] )^2 )
  dx <- sqrt(dx2) 
  over <- which(th_m < dt | maxjump < dx)
  drop <- NULL
  if(length(over)){
    n <- p$name[1]
    over <- c(0, over, nrow(p))
    for(i in seq_along(over)[-1]){
      ii <- over[i+-1:0]+c(1,0)
      p$name2[ii[1]:ii[2]] <- paste0(n, i-1)
      if(diff(ii) < atleast)
        p$rubbish[ii[1]:ii[2]] <- TRUE
    }
  }
  p
}
MTH <- .25  # months, cut if longer
JUMPTH <- 1000*100 # 100km max jump length
#
# Also filter out series with less than some amount of points, say 24 x (ave step 1h) = 24h
panta1 <- do.call(rbind, lapply(bywolf0, panta_splitter, th_m = MTH, maxjump = JUMPTH, atleast = 24) ) %>% 
  filter(!rubbish)

#
##############################################################
# Then filter faulty gps-readings that jump too far from main series.
#
# - remove noise points that appear 100s of kilometers outside the path
# 
# Do this by computing each points' probability to be like others
#
byseries <- split(panta1, panta1$name2)
filter_leaps <- function(p) {
  p <- p[order(p$timestamp), ]
  n <- nrow(p)
  dt <- p$timestamp_y[-1] - p$timestamp_y[-n]
  dt_h <- dt * 365 * 24
  xy <- p[,c("x","y")]
  cxy <- apply(xy, 2, median)
  sxy <- apply(xy, 2, sd)
  dx2 <- rowSums(  ( xy[-1,]  - xy[-n,] )^2 )
  dx_km <- sqrt(dx2) / 1000
  speed_kmh <- dx_km/dt_h
  # # attach max +-1 step speed to each point
  kmh <- pmin( c(speed_kmh, Inf) , c(Inf, speed_kmh))
  # probabilties: 
  l0 <- -log(1-.5)/20 # p = 0.5 at 20kmh i.e. 20=q.5
  p1 <- exp(-l0 * kmh)
  # distance from mass center
  v <- unname(t(t(xy)-cxy))
  dxc <- sqrt( rowSums(v^2) )
  ss <- sd( dxc[dxc < quantile(dxc, .9)] ) # jump length sd
  # inside ball safe, beyond that drop off
  R <- 5*ss
  outside <- dxc > R
  p2 <- p1
  p2[outside] <- p2[outside] * exp( -0.5 * (dxc[outside]-R)^2 / (2 * 20000)^2 )
  #
  u <- p2
  keep <- u > 0.5
  if(0){ #check plot
    plot(xy, asp = 1, pch = ".", type="l")
    symbols(cxy[1], cxy[2], circles=R, inches=F, add=T, fg=3, bg=rgb(0,1,0,.1))
    points(xy)
    points(cxy[1],cxy[2], col=2, pch=19)
    points(xy[keep,], col =  3)
    points(xy[!keep,], col=2)
  }
  #
  keep
}

# tests and checks here:
if(0){
  p <- byseries$unna
  p <- byseries$kara1
  p <- byseries$apollo1
  p <- byseries$pinna4
  p <- byseries$vellu3
  keep <- filter_leaps(p)
  xy <- p[,c("x","y")]
  plot(xy, asp = 1, pch = 19, xlim = range(xy[keep,1]), ylim = range(xy[keep, 2]), type="l")
  points(xy[keep,], col = 3, pch= 20)
  points(xy[!keep ,], col = 2, pch= 20)
  print(mean(keep))
}
###########
# filter the series
byseries_filtered1 <- lapply(byseries, function(p) p[filter_leaps(p), ])
panta2 <- do.call(rbind, byseries_filtered1)
###
# Then do the timecut second time, i.e. iterate twice.
panta2$name0 <- panta2$name
panta2$name <- panta2$name2
byseries_filtered2 <- split(panta2, panta2$name2)
panta3 <- do.call(rbind, lapply(byseries_filtered2, 
                                panta_splitter, th_m = MTH, maxjump = JUMPTH, atleast = 24) ) %>% 
  filter(!rubbish)
byseries3 <- split(panta3, panta3$name2)
byseries_filtered <- lapply(byseries3, function(p) p[filter_leaps(p), ])
panta <- do.call(rbind, byseries_filtered)
panta <- panta %>% 
  mutate(series = name2) %>% 
  dplyr::select(series, name0, x, y, timestamp, timestamp_y) %>% 
  as_tibble()
#
#


#######################################################################
# Then add some derivatives

panta <- panta %>% 
  group_by(series) %>% 
  arrange(timestamp) %>% 
  mutate(
    timestamp_utc = timestamp,
    timestamp = format(timestamp), # to avoid UTC issues
    jumpsize = c(NA, sqrt(diff(x)^2 + diff(y)^2)),
    jumptime = c(NA, diff(timestamp_utc)),
    speed = jumpsize/(jumptime * 365 * 24) )



# store
fout <- gsub(".rds", "-filtered.rds", fin)
saveRDS(panta, fout)





