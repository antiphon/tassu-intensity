# Resample Corine 
#
# * reclass first to 8 classes
# * then layerize 
# * then compute proportions in 1x1km squares 
# * keep proportions of the 8 classes in a stack.

# Might want to set .Renviron TMPDIR to to somewhere with ~100GB space for raster temp-files.

library(raster)
library(parallel)

testing <- TRUE

################################################################
#
# Corine with ad-hoc chosen Values.
# pick and mix levels, lump together not-so-interesting ones.

###############################################################
# Remap the land classes to more practical:
maa <- read.csv("corine/CorineMaanpeite2012Luokat.csv", sep =";")
maa <- maa[!is.na(maa$Value),]
# New levels
frak <- maa$Level1 == 1
fasu <- maa$Level2 == 11
ftie <- maa$Level3 == 122
fmt  <- maa$Level1 == 2
fmet <- maa$Level1 == 3
fsum <- maa$Level2 == 31
mapl <- list("Asuinalueet" = which(fasu),
             "Tiet" = which(ftie),
             "Muu rakennettu" = which(frak & !fasu & !ftie),
             "Maatalous" = which(fmt),
             "Sulk. metsät" = which(fsum),
             "Harvat metsät yms." = which(fmet & !fsum),
             "Suot ja kosteikot" = which(maa$Level1 == 4),
             "Vesialueet" = which(maa$Level1 == 5)
)
map <- do.call(rbind, lapply(names(mapl), function(i) data.frame(class=i, idx=mapl[[i]], stringsAsFactors = FALSE)))
map$Val <- maa$Value[map$idx]
map <- map[order(map$Val),]
# check we have all in nice order and we don't need a match-call:
all(map$idx == map$Val & map$idx == 1:nrow(map))
# Integer coding for new raster
map$class_f <- factor(map$class, levels = names(mapl))
map$class_i <- as.integer(map$class_f)
# This map is important, store:
saveRDS(map, "corine/corine_custom_maaluokka_map-v1.rds")
#
############################################################
# Then reclassify the full Corine'18 raster.
#
# Get if from Syke:
# https://www.syke.fi/fi-FI/Avoin_tieto/Paikkatietoaineistot/Ladattavat_paikkatietoaineistot, "CORINE maanpeite 2018, 20m"
# 
# Assuming here the loaded geotiff is in
corine_big_path <- "~/temp/corine/"
corine_full_file <- paste0(corine_big_path, "Clc2018_FI20m.tif") #
if(!file.exists(corine_full_file)) stop("full corine file not found.")

# Load.
corine0 <- raster(corine_full_file)
map <- readRDS("corine/corine_custom_maaluokka_map-v1.rds")


# for testing
if(testing) {
  subw <- extent(corine0)*0.05
  corine0 <- crop(corine0, subw)
  message("Corine Test run.\n")
}

# remap classes 
# in full rez, takes a while.
corine8 <- calc(corine0, function(v) map$class_i[v])


# Intermediate file? Write takes a while.
#writeRaster(corine8, paste0(corine_big_path, "clc2018_fi20m_reclassed1.tif") )


# Run the above only once.
################################################
#
## The resampling of the 8 classes starts here.

#corine8 <- raster( paste0(corine_big_path, "clc2018_fi20m_reclassed1.tif") )
#
# for checking.
subw <- extent(corine8) * 0.05
cor8sub <- crop(corine8, subw)

###########################################################
# Resample

# Explozde
message("layerising")

bsl <- layerize(cor8sub, classes = 1:8) # test first

# go!
beginCluster(n = min(detectCores(), 8) ) 
corine8_s <- layerize(corine8, classes = 1:8)
# Attach informative names
land_classes <- levels( map$class_f )
#names(bsl) <- land_classes
names(corine8_s) <- land_classes
#
# Then we aggregate to get fraction of each class type:
message("aggregating")
# target
W <- raster(extent(corine8), resolution = c(1000, 1000), crs = crs(corine8))

corine8_agg <- aggregate(corine8_s, 
                     fact = (dim(corine8_s)/dim(W))[-3], # careful not to aggregate over layers
                     fun = mean, na.rm=TRUE) 

endCluster()
# Done.


#store
message("writing")
writeRaster(corine8_agg, 
            file = paste0("corine/clc2018_fi_reclassed1_1km", ifelse(testing, "_testing", "")), 
            format = "raster")  # use native to keep layer names. Remember to copy both files!
message("done.")




# 
