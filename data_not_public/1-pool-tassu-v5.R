# Pool all Tassu-sources. 
# 
# Also homogenise dates etc. 
#
# v.5. 19.11.2020: 
#  - Pool all so far
#  - Mark each row by how many times it is present (identical) 
#  - do nothing else
#
# v.5.1: 11.12.2020:
#  - completely new file from Samu , Henri Näpärä
#
# v5.2: 14.1.2020:
#  - add subindex for havainto_ids
#
# v5.3: 3.8.2021:
#  - make git-friendly
#
####################################################################

library(dplyr)
library(tidyr)

####################################################################
# Helper:
# Turn a date into a running number in years. Careful with leap years etc.

posix2year <- function(x # as POSIX.clt
                       ) {
  y <- as.numeric( format(x, "%Y") )
  i1 <- as.numeric( as.POSIXct(paste0(y, "/01/01 00:00:00"), tz = "UTC") )
  i2 <- as.numeric( as.POSIXct(paste0(y, "/12/31 23:59:59"), tz = "UTC") )
  s <- as.numeric(x)
  pmax(y, y + (s-i1)/(i2-i1))
}
#
#################################################################
# Go
#################################################################
# Read raw files in.
#fils <- dir("tassu_raw/", "csv$", full.names = TRUE)
fils <- "tassu_raw/data-1607608115908.csv"
traw <- NULL
for(infi in fils) {
  cat("File:", infi, "\n")
  traw <- rbind( traw, read.csv2(infi, sep=";", 
                                 stringsAsFactors = FALSE , 
                                 #encoding = "latin1", 
                                 na.strings = "NULL") %>% mutate(file = infi))
}
#
############################################
# Make sense of the timestamps
ts1 <- traw$pvm
ts2 <- traw$aika
# Some stamps missing leading 2.
ts1 <- gsub("000", "200", ts1)
ts1 <- gsub("001", "201", ts1)
ts2 <- gsub("000", "200", ts2)
ts2 <- gsub("001", "201", ts2) # remove this if time goes back to 2001.
#
# check if just date available (i.e. missing 'aika'), and replace with midday (<-!!!! Choice! )
i2 <- nchar(ts2) == 0 | is.na(ts2)
tts2_2 <- sapply(which(i2), function(i){ 
  out <- ts1[i]
  paste(out, "12:00:00+2")
})
tts2 <- ts2
if(any(i2)) tts2[i2] <- tts2_2
# proper label with clock
tts2 <- as.POSIXct(tts2, format = "%F %T")
# Check
tok <- !is.na(tts2)
cat("Dropping strange date rows, n =", sum(!tok), "(", round(mean(!tok)*100, 4), "%)\n")
#
traw_okdate <- traw[tok,]
traw_okdate$timestamp <- tts2[tok]
# Also time in running years
traw_okdate$timestamp_y <- posix2year(traw_okdate$timestamp)

# Timestamp cleaning done.

tassu <- traw_okdate %>% 
# Location
  mutate(st_x = as.numeric(st_x),
         st_y = as.numeric(st_y),
# some minor tweaks 
         pihapiiri = pihapiiri == "True",
         ryhmassa = ryhmassa == "True",
         tarkastettu = tarkastettu == "True",
         laji = laji_nimi,
         laji_nimi = NULL)

################################################
## then count duplication space-time-coordinatewise
tassu <- tassu %>% 
  mutate(sid = 1:n()) %>% 
  group_by(laji, havaintotyyppi, havaitsija_id, timestamp, st_x, st_y) %>%
  mutate(    similar_rows_id = sid[1],
         similar_rows_subidx = 1:n(),
         similar_rows_count  = n()) %>% ungroup() %>% mutate(sid = NULL)

# Add info on the possible number of animals
tassu <- tassu %>% 
  group_by(havainto_id) %>%
  mutate(havainto_id_subidx = 1:n(),
         havainto_id_count = n()) %>% ungroup()



# Tassu tables compiled.
message("Tassu basic table ready.")
############################################
# Done. Save.
saveRDS(tassu, 
        paste0("processed/tassu_", format(Sys.Date(), "%d%m%y"), ".rds"))
#####################################
#



