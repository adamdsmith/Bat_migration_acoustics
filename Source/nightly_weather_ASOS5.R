# Loading required packages
toLoad <- c("plyr", "data.table", "lubridate", "circular", "reshape")
lapply(toLoad, require, character.only = TRUE)

# Inputting years and stations of interest
# This can/should change depending on the analysis of interest
years <- c("2010", "2011")

# Loading sunrise/sunset lookup table
lookup <- convert.magic(read.csv("./Data/set_rise_lookup.csv", header = T),
					   types = c("date", "character", "character"))
lookup <- within(lookup, {
  c_rise <- hm(paste(substr(c_rise, 1, 1), substr(c_rise, 2, 3)))
  c_set <- hm(paste(substr(c_set,1,2), substr(c_set, 3, 4)))
  dtset <- date + c_set
}) 

join.csvs <- function(file.names, ...) {
  ldply(file.names, function(fn) data.frame(read.csv(fn, header = T)))
}  

# Loading in data missing from ASOS5 and pulled from NCDC QCLCD
# Roughly hourly observations rather than every 5 minutes
missingASOS5 <- vector(mode = "list", length = length(years))
names(missingASOS5) <- years
for (i in 1:length(years)) {
  missingASOS5[[i]] <- convert.magic(read.csv("./Data/missingASOS5.csv", header = T), 
								    types = c("factor", rep("character", 5), "numeric", "factor", "numeric"))
  missingASOS5[[i]] <- subset(missingASOS5[[i]], year == years[i])
}

ASOS5 <- vector(mode = "list", length = length(years))
names(ASOS5) <- years
nightSky <- vector(mode = "list", length = length(years))
names(nightSky) <- years

# Creating lookup dataframe for sky conditions
skycode <- data.frame(code = c("CLR", "FEW", "SCT", "BKN", "OVC"), cover = c(0, 3/16, 7/16, 3/4, 1))

for (i in 1:length(years)) {
    year <- years[i]
    # Generating list of files for current year
    yearlist <- grep(year, list.files(pattern = ".csv$"), value = T)
    
    # Joining station files into year file using function created above
    ASOS5[[i]] <- join.csvs(yearlist)
    ASOS5[[i]] <- convert.magic(ASOS5[[i]], c("factor", rep("character", 5), "numeric", "factor", "numeric"))
    
    # Joining missing data
    ASOS5[[i]] <- rbind(ASOS5[[i]], missingASOS5[[i]])
  
    # First, create variable "night" that indicates date of sunset (i.e., date that recording night started)
    ASOS5[[i]] <- within(ASOS5[[i]], {
      month <- ifelse(nchar(month) == 1, paste0(0, month), month)
      day <- ifelse(nchar(day) == 1, paste0(0, day), day)
      hr <- ifelse(nchar(hr) == 1, paste0(0, hr), hr)
      min <- ifelse(nchar(min) == 1, paste0(0, min), min)
      date <- ymd(paste0(year, month, day))
      time <- hms(paste(hr, min, "00", sep = ":"))
      night <- as.POSIXct(ifelse(time < hours(12), date - ddays(1), date), origin = origin)
      keep <- ifelse(time <= lookup$c_rise[match(date,lookup$date)] | 
        time >= (lookup$c_set[match(date,lookup$date)] - hours(1)), "Y", "N")
      datetime <- ymd_hms(paste(year, month, day, hr, min, "00", sep = "-"))
    })
	
    ASOS5[[i]] <- within(arrange(subset(ASOS5[[i]], month(night) > 7 & keep == "Y"), station, datetime), {
      sinceset <- round((decimal_date(datetime) - decimal_date(lookup$dtset[match(night,lookup$date)])) * 8760, 3)
      hrset <- floor(sinceset)
      keep <- month <- day <- year <- hr <- min <- time <- date <- NULL
      skycov <- skycode$cover[match(sky,skycode$code)]
      vis <- ifelse(vis > 10, NA, vis) #checking for invalid values
    })
  
    # In some instances, duplicate records are present for a given date and time.
    # We have always found it to be the case that the first record is correct,
    # i.e., it corresponds to the data reported by the NCDC QCLCD.
    # The following line keeps the first record of duplicates
    ASOS5[[i]] <- subset(ASOS5[[i]], !duplicated(interaction(station, datetime)))
  

    # Calculating average sky cover, visibility, and cloud ceiling for each night
    nightSky[[i]] <- as.data.frame(cast(melt(ASOS5[[i]], id = "night", 
                                            measure = c("skycov", "vis", "ceiling"), na.rm = T),
                                       night ~ variable, mean))
    
}

save(ASOS5, file = "./Data/ASOS5.rda")
save(nightSky, file = "./Data/nightSky.rda")
