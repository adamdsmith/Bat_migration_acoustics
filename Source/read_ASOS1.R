# Loading required packages
toLoad <- c("plyr", "data.table", "lubridate", "circular", "ggplot2")
instant_pkgs(toLoad); rm(toLoad)

# Inputting years and stations of interest
# This can/should change depending on the analysis of interest
years <- c("2009", "2010", "2011", "2012")
stations <- c("KGON", "KMTP", "KPVD", "KUUU", "KWST")

# Creating empty data frame to store sample size summaries of stations
sample_summary <- data.frame()

# Creating functions used to read and combine PAGE 1 and PAGE 2 files
read.tables_p1 <- function(file.names, ...) {
  ldply(file.names, function(fn) data.frame(read.fwf(fn, widths=page1widths, col.names=page1names, stringsAsFactors=F)))
}  
read.tables_p2 <- function(file.names, ...) {
  ldply(file.names, function(fn) data.frame(read.fwf(fn, widths=page2widths, col.names=page2names, stringsAsFactors=F)))
}  
# Spacings for fixed width files.  Negative values drop columns, positive values keep columns
page1widths <- c(-5, 4, -4, 8, 4, -7, 5, -1, 1, -29, 3, -4, 2, -18, 2, -15)
page2widths <- c(-13, 8, 4, -19, 4, -21, 6, -2, 6, -2, 6, -3, 3, -3, 2, -10)
# Names of retained variables (columns)
page1names <- c("station", "date", "time", "extcoef", "sensor", "wdir", "wsp", "rvr")
page2names <- c("date", "time", "precip", "press1", "press2", "press3", "temp", "dewpt")

# Loading civil sunrise/sunset lookup table and converting date to useful format
# The sunrise/sunset lookup is a location-specific file created by the user
# The lookup used in this code is available for comparison: http://goo.gl/vvMPs
lookup <- read.csv("./Data/set_rise_lookup.csv", header=T)
lookup$date <- mdy(lookup$date)

# Create progress bar
total <- length(years) * length(stations)
pb <- winProgressBar("Processing ASOS data", "Progress updates after completion of first iteration",
                    min=0, max=total, 300)

for (year in years) {
  i <- which(years == year) # for use in progress bar

  # Loading in data missing from ASOS and pulled from NCDC QCLCD, if necessary
  # Roughly hourly observations rather than every minute
  missingASOS <- paste0("./Data/missingASOS", year, ".csv")
  checkFile <- file.exists(missingASOS)
  if (checkFile) {missingData <- read.csv(missingASOS, header = T)}
  
  for (station in stations) {
    j <- which(stations == station) # for use in progress bar

    ASOSdir <- paste("D:/ASOS 1-minute data", year, station, sep="/")
    
    # Generating list of PAGE 1 and PAGE 2 files for current year and station
    page1list <- paste(ASOSdir, grep("64050[K]", list.files(path = ASOSdir, pattern=".dat$"), value=T), sep = "/")
    page2list <- paste(ASOSdir, grep("64060[K]", list.files(path = ASOSdir, pattern=".dat$"), value=T), sep = "/")

    # Reading monthly (individual) files into PAGE 1 and PAGE2 files using functions created above
    page1 <- read.tables_p1(page1list)
    page2 <- read.tables_p2(page2list)
    
    # Joining two pages, and removing irrelevant observations (between civil rise and an hour pre-civil set)
    pages12 <- join(page1, page2, by=c("date", "time"))
    
    # Adding data for missing night(s); not ideal because data from multiple stations get added to first 
    # created file in each year, but not going to toil with making separate files, etc...
    if (checkFile == T & j == 1) {pages12 <- rbind(pages12, missingData)}
    
    # First, create variable "night" that indicates date of sunset (i.e., date that recording night started)
    pages12 <- within(pages12, {
      date <- ymd(date)
      night <- as.POSIXct(ifelse(time < 1000, date - ddays(1), date), origin = origin)
      keep <- ifelse(time <= lookup$c_rise[match(date,lookup$date)] | 
        time >= lookup$c_set[match(date,lookup$date)] - 700, "Y", "N")
    })
    pages12 <- subset(pages12, keep == "Y")  
    pages12$keep <- NULL
    
    # In some instances, duplicate records are present for a given date and time.
    # We have always found it to be the case that the first record is correct,
    # i.e., it corresponds to the data reported by QCLCD
    # The following line keeps the first record of duplicates
    pages12 <- subset(pages12, !duplicated(interaction(date,time)))
    
    # Coaxing variables into their proper data type, and replace missing values with NAs
    # Note, this will generate warnings that "NA" values were coerced.  These can be ignored,
    # although the user should check that other warnings are not hidden amongst them.
    pages12 <- within(pages12, {
      extcoef <- as.numeric(extcoef)
      sensor[sensor != "D" & sensor != "N"] <- NA
      wdir <- as.integer(wdir); wsp <- as.integer(wsp); rvr <- as.integer(rvr); temp <- as.integer(temp)
      dewpt <- as.integer(dewpt); precip <- as.numeric(precip); press1 <- as.numeric(press1)
      press2 <- as.numeric(press2); press3 <- as.numeric(press3)
      #The following section is quality control.  The settings are application-specific, based on
      # an thorough review of each input ASOS file for each station & year.
      precip[precip > 0.5] <- NA 
      temp[temp < 20] <- NA 
      temp[temp > 100] <- NA
      press1[press1 > 30.6] <- NA # Corresponds to a pressure of about 1035 mb
      press3[press2 > 30.6] <- NA
      press3[press3 > 30.6] <- NA
      dewpt <- ifelse(dewpt > temp, NA, dewpt) 
    })
    
    # VARIABLE CREATION OR MODIFICATION

    # First, remove remnant of unsampled night at beginning of records (in this case, night of July 31)
    pages12 <- subset(pages12, month(night) > 7)
    
    # Now, some tidying and preparation
    # Modifying sunset lookup table to facilitate date match
    lookup$dtset <- lookup$date + hm(formatC(lookup$c_set/100, digits=2, format="f"))
    
    # set reference direction for wind profit calculation below
    ref_dir <- 135 
    
    # Relative humidity function
    # August-Roche-Magnus formula (Lawrence 2005, Bull. Amer. Meteor. Soc. 225-233)
    # Values of A and B from Alduchov and Eskridge (1996, J App Meteor 35:601-609)
    relhum <- function(tempC, dewptC) {
      A <- 17.625; B <- 234.04
      e <- exp((A * dewptC)/(B + dewptC))
      es <- exp((A * tempC)/(B + tempC))
      RH <- 100 * (e / es)
      return(RH)
    }
    
    
    # Finally, create some new variables, etc.
    pages12 <- within(pages12, {
      
      # Calculate time since sunset...
      datetime <- date + hm(formatC(time/100, digits=2, format="f"))
      sinceset <- round((decimal_date(datetime) - decimal_date(lookup$dtset[match(night,lookup$date)])) * 8760, 3)
      hrset <- floor(sinceset)
      
      # Following NWS convention, we set wind direction to 0 for observations with wind speed = 0
      wdir <- ifelse(wsp == 0, 0, wdir)
      
      # converting wind speed from knots to m/s
      wsp <- round(wsp * 0.514444, 2)
      
      # calculating wind profit according to Erni et al. (2002, Ardea 90:155-166), desired heading to SE (135 degrees)
      wprof12 <- round(ifelse(wdir >= 180,
                           12 - sqrt(12^2 + wsp^2 - 2 * 12 * wsp * cos(rad(wdir - 180 - ref_dir))),
                           12 - sqrt(12^2 + wsp^2 - 2 * 12 * wsp * cos(rad(wdir + 180 - ref_dir)))), 2)
            
      # converting temperature and dewpoint to Celsius, calculating relative humidity
      tempC <- round((temp - 32) / 1.8, 1) 
      dewptC <- round((dewpt - 32) / 1.8, 1)
      rh <- round(relhum(tempC, dewptC), 1)

      #convert average pressure (from up to 3 sensors) into single pressure measurement "mb"
      mb <- round(rowMeans(pages12[,c("press1", "press2", "press3")], na.rm = T) * 33.863886667, 1)
      
      # Dropping some unused variables
      press1 <- press2 <- press3 <- rvr <- temp <- dewpt <- dewptC <- datetime <- NULL
    })

    # Creating summary data frame of samples/day; appending to summary table
    temp_sum <- ddply (pages12, .(date), summarise, page1 = sum(!is.na(wsp)), page2 = sum(!is.na(tempC)))
    temp_sum$station <- station
    sample_summary <- rbind(sample_summary, temp_sum)
    
    # Writing individual stations_year *.csv file
    export <- paste0("./Data/", station, year, ".csv")
    write.table(pages12, export, sep=",", row.names = F, qmethod = "d")

    # adding station data to single *.csv file for each year
    export <- paste0("./Data/ASOS", year, ".csv")
    if (j == 1) {
      write.table(pages12, export, sep=",", row.names=F, qmethod="d")
    } else {
      write.table(pages12, export, sep=",", col.names=F, row.names=F, qmethod="d", append=T)
    }
           
    Sys.sleep(0.1)
    k <- (i - 1) * length(stations) + j
    info <- sprintf("%d%% done", round(k/total*100))
    setWinProgressBar(pb, k, paste("Processing ASOS data for", station, "in", year), info)
  }
}
close(pb)

# Creating combined ASOS1 data set
years <- c("2009", "2010", "2011", "2012") # Change as necessary
ASOS <- vector(mode="list", length=length(years))
names(ASOS) <- years
for (i in 1:length(years)) {
  ASOS1[[i]] <- read.csv(paste0("./Data/ASOS", years[i], ".csv"), header=T)
  ASOS1[[i]] <- subset(ASOS1[[i]], sinceset < -6 | sinceset >= -1)
  ASOS1[[i]] <- within(ASOS1[[i]], {
    date <- ymd(date)
    night <- ymd(night)
  })
  ASOS1[[i]] <- arrange(subset(ASOS1[[i]], month(night) < 11), date, time, station)
}

# Manualy adding pre-sunset wind speed and wind direction (from NCDC QCLCD) for 29 Oct 2012
ASOS1[[4]][336046, ] <- ASOS1[[4]][336045, ]
ASOS1[[4]][336046,] <- data.frame("KGON", ymd("20121029"), 1616, NA, NA, 90, 16.09, NA, ymd("20121029"), NA, NA, NA, NA, -1, -0.947)
ASOS1[[4]][336047,] <- data.frame("KUUU", ymd("20121029"), 1616, NA, NA, 100, 16.54, NA, ymd("20121029"), NA, NA, NA, NA, -1, -0.947)

save(ASOS1, file = "./Data/ASOS1.rda")
