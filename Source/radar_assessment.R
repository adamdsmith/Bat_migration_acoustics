toLoad <- c("sp", "fields", "sp", "lubridate", "plyr", "reshape")
instant_pkgs(toLoad); rm(toLoad)

# Place all reflectivity files of interest (i.e., those in the just created text file), 
# un-tarred, in a directory.  Then, batch export them to ASCII grid files using the NOAA 
# Weather & Climate Toolkit (WCT) 

# gzfiles contains information about *.gz files from which
# ascii grid files (*.asc) were created via export from the WCT.  
# This object contains, for each night:
# - time of civil sunset (c_set)
# - corresponding UTC time (i.e., radar time), plus 32 minutes to put it one hour post-sunset
# - the name of the *.gz file closest to one hour post-sunset
# - the date (night)
# - the number of minutes difference between the actual radar data and one hour post-sunset (diff)

# gzfiles should be placed in the directory that contains the actual files
# This directory, in turn, should be set as the working directory
setwd("D:/NEXRAD/Sunset_radar/Analyzed")

# Now, load gzfiles
load("./Data/NEXRAD_scan_metadata.rda")

# Loading 150 km mask
# For whatever reason, NOAA's WCT would not apply the range limits (0 - 150 km)
# that we attempted to specifiy in the config xml file, so we created our own ascii grid mask.
# We multiply the mask's cell values (1 within 150 km, NA outside of 150 km) times
# the corresponding reflectivity value in the NEXRAD image to mask out reflectivity
# beyond 150 km from the radar site; a mask will have to be created separately for each
# distinct NEXRAD site
NEXRAD_mask <- read.asciigrid("./Data/kbox_150_mask.asc", as.image=T)
NEXRAD_mask$z[NEXRAD_mask$z == 0] <- NA

# Loading radar volume coverage patterns (VCP) data
# VCP indicate whether the radar is operating in "clear air" mode (VCPs 31 and 32)
# or precipitation mode (VCPS 11, 21, 12, or 121)
# This information must be gathered manually be inspecting the base reflectivity images
load("./Data/NEXRAD_precip_mode.rda")
gzfiles <- join(gzfiles, precip_mode, by = "night")

# This function iterates through the list of ASCII grid files, applies the 150 km mask, 
# filters unexpectedly high dBz values when precipitation is not present,
# and calculates varies measures of reflectivity
radarAnalyze <- function(filelist, plot = FALSE) {
  # Input filelist is pathname to a text file containing the ascii grid files to be analyzed
  # Loading required packages
  toLoad <- c("sp", "fields", "lubridate", "plyr", "reshape")
  lapply(toLoad, library, character.only = TRUE)
  
  filelist <- read.table(filelist, header=F)
  n <- nrow(filelist)
  output <- data.frame(night = rep(ymd("20000101"), n), 
                       dbZ = rep(0, n), Z = rep(0, n), 
                       Eta = rep(0, n), logEta = rep(0, n))

  # Iterating through the files, calculating what we need, 
  # putting it into a table and moving on...
  for (i in 1:n) {
    fn <- paste(getwd(), filelist[i, 1], sep = "/")
    grd <- read.asciigrid(fn, as.image=T)
    
    # Gathering night of radar image and VCP status for later use
    radarnight <- gzfiles$night[match(substr(as.character(filelist[i, 1]), 1, 19),
                                      substr(gzfiles$file, 1, 19))]
    VCP <- ifelse(gzfiles$VCP[match(substr(as.character(filelist[i, 1]), 1, 19),
                                      substr(gzfiles$file, 1, 19))] == -1, "Clear air", "Precip")
    precipYN <- ifelse(gzfiles$precip150[match(substr(as.character(filelist[i, 1]), 1, 19),
                                    substr(gzfiles$file, 1, 19))] == -1, "no", "yes")
    
    # Creating temporary data frame summarizing dbZ 
    dbZ_vect <- as.vector(grd$z) * as.vector(NEXRAD_mask$z) # Applying 150km mask
    dbZ_df <- data.frame(table(dbZ_vect))
       # A bit of data frame cleanup
      names(dbZ_df) <- c("dbZ", "freq")
      dbZ_df$dbZ <- as.numeric(as.character(dbZ_df$dbZ))
    
    # Filtering values above 32.5 dbZ (maximum reflectivity for migrating birds;
    # Gauthreaux and Belser 1998, Gauthreaux et al. 2008) when radar operated in "clear air" mode 
    dbZ_df <- if(precipYN == "no") {
      subset(dbZ_df, dbZ <= 32.5)
      } else {
        dbZ_df
      }
    
    # Creating new variables to be summarized
    dbZ_df <- within(dbZ_df, {
      # linearizing dbZ to Z
      Z <- 10^(dbZ/10)
      # Convert to radar reflectivity (Eta) of bioscatter (Chilson et al. 2012, Ecosphere 3:72)
      logEta <- dbZ + 10*log10(10^3 * pi^5 * 0.93 / 10.7^4) # Eq. 19 in Chilson et al. 2012
      Eta <- 10^(logEta/10) # linear units
    })
    
    # Plotting, if requested; this will add some processing time
    if (plot) {
      grd$z <- matrix(dbZ_vect, nrow = 1000)
      if (i == 1) {
        x11()
        image.plot(grd)
        text(radarnight, x = min(grd$x), y = max(grd$y) - (max(grd$y) - min(grd$y))/20, pos = 4, cex = 2)
        text(VCP, x = min(grd$x), y = max(grd$y) - (max(grd$y) - min(grd$y))/10, pos = 4, cex = 2)
        Sys.sleep(0.5)
      } else {
        image.plot(grd)
        text(radarnight, x = min(grd$x), y = max(grd$y) - (max(grd$y) - min(grd$y))/20, pos = 4, cex = 2)
        text(VCP, x = min(grd$x), y = max(grd$y) - (max(grd$y) - min(grd$y))/10, pos = 4, cex = 2)
        Sys.sleep(0.5)
      }
    }
    
    # Adding sum of these variables, weighted by their frequency, into the output data frame
    weights <- dbZ_df$freq
    output[i, "dbZ"] <- sum(dbZ_df$dbZ * weights)
    output[i, "Z"] <- sum(dbZ_df$Z * weights)
    output[i, "Eta"] <- sum(dbZ_df$Eta * weights)
    output[i, "logEta"] <- sum(dbZ_df$logEta * weights)
    # Adding night of radar image
    output[i, "night"] <- radarnight
    cat("\nProcessing:", as.character(radarnight))
  }

  return(output)
}

reflectivity <- radarAnalyze(filelist)
  
# Creating the final combined data set
reflectivity <- join(reflectivity, precip_mode, by = "night")
years <- c("2009", "2010", "2011", "2012")
reflect_list <- vector(mode="list", length=length(years))
names(reflect_list) <- years
for (i in 1:length(years)) {
  reflect_list[[i]] = subset(reflectivity, year(night) == years[i])
}
reflectivity <- reflect_list
save(reflectivity, file = "./Data/NEXRAD_reflectivity.rda")

# Use this function to plot the filtered radar image (i.e., within 150 km) for 
# a given night
# Requires "sp" and "fields" packages to be loaded

plotRadar <- function(stringymd) {

  if (nchar(stringymd) != 10 | substr(stringymd, 5, 5) != "-") {
    stop("Date (stringymd) must in the format YYYY-MM-DD")
  }
  
  # Get list of files containing used radar scans (filelist)
  load("./Data/NEXRAD_file_list.rda")

  # Get metadata associated with used radar scan (e.g., minutes from sunset; gzfiles)
  load("./Data/NEXRAD_scan_metadata.rda")

  # Radar precip mode status (precip_mode)
  load("./Data/NEXRAD_precip_mode.rda")
  
  if (!exists("NEXRAD_mask")) {
    NEXRAD_mask = read.asciigrid("./Data/kbox_150_mask.asc", as.image=T)
    NEXRAD_mask$z[NEXRAD_mask$z == 0] <- NA
  }

  gzfiles <- join(gzfiles, precip_mode, by = "night")
  # Find the file containing the correct radar scan
  findFile <- gzfiles$file[which(as.character(gzfiles$night) == stringymd)]
  radarScan <- which(filelist[,1] == paste0(findFile, ".asc"))
  
  # Set directory of used radar scans
  scan_loc <- "E:/NEXRAD/Sunset_radar/Analyzed"
  fn <- paste(scan_loc, filelist[radarScan, 1], sep = "/")
  grd <- read.asciigrid(fn, as.image=T)
  dbZ_vect <- as.vector(grd$z) * as.vector(NEXRAD_mask$z)
  grd$z <- matrix(dbZ_vect, nrow = 1000)
  image.plot(grd, asp=1)
  radarnight <- gzfiles$night[match(substr(as.character(filelist[radarScan, 1]), 1, 19),
                                    substr(gzfiles$file, 1, 19))]
  VCP <- ifelse(gzfiles$VCP[match(substr(as.character(filelist[radarScan, 1]), 1, 19),
                                    substr(gzfiles$file, 1, 19))] == -1, "clear air", "precip")
  text(radarnight, x = min(grd$x), y = max(grd$y), pos = 4, cex = 2)
  text(paste("mode:", VCP), x = max(grd$x), y = max(grd$y), pos = 2, cex = 2)
}
