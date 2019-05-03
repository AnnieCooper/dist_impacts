# -----------------------------------------------------------------------------------
# BFAST_TEST_v0.R
# Z Holden, E Landguth
# September 2012
# v0 - Initial BFAST testing
# v1 - Reads in Time series pixel locations data, runs BFAST, gets breakpoint,
# recreates a new file with time series pixel locations of 0 and 1 (for breakpoint
# locations)
# v2 - Interpolate missing data values if any exist in time series data
# v4 - Store break points to file
# V5 - Reads in break point data file and creates time series data from that (mid way through script)
# v6 - Reads in raster stack, runs bfastmonitor.
# v7 - "" Cleaned up...added time stamps - NA checks.
# v8 - NA checks more...
# v9 - Move out of apply function calc and loop through each pixel.
# v10 - Parallel code.
# v11 - Fixed NA code - error in loess.
# v12 - Change to run bfast() - includes filtered MODIS imager (forest no forest)
# - includes filter downtrend at end - includes cummulative raster at end.
# -----------------------------------------------------------------------------------
# BFAST_TEST_cooper.R
# L Cooper
# Last edited: March 25, 2015.
# v13 - Changed the masking/cropping process. First crops to extent of 		region's ads surveys (bbox). Then masks out non-forested areas.
# v14 - Added another set of output rasters to give me percent change and difference in EVI.
# v15 - Changed Pixel.List to only include indices that are not masked.
# v16 - New code to allow for NAs at beginning and end of ts. 
#------------------------------------------------------------------------------------

# Packages Load
# ---------------
library(bfast)
library(wq)
library(raster)
library(date)
start.time <- Sys.time()

# Location files
# ---------------
datadir <- "/home/um/lcooper/data/lcooper/BFASTmodis/EVI/R1/" #Make sure to load and delete the files each time.
outdir <- "/home/um/lcooper/data/lcooper/BFASTmodis/output/R1/" #change each time

#change tmp directory.
rasterOptions(tmpdir="/home/um/lcooper/data/lcooper/BFASTmodis/tmp/")

# Read in rasters
# ---------------

# List files
fileList <- list.files(datadir, pattern="EVI.tif$", full.names=TRUE,recursive=FALSE,include.dirs=FALSE)

# Convert to brick. #Clouds, etc will already be masked out and files clipped.
#Convert to brick
for(i in 1:length(fileList)) 
{
  if(i==1) 
  {
    tempr <- raster(fileList[i])
    newbrick <- tempr
  }
  if(i > 1) 
  {
    r <- raster(fileList[i])
    newbrick <- stack(newbrick, r)
  }
}
modisraster <- brick(newbrick, values=T)
remove(newbrick)
#modisraster <- stack(fileList)
gc()

# Mask modisraster here
# ---------------------
maskfile <- "/home/um/lcooper/data/lcooper/BFASTmodis/maskFiles/R1_forestArea.tif" #Land Cover that is at least 15% forest (according to modis vcf product).
mpbrast <- raster(maskfile) 
NAvalue(mpbrast) <- 0 #Tell R what the NA value is for the mask file.
print("Resampling...")
mpbrast<-resample(mpbrast,modisraster,method="ngb") #Make the extent, resolution, etc. the same. 
print("Masking...")
gc()
m<-mask(modisraster[[1]],mpbrast) #Work-around to mask raster stack (mask doesn't work on stacks).
modisraster<-mask(modisraster,m)
gc()

# Print out statement
print(paste("Brick read in time:"))
print(Sys.time() - start.time)

# Grab dates. 
dates<-as.numeric(substr(fileList,77,80))+as.numeric(substr(fileList,81,83))*(1/1000)

# ------------------------

# Parallel Code
#
#This is used for the western US run.
# -------------

# Timing start
parallel.time <- Sys.time()

# Create a list of your MODIS pixels here. 
Pixel.List <- which(((getValues(mpbrast))==1),arr.ind=TRUE)
Pixel.List <- as.list(as.vector(Pixel.List))

# Now write function that operates on a single cell/vector from your MODIS brick
# every instance of ipix from your loop was replaced with ModisPixel.
bfast_funct <- function(ModisPixel,modisraster,dates)
{  
  # Remove the loop through pixels.
  pixdrill <- as.integer(modisraster[ModisPixel])
  
  # If there are NAs in the entire string - skip and return NAs
  if(sum(is.na(pixdrill)) == length(pixdrill))
  {
    breakpoint <- NA
  }
  # If less than 1/4 of the pixels are non-NA values, return NA.
  if (sum(!is.na(pixdrill)) <= (length(pixdrill)/4))
  {
    breakpoint <- NA
  }
  #If all of the above are false, continue.
  if(sum(is.na(pixdrill)) != length(pixdrill) & sum(!is.na(pixdrill))>(length(pixdrill)/4))    
  {
    tt <- 1:length(pixdrill)
    resid <- residuals(loess(pixdrill~tt))
    resid.q <- quantile(resid,prob=c(0.25,0.75))
    iqr <- diff(resid.q)
    limits <- resid.q + 1.5*iqr*c(-1,1)
    score <- abs(pmin((resid-limits[1])/iqr,0) + pmax((resid - limits[2])/iqr,0))
    pixdrill[score > 0.0] = NA
    
    # Remove the NAs with interpTS package
    pixdrill <- interpTs(pixdrill,gap=20)
    
    # If there are more than 20 NA times left in the series, return NA.
    if(sum(is.na(pixdrill)) > 20){     
      breakpoint <- NA
    }
    
    # If there are fewer than 20 NA values left (they will be at the end or beginning), continue with revised start/end dates.
    if(sum(is.na(pixdrill)) <=20){
      firstDate = which(!is.na(pixdrill))[1]
      lastDate = tail(which(!is.na(pixdrill)))[6]
      tempDate = dates[firstDate:lastDate]
      pixdrill = pixdrill[firstDate:lastDate]
      df <- data.frame(time=tempDate,test=pixdrill)
      if (floor(tempDate[1]==2000)){
        startInd <- (20-(sum(floor(tempDate)==2000))) + 4 #dates in 2002, subtract from 20, add 4
        endInd <- (23 - (sum(floor(tempDate)==floor(tail(tempDate)[6])))) + 1
      }
      else{
        startInd <- (23 - (sum(floor(tempDate)==floor(tempDate[1])))) + 1#dates in final year, subtract from 23, add 1
        endInd <- (23 - (sum(floor(tempDate)==floor(tail(tempDate)[6])))) + 1
      }
      df.ts=ts(df$test, start=c(floor(tempDate[1]),startInd), end=c(floor(tail(tempDate)[6]),endInd), frequency=23)
      
      # Run BFAST - defaults from manualS
      fit <- bfast(df.ts,h=0.15, season="harmonic", max.iter=1, breaks = 1)
      breakpoint <- dates[floor(fit$Time)] #just returns the year - could change to return the date too
    }
  }
  
  # Return
  return(breakpoint)
  
} # End of function

################################################# 

# Test to make sure this works on one pixel in your brick:
(res <- bfast_funct(Pixel.List[[1]],modisraster,dates))

# Now try parallel (careful of hanging processors)

library(snowfall)
sfStop() # just in case
sfInit(parallel=T,cpus=32) # Open number of processors...
sfLibrary(bfast) # pass libraries...
sfLibrary(wq)
sfLibrary(date)

# actual call: This loads up list of runs and fires them off to procs as they finish
cat(system.time(result <- sfClusterApplyLB(Pixel.List,bfast_funct,modisraster,dates)))

# Stop snowfall.
sfStop()

# Unlist results and write to the empty raster. 
one<-data.frame(x=which((getValues(mpbrast==1)),arr.ind=TRUE))
one$x2<-unlist(rbind(result))
two<-data.frame(y=seq(1,ncell(modisraster),1)) #List of all pixels in modisraster.
merged<-merge(one,two,by.x="x",by.y="y",all=TRUE) #Merge unmasked pixel values with entire list of possible pixels.

p.bfast.year <- raster(modisraster[[1]]) #Create random correctly-sized raster.
p.bfast.year <- setValues(raster(p.bfast.year),merged$x2) #Change all values of random raster to bfast values.

# Print out statement
print(paste("Finished Loop time: "))
print(Sys.time() - parallel.time)


# Split up into yearly rasters.
# ----------------------------

minVal <- as.integer(cellStats(p.bfast.year,min))
maxVal <- as.integer(cellStats(p.bfast.year,max))
for(i in minVal:maxVal)
{
  tempr <- floor(p.bfast.year)
  # Where values are i make 1, else make 0
  tempr[tempr == i] <- 1
  tempr[tempr != i & tempr != 1] <- 0
  
  # Write out raster
  writeRaster(tempr,file=paste(outdir,'bfast',as.character(i),'.tif',sep=""))
  
}

# Filter out up trends
# --------------------
# List files - TIME SERIES
fileList <- list.files(outdir, pattern=".tif", full.names=TRUE)
# Loop through time series data
for (ifile in 1:length(fileList))  
{
  
  # Read in as raster - here assume in order, careful!
  inrast <- raster(fileList[ifile])
  
  # Get name of data for output
  tempname <- strsplit(fileList[ifile],"/")
  tempname <- tempname[[1]][length(tempname[[1]])]
  tempname <- strsplit(tempname,".tif")
  
  # Extract year 
  tempyear <- as.numeric(substr(tempname,6,9))
  
  # If it not the last file
  if(ifile != length(fileList))
  {
    # locate index for dates - Preceeding year    
    dates.notbefore <- which((floor(dates) != (tempyear-1)))    
    r.before <- dropLayer(modisraster,dates.notbefore)                   
    meanEVI.before <- calc(r.before,fun=function(x) mean(x,na.rm = TRUE))
    
    # locate index for dates - Post year
    dates.notafter <- which((floor(dates) != (tempyear+1)))
    r.after <- dropLayer(modisraster,dates.notafter)
    meanEVI.after <- calc(r.after,fun=function(x) mean(x,na.rm = TRUE))
    
    # Create raster with EVI difference (average of entire year)
    diffr <- meanEVI.after - meanEVI.before
    one<-data.frame(x=which((getValues(inrast[[1]]==1)),arr.ind=TRUE))
    two<-data.frame(y=seq(1,ncell(diffr),1))
    two$y2<-getValues(diffr)
    merged<-merge(two,one,by.x="y",by.y="x",all=FALSE)
    merged2<-merge(merged,two[,-2],by.x="y",by.y="y",all=TRUE)
    merged2$y2[merged2$y2>0]<-NA
    diffr <- raster(inrast[[1]])
    diffr <- setValues(raster(diffr),merged2$y2)
    
    # Create raster with percent EVI difference (may/june 1 year before, may/june year of)
    percent <- diffr/meanEVI.before
    one<-data.frame(x=which((getValues(inrast[[1]]==1)),arr.ind=TRUE))
    two<-data.frame(y=seq(1,ncell(percent),1))
    two$y2<-getValues(percent)
    merged<-merge(two,one,by.x="y",by.y="x",all=FALSE)
    merged2<-merge(merged,two[,-2],by.x="y",by.y="y",all=TRUE)
    merged2$y2[merged2$y2>0]<-NA
    percent <- raster(inrast[[1]])
    percent <- setValues(raster(percent),merged2$y2)
    
    # Remove upward trends
    uptrendMask <- Which(diffr < 0)
    inrast <- inrast + uptrendMask 
    inrast[inrast == 1] <- 0 
    inrast[inrast == 2] <- 1 
    
    writeRaster(percent,filename = paste(outdir,"percentDiff",tempname,".tif",sep=""),overwrite=T)
    writeRaster(diffr,filename = paste(outdir,"rawDiff",tempname,".tif",sep=""),overwrite=T)

    
    # Write out new raster.
    writeRaster(inrast, filename=paste(outdir,"DownTrend",tempname,".tif",sep=""), overwrite=T)    
  }
  # If it is the last file
  if(ifile == length(fileList))
  {
    
    # locate index for dates - Preceeding year 
    dates.notbefore <- which((floor(dates) != (tempyear-1)))
    r.before <- dropLayer(modisraster,dates.notbefore)                   
    meanEVI.before <- calc(r.before,fun=function(x) mean(x,na.rm = TRUE))
    
    # locate index for dates - That year
    dates.notafter <- which((floor(dates) != (tempyear)))
    r.after <- dropLayer(modisraster,dates.notafter)
    meanEVI.after <- calc(r.after,fun=function(x) mean(x,na.rm = TRUE))
    
    # Create raster with EVI difference (may/june 1 year before, may/june year of)
    diffr <- meanEVI.after - meanEVI.before
    one<-data.frame(x=which((getValues(inrast[[1]]==1)),arr.ind=TRUE))
    two<-data.frame(y=seq(1,ncell(diffr),1))
    two$y2<-getValues(diffr)
    merged<-merge(two,one,by.x="y",by.y="x",all=FALSE)
    merged2<-merge(merged,two[,-2],by.x="y",by.y="y",all=TRUE)
    merged2$y2[merged2$y2>0]<-NA
    diffr <- raster(inrast[[1]])
    diffr <- setValues(raster(diffr),merged2$y2)
    
    # Create raster with percent EVI difference (may/june 1 year before, may/june year of)
    percent <- diffr/meanEVI.before
    one<-data.frame(x=which((getValues(inrast[[1]]==1)),arr.ind=TRUE))
    two<-data.frame(y=seq(1,ncell(percent),1))
    two$y2<-getValues(percent)
    merged<-merge(two,one,by.x="y",by.y="x",all=FALSE)
    merged2<-merge(merged,two[,-2],by.x="y",by.y="y",all=TRUE)
    merged2$y2[merged2$y2>0]<-NA
    percent <- raster(inrast[[1]])
    percent <- setValues(raster(percent),merged2$y2)
    
    # Remove upward trends
    uptrendMask <- Which(diffr < 0)
    inrast <- inrast + uptrendMask
    inrast[inrast == 1] <- 0
    inrast[inrast == 2] <- 1
    
    writeRaster(percent,filename = paste(outdir,"percentDiff",tempname,".tif",sep=""),overwrite=T)
    writeRaster(diffr,filename = paste(outdir,"rawDiff",tempname,".tif",sep=""),overwrite=T)
    
    # Write out new raster
    writeRaster(inrast, filename=paste(outdir,"DownTrend",tempname,".tif",sep=""), overwrite=T)    
  }    
}

##### NOT done yet - filter out the DOWNTREND ones to make cum maps

# Get Cumulative years
# ----------------------

# Probability cutoff option 
cutoff <- 0.8

# Mask files for new rasters
maskans <- FALSE

for (ifile in 1:length(fileList))  
{
  
  # Read in as raster - here assume in order, careful!
  inrast <- raster(fileList[ifile])
  
  # Get name of data for output
  tempname <- strsplit(fileList[ifile],"/")
  tempname <- tempname[[1]][length(tempname[[1]])]
  tempname <- strsplit(tempname,".tif")
  tempname <- paste(tempname,"C.tif",sep="")
  
  # If the first file - create a write raster
  if(ifile == 1)
  {
    writerast <- inrast
  }
  # If not the first file - then add to stored write raster
  if(ifile != 1)
  {
    writerast <- inrast+writerast
  }
  
  # Then apply probability cutoff and make binary
  writerast[writerast > cutoff] <- 1
  writerast[writerast < cutoff] <- 0
  
  # Make forest mask of 1 vs. NA values - should already be, but just in case
  if (maskans==TRUE)
  {
    maskrast <- raster(maskfile)
    writerast <- mask(writerast,maskrast)
  } 
  # WRite our raster
  writeRaster(writerast, filename=paste(outdir,tempname,sep=""), overwrite=T)  
}
