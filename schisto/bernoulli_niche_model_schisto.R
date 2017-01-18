######################################
#Original code author: Nick Golding
#Modified and updated by: Erin Stearns
#Date: 1/5/17
#Code intent: Bernoulli niche model
#####################################

## notes.
#1st iteration, use only point data
#Ellen cleaned data in such a way to only include point data
########################################################################################

# clear workspace
#rm(list = ls())

########################################################################################
#Setting up

# - Initially examining environment and indicating file paths dependent on OS - Linux vs Windows
# - Calling 'parallel' and specifying cores 
# - Setting data location, output path,  working directory and library
# - Sourcing external function scripts
# - Loading packages
# - Creating a time stamp for model run
########################################################################################

if (Sys.info()[1] == "Linux"){
  j <- "/home/j"
  h <- paste0("/home/",Sys.info()[6]) # what is this 6?
  package_lib <- paste0(j,'/temp/stearns7/packages_cl') # Library for packages on cluster. Ensures that none of this code is dependent on the machine where the user runs the code.
  repo <- '/share/code/geospatial/stearns7/eco_niche' ## Set repo location on cluster
}else{
  j <- "J:"
  h <- "H:"
  package_lib <- paste0(j,'/temp/stearns7/packages') #library for packages locally
  repo <- paste0(j, '/temp/stearns7/eco_niche') ## Set repo location locally
}

library('parallel')
slots = 40  
#if i requested 40 slots, how many cores does this mean?
cores_to_use = ifelse(grepl('Intel', system("cat /proc/cpuinfo | grep \'name\'| uniq", inter = T)), floor(slots * .86), floor(slots*.64))


## Set data location
data_loc <- (paste0(j, '/temp/stearns7/schisto/data/eco_niche_data'))

# Set output path
outpath <- (paste0(data_loc, '/output'))

## Load libraries
setwd(repo)

# Library for packages. Ensures that none of this code is dependent on the machine where the user runs the code.
.libPaths(package_lib)# Ensures packages look for dependencies here when called with library().

# Load functions files
source(paste0(repo, '/econiche_central/functions.R'))                   
source(paste0(repo, '/econiche_central/econiche_qsub.R'))  
source(paste0(repo, '/econiche_central/check_loc_results.R'))  

# Load packages
package_list <- c('car', 'MASS', 'stringr', 'reshape2', 'ggplot2', 'plyr', 'dplyr', 'rgeos', 'data.table','raster','rgdal', 'seegSDM','sp')
for(package in package_list) {
  library(package, lib.loc = package_lib, character.only=TRUE)
}


## Create run date in correct format - calls make_time_stamp function from 'functions' - copied from Nick Graetz's in 'prep_functions' for MBG code
time_stamp <- TRUE 
run_date <- make_time_stamp(time_stamp)

########################################################################################
# Preparing the data

# - Loading the covariate raster brick - created prior to running this script
# - Loading occurrence data - cleaned and formatted prior to running this script
# - Generating background data using the pixel grid of the aridity layer
# - Coercing background data into a data frame and assigning each set of coordinates a value of '0' for outbreak_id/prev value
# - Joining the occurrence and background data into single data frame and creating new binary var to indicate presence/absence (PA) then dropping 'outbreak_id'
# - Extracting values for each covariate at every data point (occurrence and background)
# - Joining these covariate vals to the greater data frame (w/ occurrence and background data)
# - Omitting all null values - coastlines issues likely
########################################################################################

# Load covariate raster brick here (created ahead of time)
covs <- brick(paste0(data_loc, "/covariates/schisto_covs.grd"))
print('Loading covariate brick')

# Occurrence data - schisto point data; will need to change for each species, currently mansonia
occ <- read.csv(file = (paste0(data_loc, '/man_fin.csv')))
print('Loading occurrence data')
occ <- occ[,c(2:4)]

# Generate pseudo-absence data according to the aridity surface and suppress weighting (prob=FALSE) so as to not weight by aridity pixel values
aridity <- raster(paste0(data_loc, "/covariates/aridity_annual.tif"))
print('Loading grid for background point generation')

bg <- bgSample(aridity, # Weighting grid - population in this case, custom function defined in github 
               n = 2500, # Background data points desired
               prob = FALSE, # Set to FALSE so doesn't weight by raster specified above
               replace = TRUE,
               spatial = FALSE)
print('Generating background data')

colnames(bg) <- c('long', 'lat') 
bg <- data.frame(bg)
print('Making background data into a dataframe')

# Add an outbreak id to this
bg$outbreak_id <- 0

# Combine the occurrence and background records
dat <- rbind(cbind(PA = rep(1, nrow(occ)),
                   occ[, c('long', 'lat', 'outbreak_id')]),
             cbind(PA = rep(0, nrow(bg)),
                   bg))
print('Combining occurrence and background records')
 
#need to drop 'outbreak_id'
dat <- dat[,c(1:3)]

# Get the covariate values for every data point - pseudo and actual
dat_covs <- extract(covs, dat[, 2:3])  #this is where we will need to update to include years/subset by years, extract to those subsets, then re-merge
print('Getting covariate values for every data point - background and observed')

# Then add them
dat_all <- cbind(dat, dat_covs)
print('Adding extracted covariate values to the occurrence and background records dataframe')

# Remove NAs - may need to check this - coastline issues, especially with island data
dat_all <- na.omit(dat_all)
print('Omitting all null values from dataframe')

#write.csv(dat_all, file = (paste0(data_loc, "/dat_all.csv"))) ## may be able to remove if do not run the below

########################################################################################
# Running a BRT ensemble in parallel

# - Extract a random subset of n records from data, ensuring that there are at least minimum presence and absence records in the subset. Default assumes first 
#          column is presence/absence column
# - Running the BRT ensemble in parallel
# - Getting model validation statistics for model objects - Given an object returned by runBRT, extract devBern, rmse, auc, Kappa, sensitivity and specificity and proportion correctly classified (pcc) validation statistics - calculated using either the PresenceAbsence or seegSDM functions. Note that auc is calculated with a seegSDM clone of the auc function in PresenceAbsence but in which worse-than-random AUC scores are not inverted.
########################################################################################
#create multiple versions of the dataset via sample; right now: takes 25 samples of 800 obs with at least 30 presence and 30 absence points
data_sample = lapply(1:25, function(x) subsample(dat_all, 800, minimum= c(30,30)))

#Run the brts
models <- mclapply(data_sample, function(x) runBRT(x,
                                                   gbm.x = 4:ncol(x),
                                                   gbm.y = 1,
                                                   pred.raster = covs, #brick
                                                   gbm.coords = 2:3,
                                                   wt = function(PA) ifelse(PA == 1, 1, sum(PA) / sum(1 - PA))),mc.cores = cores_to_use ) #mc.cores defined above

########################################################################################
#Getting model validation statistics for model objects 

# - Given an object returned by runBRT, extract devBern, rmse, auc, Kappa, sensitivity and specificity and proportion correctly classified (pcc) 
#   validation statistics - calculated using either the PresenceAbsence or seegSDM functions. Note that auc is calculated with a seegSDM clone of 
#   the auc function in PresenceAbsence but in which worse-than-random AUC scores are not inverted.

#Summarizing the BRT ensemble 
#There are three helper functions to help us summarize the BRT ensemble: 
#1. getRelInf
#2. getEffectPlots
#3. combinePreds
########################################################################################
#get model stats
model_stats <- suppressWarnings(lapply(models, function(x) getStats(x)))

#get all the prediction results - lapply to extract the predictions into a list then coerce the list into a rasterbrick
preds <- brick(lapply(models, '[[', 4)) #4th component likely the prediction raster layer

#3. combinePreds: combines the prediction maps (on the probability scale) from multiple models and returns rasters giving the mean, median and quantiles
#                 of the ensemble predictions. unlike the previous two functions, combinePreds needs a RasterBrick or RasterStack object with each layers 
#                 giving a single prediction. So we need to create one of these before we can use combinePreds. Note that we can also run combinePreds in 
#                 parallel to save some time if the rasters are particularly large.
# Run combinePreds - could summarise the predictions in parallel
preds_sry <- combinePreds(preds)

# plot the resulting maps - interactive sessions
plot(preds_sry, zlim = c(0, 1))

# convert the stats list into a matrix using the do.call function
stats <- do.call("rbind", model_stats)
print('Converting stats list into a matrix')
head(stats)

# and produce a boxplot of a few imnportant statistics
boxplot(stats[, 3:7], col = 'grey', ylim = c(0, 1))

# save them
write.csv(stats,
          paste0(outpath, '/stats_', run_date, ".csv"))


names(preds_sry) <- c('mean',
                      'median',
                      'lowerCI',
                      'upperCI')

# save the prediction summary
writeRaster(preds_sry,
            file = paste0(outpath,
                          '/Schisto_', run_date),
            format = 'GTiff',
            overwrite = TRUE)

#what's important after this?


############################################################################################################################
#We now have a list of model outputs, each of which contains the fitted model, predictions and information for plotting. 
#We can pull out individual model runs and plot the predictions.
############################################################################################################################
# par(mfrow = c(1, 2))
# plot(models[[1]]$pred, main = 'run 1', zlim = c(0, 1))
# plot(models[[25]]$pred, main = 'run 25', zlim = c(0, 1))

############################################################################################################################
#Summarizing the BRT ensemble 
#There are three helper functions to help us summarize the BRT ensemble: 

#1. getRelInf: combines the models to get summaries of the relative influence of the covariates across all the models. 
#              It returns a matrix of means and quantiles, and optionally produces a boxplot.
relinf <- getRelInf(models, plot = TRUE)
print("Summarizing relative influence of covariates across all models")

#Save the relative influence scores
write.csv(relinf,
          file = paste0(outpath,
                        '/relative_influence', run_date, '.csv'))
print("saving relative influence scores as a csv file")

# plot the risk map
png(paste0(outpath,
           '/Schisto_', run_date, '.png'),
    width = 2000,
    height = 2000,
    pointsize = 30)

par(oma = rep(0, 4),
    mar = c(0, 0, 0, 2))

plot(preds_sry[[1]],
     axes = FALSE,
     box = FALSE)

points(dat[dat$PA == 1, 2:3],
       pch = 16,
       cex = 1,
       col = 'blue')

dev.off()

#2. getEffectPlots: performs a similar operation for the effect plots, optionally plotting mean effects with uncertainty intervals.
par(mfrow = c(1, 3))
effect <- getEffectPlots(models, plot = FALSE)

###########################################################################################################################
# Optional - We can create a simple map of prediction uncertainty by subtracting the lower from the upper quantile.
##########################################################################################################################

# calculate uncertainty
preds_sry$uncertainty <- preds_sry[[4]] - preds_sry[[3]]

# plot mean and uncertainty
par(mfrow = c(1, 2))

# plot mean
plot(preds_sry$mean,
     zlim = c(0, 1),
     main = 'mean')

# and uncertainty
plot(preds_sry$uncertainty,
     col = topo.colors(100),
     main = 'uncertainty')

# write the mean prediction and uncertainty rasters as Geo-Tiffs
writeRaster(preds_sry$mean, paste0(outpath, '/prediction_map_', run_date, '.tif'), format = 'GTiff')
writeRaster(preds_sry$uncertainty, paste0(outpath, '/uncertainty_map_', run_date, '.tif'), format = 'GTiff')
###########################################################################################################################



##########################################################################################################################3
# Plot marginal effect curves
########################################################################################


# get the order of plots (all except relinf)! - customize list to covs of interest
order <- match(rownames(relinf), names(covs))

# Set up x axis labels and titles
short_names <- c(
  'precipitation',
  'temperature',
  'aridity',
  'elevation',
  'evi',
  'land_cover',
  'ses',
  'urban',
  'urban_access',
  'dist_fresh',
  'tcb',
  'tcw',
  'Schistosomiasis distribution')


units <- c(
  'mm',
  'degrees celsius',
  'index',
  'meters',
  'index (0-1)',
  '16 classes',
  'Gross cell product',
  'Urban/Rural(3 classes)',
  'Travel time to nearest settlement',
  'meters',
  'dimensionless',
  'dimensionless',
  'suitability index')



# Set up device 
png(paste0(outpath,
           '/effects_', run_date, '.png'),
    width = 3000,
    height = 3000,
    pointsize = 60)

# set up multi panels and margins
par(mfrow = c(3, 3),
    mar = c(5, 2, 4, 2) + 0.1,
    oma = c(0, 3, 0, 0))

# loop through plots
for (i in 1:length(effect)) {
  
  # extract summary stats
  df <- effect[[order[i]]][, 1:4]
  
  # pick y axis
  if (i %% 3 == 1) {
    ylab = 'marginal effect'
  } else {
    ylab = ''
  }
  
  
  # set up empty plotting region
  plot(df[, 2] ~ df[, 1],
       type = 'n',
       ylim = c(-3.5, 2.5),
       ylab = '',
       xlab = '')
  
  # add the 95% CIs
  polygon(x = c(df[, 1], rev(df[, 1])),
          y = c(df[, 3], rev(df[, 4])),
          border = NA,
          col = grey(0.7))
  
  # add the mean line
  lines(df[, 2] ~ df[, 1],
        lwd = 5,
        col = grey(0.2))
  
  # y axis lable (only on left hand column)
  title(ylab = ylab,
        cex.lab = 1.2,
        col.lab = grey(0.3),
        xpd = NA,
        line = 2.5)
  
  # x-axis label
  title(xlab = units[order[i]],
        cex.lab = 1.2,
        col.lab = grey(0.3),
        line = 2.5)
  
  # title
  title(main = short_names[order[i]],
        line = 1.5,
        cex.main = 1.2)
  
  # relative contribution inset
  mtext(text = round(relinf[i, 1] / 100, 2),
        side = 3,
        line = -2,
        adj = 0.07,
        col = grey(0.5))
  
}

dev.off()