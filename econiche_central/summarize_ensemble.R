########################################################################################
#Summarizing the BRT ensemble 
#There are three helper functions to help us summarize the BRT ensemble: 
#1. getRelInf
#2. getEffectPlots
#3. combinePreds
########################################################################################
repo <-  commandArgs()[3]
outpath <- commandArgs()[4]
data_loc <- commandArgs()[5]
run_date <-  commandArgs()[6]
package_lib <- commandArgs()[7]
data_dir_model <- commandArgs()[8]
data_dir_stats <- commandArgs()[9]

## Load libraries
setwd(repo)

# Library for packages. Ensures that none of this code is dependent on the machine where the user runs the code.
.libPaths(package_lib)# Ensures packages look for dependencies here when called with library().

# Load packages
package_list <- c('car', 'MASS', 'stringr', 'reshape2', 'ggplot2', 'plyr', 'dplyr', 'rgeos', 'data.table','raster','rgdal', 'seegSDM','sp')
for(package in package_list) {
  library(package, lib.loc = package_lib, character.only=TRUE)
}

# Load functions files
source(paste0(repo, '/econiche_central/functions.R'))                   

#########################################################################################
#Bring in model and stats and make into 2 lists
model_filenames <- list.files(path = (data_dir_model))
stat_filenames <- list.files(path = (data_dir_stats))

#Need something else for reading in Rdata files - from Daniel
#Function for loading data from a specific path

fetch_from_rdata = function(file_location, item_name, use_grep = F){
  load(file_location)
  if(use_grep){
    return(mget(grep(item_name, ls(), value=T)))
  }else{
    return(get(item_name))
  }
}

model_list <- lapply(model_filenames, function(x) fetch_from_rdata(x, 'model', F)) #change file location



#Reading in files and making into lists
model_list <- lapply(paste0(data_dir_model, model_filenames), fread)
stat_lis <- lapply(paste0(data_dir_stats, stat_filenames), fread)

#get all the prediction results - lapply to extract the predictions into a list then coerce the list into a rasterbrick
preds <- brick(lapply(model_list, '[[', 4)) #4th component likely the prediction raster layer

#3. combinePreds: combines the prediction maps (on the probability scale) from multiple models and returns rasters giving the mean, median and quantiles
#                 of the ensemble predictions. unlike the previous two functions, combinePreds needs a RasterBrick or RasterStack object with each layers 
#                 giving a single prediction. So we need to create one of these before we can use combinePreds. Note that we can also run combinePreds in 
#                 parallel to save some time if the rasters are particularly large.
# Run combinePreds - could summarise the predictions in parallel
preds_sry <- combinePreds(preds)

names(preds_sry) <- c('mean',
                      'median',
                      'lowerCI',
                      'upperCI')

# save the prediction summary
writeRaster(preds_sry,
            file = paste0(outpath,
                          '/schisto_pred_summary_', run_date),
            format = 'GTiff',
            overwrite = TRUE)

# convert the stats list into a matrix using the do.call function
stats <- do.call("rbind", model_stats)

# save them
write.csv(stats,
          paste0(outpath, '/stats_', run_date, ".csv"))


#1. getRelInf: combines the models to get summaries of the relative influence of the covariates across all the models. 
#              It returns a matrix of means and quantiles, and optionally produces a boxplot.
relinf <- getRelInf(model_list, plot = TRUE)

#Save the relative influence scores
write.csv(relinf,
          file = paste0(outpath,
                        '/relative_influence', run_date, '.csv'))

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
effect <- getEffectPlots(model_list, plot = FALSE)

# calculate uncertainty
preds_sry$uncertainty <- preds_sry[[4]] - preds_sry[[3]]

# write the mean prediction and uncertainty rasters as Geo-Tiffs
writeRaster(preds_sry$mean, paste0(outpath, '/prediction_map_', run_date, '.tif'), format = 'GTiff')
writeRaster(preds_sry$uncertainty, paste0(outpath, '/uncertainty_map_', run_date, '.tif'), format = 'GTiff')

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
par(mfrow = c(3, 4),
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
