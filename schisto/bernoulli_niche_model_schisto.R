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
########################################################################################
if (Sys.info()[1] == "Linux"){
  j <- "/home/j"
  h <- paste0("/home/",Sys.info()[6]) # what is this 6?
}else{
  j <- "J:"
  h <- "H:"
}
## Set repo location 
repo <- (paste0(j,'/temp/stearns7/eco_niche')) 

## Set data location
data_loc <- (paste0(j, '/temp/stearns7/schisto/data/eco_niche_data'))

## Load libraries
setwd(repo)

package_lib <- paste0(j,'/temp/stearns7/packages') # Library for packages. Ensures that none of this code is dependent on the machine where the user runs the code.
.libPaths(package_lib)# Ensures packages look for dependencies here when called with library().

# Load functions files
source(paste0(repo, '/econiche_central/functions.R'))                   
source(paste0(repo, '/econiche_central/econiche_qsub.R'))  
source(paste0(repo, '/econiche_central/check_loc_results.R'))  

# Load packages
package_list <- c('car', 'MASS', 'seeg', 'stringr', 'reshape2', 'ggplot2', 'dplyr', 'Amelia', 'rgeos', 'data.table','raster','rgdal','INLA','seegSDM','seegMBG','plyr','sp')
for(package in package_list) {
  library(package, lib.loc = package_lib, character.only=TRUE)
}


## Create run date in correct format - calls make_time_stamp function from 'functions' - copied from Nick Graetz's in 'prep_functions' for MBG code
time_stamp <- TRUE 
run_date <- make_time_stamp(time_stamp)

# Set output path
outpath <- (paste0(data_loc, '/output'))

########################################################################################
# Preparing the data
########################################################################################

# Load covariate raster brick here (created ahead of time)
covs <- brick(paste0(data_loc, "/covariates/schisto_covs.grd"))
print('Loading covariate brick')

# Occurrence data - schisto point data; will need to change for each species, currently mansonia
occ <- read.csv(file = (paste0(data_loc, '/man_fin.csv')))
print('Loading occurrence data')


# Generate pseudo-absence data according to the aridity surface and suppress weighting (prob=FALSE) so as to not weight by aridity pixel values
aridity <- raster(paste0(data_loc, "/covariates/aridity_annual.tif"))
print('Loading grid for background point generation')

bg <- bgSample(aridity, # Weighting grid - population in this case, custom function defined in github 
               n = 10000, # Background data points desired
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

# Get the covariate values for every data point - pseudo and actual
dat_covs <- extract(covs, dat[, 2:3])  #this is where we will need to update to include years/subset by years, extract to those subsets, then re-merge
print('Getting covariate values for every data point - background and observed')

# Then add them
dat_all <- cbind(dat, dat_covs)
print('Adding extracted covariate values to the occurrence and background records dataframe')

# Remove NAs - may need to check this - coastline issues, especially with island data
dat_all <- na.omit(dat_all)
print('Omitting all null values from dataframe')

write.csv(dat_all, file = (paste0(dat_loc, "dat_all.csv")))
###output as ref csv for random permutations and create new script to randomly sample and call from qsub; and set seed in qsub call

########################################################################################
# Preparing to run models
########################################################################################
njobs <- 50 #no. of bootstraps; determines number of model runs - reduced to 50 for first run# 
########################################################################################
#Parallelizing
########################################################################################
parallel_script <- (paste0(repo,"/econiche_central/brt_model.R"))

for(jobnum in 1:njobs) {
  qsub(paste0("jobname_",jobnum), parallel_script, pass=list(jobnum, repo, outpath, data_loc, run_date, package_lib, covs), proj="proj_geospatial", log=T, slots=1)
}

## Check for results - makes sure models running and allows time for them to run
Sys.sleep(600)
check_loc_results(c(1:njobs),data_dir,prefix="results_",postfix=".csv")

########################################################################################
#Bring in model and stats and make into 2 lists; run a loop around job_num from 1:njobs and load, then kick into a list
model_filenames <- list.files(path=paste0(outpath, "/"), pattern = "model")
stat_filenames <- list.files(path=paste0(outpath, "/"),  pattern = "stat")

#Reading in files and making into lists
model_list <- lapply(paste0(outpath,"/",model_filenames), fread)
stat_lis <- lapply(paste0(outpath,"/",stat_filenames), fread)

# summarise all the ensembles
preds <- stack(lapply(model_list, '[[', 4)) #4th component likely the prediction raster layer

# summarise the predictions in parallel
preds_sry <- combinePreds(preds)

########################################################################################

# convert the stats list into a matrix using the do.call function
stats <- do.call("rbind", stat_lis)

# save them
write.csv(stats,
          paste0(outpath, '/stats.csv'))

names(preds_sry) <- c('mean',
                      'median',
                      'lowerCI',
                      'upperCI')

# save the prediction summary
writeRaster(preds_sry,
            file = paste0(outpath,
                          'Schisto'),
            format = 'GTiff',
            overwrite = TRUE)

# save the relative influence scores - look up function; functional form responses,
relinf <- getRelInf(model_list)
write.csv(relinf,
          file = paste0(outpath,
                        'relative_influence.csv'))


# plot the risk map
png(paste0(outpath,
           'Schisto.png'),
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


########################################################################################
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


effect <- getEffectPlots(model_list)

# Set up device 
png(paste0(outpath,
           'effects.png'),
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