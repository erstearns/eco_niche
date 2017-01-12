
jobnum <- commandArgs()[3]
repo <- 'J:/temp/stearns7/eco_niche/'  
outpath <- (paste0(data_loc, 'output/'))
source(paste0(repo, '/econiche_central/functions.R'))                     

## Set repo location 
repo <- 'J:/temp/stearns7/eco_niche/'  

## Set data location
data_loc <- 'J:/temp/stearns7/schisto/data/eco_niche_data/'

## Load libraries
setwd(repo)

####
root <- paste0(j, "temp/stearns7/eco_niche/")  ###################  ?????????????????
OR
root <- ifelse(Sys.info()[1]=="Windows", "J:/", "/home/j/") ###################### ??????????????????
####



package_lib <- paste0(root,'/temp/stearns7/packages') # Library for packages. Ensures that none of this code is dependent on the machine where the user runs the code.
.libPaths(package_lib)# Ensures packages look for dependencies here when called with library().

# Load packages
package_list <- c('seeg', 'stringr', 'reshape2', 'ggplot2', 'dplyr', 'Amelia', 'rgeos', 'data.table','raster','rgdal','INLA','seegSDM','seegMBG','plyr','sp')
for(package in package_list) {
  library(package, lib.loc = package_lib, character.only=TRUE)
}

# Load functions files
source(paste0(repo, '/econiche_central/functions.R'))                   
  
## Create run date in correct format - calls make_time_stamp function from 'functions' - copied from Nick Graetz's in 'prep_functions' for MBG code
run_date <- make_time_stamp(time_stamp)

# Set output path
outpath <- (paste0(data_loc, 'output/'))

#########################################################################################

# create a list with random permutations of dat_all, 
# This is like k-folds stuff in mbg; subsample custom function
dat_all <- read.csv((paste0(data_loc, "dat_all.csv"))) 

# Set the RNG seed
set.seed(1)

data_sample <- subsample(dat_all,
                          n = 100, #random choice
                          minimum = c(30, 30), #half of actual data points
                          simplify = FALSE)

model <- runBRT(data_sample,
          gbm.x = 4:ncol(data_sample),
          gbm.y = 1,
          pred.raster = covs, #brick
          gbm.coords = 2:3,
          wt = function(PA) ifelse(PA == 1, 1, sum(PA) / sum(1 - PA)))

stats <- getStats(model)

# Output model results
save(model, stats, file = paste0(outpath, run_date, jobnum,".RData"))




#things to be reset - all objects called
#repo
#outpath
#data_loc
#covs

