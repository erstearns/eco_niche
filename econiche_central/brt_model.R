jobnum <- commandArgs()[3]
repo <-  commandArgs()[4]
outpath <- commandArgs()[5]
data_loc <- commandArgs()[6]
run_date <-  commandArgs()[7]
package_lib <- commandArgs()[8]
covs <- commandArgs()[9]

## Load libraries
setwd(repo)

# Library for packages. Ensures that none of this code is dependent on the machine where the user runs the code.
.libPaths(package_lib)# Ensures packages look for dependencies here when called with library().

# Load packages
package_list <- c('seeg', 'stringr', 'reshape2', 'ggplot2', 'dplyr', 'Amelia', 'rgeos', 'data.table','raster','rgdal','INLA','seegSDM','seegMBG','plyr','sp')
for(package in package_list) {
  library(package, lib.loc = package_lib, character.only=TRUE)
}

# Load functions files
source(paste0(repo, '/econiche_central/functions.R'))                   
  
#########################################################################################

# create a list with random permutations of dat_all, 
# This is like k-folds stuff in mbg; subsample custom function
dat_all <- read.csv((paste0(data_loc, "/dat_all.csv"))) 

# Set the RNG seed
set.seed(1)

data_sample <- subsample(dat_all,
                          n = 100, #random choice
                          minimum = c(30, 30)) #half of actual data points; changed from 30,30
                          #simplify = FALSE)

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

