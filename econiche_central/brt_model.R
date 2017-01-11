
jobnum <- commandArgs()[3]
repo <- 'C:/Users/stearns7/OneDrive - UW Office 365/ecological_niche_models'
outpath <- (paste0(repo, '/code/schisto/output/'))
source('econiche_central/functions.R')                   



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
save(model, stats, file = paste0(outpath, jobnum,".RData"))




#things to be reset - all objects called
#repo
#outpath
#data_loc
#covs

