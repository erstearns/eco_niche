########## Within the child script
jobnum <- commandArgs()[3]
data <- fread(paste0(data_dir,"/simulation_",jobnum,".csv"))

sfInit(parallel = TRUE, cpus = ncpu)
sfLibrary(seegSDM)

model_list <- sfLapply(data_list,
                       runBRT,
                       gbm.x = 4:ncol(data_list[[1]]),
                       gbm.y = 1,
                       pred.raster = covs,
                       gbm.coords = 2:3,
                       wt = function(PA) ifelse(PA == 1, 1, sum(PA) / sum(1 - PA)))

# get cv statistics in parallel
stat_lis <- sfLapply(model_list, getStats)

# summarise all the ensembles
preds <- stack(lapply(model_list, '[[', 4))

# summarise the predictions in parallel
preds_sry <- combinePreds(preds)

# stop the cluster
sfStop()

write.csv(paste0(data_dir,"/results_",jobnum,".csv"))



