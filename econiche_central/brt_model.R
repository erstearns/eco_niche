########## Within the child script
jobnum <- commandArgs()[3]
data <- fread(paste0(data_dir,"/simulation_",jobnum,".csv"))

.
.
.

write.csv(paste0(data_dir,"/results_",jobnum,".csv"))