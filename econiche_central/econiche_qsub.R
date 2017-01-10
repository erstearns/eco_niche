############
## Define qsub function
############

## jobname: Name for job, must start with a letter
## code:    Filepath to code file
## hold:    Comma-separated list of jobnames to hold the job on
## pass:    List of arguments to pass on to receiving script
## slots:   Number of slots to use in job
## submit:  Should we actually submit this job?
## log:     Should this job create a log in /share/temp/sgeoutput/user/   output and errors?
## proj:    What is the project flag to be used?

## In code_dir, there must be an r_shell, python_shell, and stata_shell available for use.

## Example submission: 
# setwd(code_dir)
# pn <- "nothing"
# jobids <- NULL
# for (cc in c("AFG","USA","VIR")) {
#   if (start<=6 & end >=6) qsub(paste("am06", cc, sep="_"), paste0(code_dir,"/06_select_parameters.r"), hold = paste(pn, collapse=","), pass=list(cc), slots=4, submit=!test)
#   jobids[count] <- paste("am06", cc, sep="_")
# }
# 
# if (start<=6 & end >=6) qsub("am06b", paste0(code_dir,"/06b_combine_selected_parameters.R"),  hold=paste(jobids, collapse=","), submit =!test) #job name should be 



####### In the main script
## Make sure you setwd to a folder that has a r_shell.sh

## Delete existing results
system(paste0("rm ",data_dir,"/results_*.csv"))

njobs <- 50
for(jobnum in 1:njobs) {
  qsub(paste0("jobname_",jobnum), paste0(code_dir,"/code_file.R"), pass=list(jobnum), proj="proj_geospatial", log=T, slots=1)
}

## Check for results
check_loc_results(c(1:njobs),data_dir,prefix="results_",postfix=".csv")


########## Within the child script
jobnum <- commandArgs()[3]
data <- fread(paste0(data_dir,"/simulation_",jobnum,".csv"))

.
.
.

write.csv(paste0(data_dir,"/results_",jobnum,".csv"))



qsub <- function(jobname, code, hold=NULL, pass=NULL, slots=1, submit=F, log=T, intel=F, proj = "proj_mortenvelope") { #do not need to define hold unless multi-step; pass parameters
  user <- Sys.getenv("USER") # Default for linux user grab. "USERNAME" for Windows
  # choose appropriate shell script 
  if(grepl(".r", code, fixed=T) | grepl(".R", code, fixed=T)) shell <- "r_shell.sh" else if(grepl(".py", code, fixed=T)) shell <- "python_shell.sh" else shell <- "stata_shell.sh" 
  # set up number of slots
  if (slots > 1) { 
    slot.string = paste(" -pe multi_slot ", slots, sep="")
  } 
  # set up jobs to hold for 
  if (!is.null(hold)) { 
    hold.string <- paste(" -hold_jid \"", hold, "\"", sep="")
  } 
  # set up arguments to pass in 
  if (!is.null(pass)) { 
    pass.string <- ""
    for (ii in pass) pass.string <- paste(pass.string, " \"", ii, "\"", sep="")
  }  
  # construct the command 
  sub <- paste("qsub",
               if(log==F) " -e /dev/null -o /dev/null ",  # don't log (if there will be many log files)
               if(log==T) paste0(" -e /share/temp/sgeoutput/",user,"/errors -o /share/temp/sgeoutput/",user,"/output "),
               if(intel==T) paste0(" -l hosttype=intel "),
               if(proj != "") paste0(" -P ",proj," "),
               if (slots>1) slot.string, 
               if (!is.null(hold)) hold.string, 
               " -N ", jobname, " ",
               shell, " ",
               code, " ",
               if (!is.null(pass)) pass.string, 
               sep="")
  # submit the command to the system
  if (submit) {
    system(sub) 
  } else {
    cat(paste("\n", sub, "\n\n "))
    flush.console()
  } 
} 