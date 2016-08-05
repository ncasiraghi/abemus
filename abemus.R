#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=1){
  message("\nError!\nusage: Rscript abemus.R abemus_configure.R\n")
  quit()
}
source(args[1])

abemus_ws_1 = file.path(abemusdir,"pbem.R")
abemus_ws_2 = file.path(abemusdir,"afthreshold.R")
abemus_ws_3 = file.path(abemusdir,"callsnvs.R")

if(1 %in% aws){
  cat(paste("[",Sys.time(),"]\tabemus step 1","\n"))
  system(paste("Rscript",abemus_ws_1,args[1]))
}

if(2 %in% aws){
  cat(paste("[",Sys.time(),"]\tabemus step 2","\n"))
  system(paste("Rscript",abemus_ws_2,args[1]))
}

if(3 %in% aws){
  cat(paste("[",Sys.time(),"]\tabemus step 3","\n"))
  system(paste("Rscript",abemus_ws_3,args[1]))
}