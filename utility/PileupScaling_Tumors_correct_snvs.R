library(data.table)

args <- commandArgs(trailingOnly = TRUE)

if(length(args)!=1){
  message("\nERROR! Required arguments:\n\n1. main folder i.e. /scratch/sharedCO/Casiraghi/InSilicoData_Tumors/HALO_byChrom_c10_model\n")
  quit()
}

#main.folder = "/scratch/sharedCO/Casiraghi/InSilicoData_Tumors/HALO_byChrom_c10_model"
main.folder = args[1]
samples = list.dirs(main.folder,full.names = T,recursive = F)
samples = grep(samples,pattern = "adm",value = T)

get.ALT = function(m){
  REF = m["ref"]
  cdz = m[setdiff(c("A","C","G","T"),REF)]  
  ALT = names(cdz[which(cdz == max(cdz))])
  if(length(ALT)>1){
    ALT = "N"
  }
  return(as.character(ALT))
}

#s<-samples[1]
for(s in samples){
  message(basename(s))
  p = list.files(path = file.path(s,"pileup"),full.names = T)
  # generate and write snvs file
  setwd(file.path(s,"snvs"))
  for(x in p){
    message(x)
    p.scaled = fread(input = x,stringsAsFactors = F,showProgress = F,header = T,na.strings = "",colClasses = list(character=10),data.table = F)
    s.scaled = p.scaled[which(p.scaled$af > 0),,drop=F]
    s.scaled$Ars = 0 
    s.scaled$Crs = 0
    s.scaled$Grs = 0
    s.scaled$Trs = 0
    s.scaled = s.scaled[,c("chr","pos","ref","A","C","G","T","af","cov","Ars","Crs","Grs","Trs","dbsnp")]
    s.scaled$alt = apply(s.scaled,MARGIN = 1,FUN = get.ALT)
    s.scaled$af[which(s.scaled$alt == "N")] <- 0
    s.scaled = s.scaled[,c("chr","pos","ref","alt","A","C","G","T","af","cov","Ars","Crs","Grs","Trs","dbsnp")]
    write.table(x = s.scaled,file = gsub(basename(x),pattern = "\\.pileup$",replacement = ".snvs"),col.names = T,row.names = F,sep = "\t",quote = F,na = "")
  }
}
