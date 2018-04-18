library(data.table)
wd = "/scratch/sharedCO/Casiraghi/Abemus_data/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/PLASMA_DEDUP/Results_new_af_0.01"
setwd(wd)

od = "/scratch/sharedCO/Casiraghi/Abemus_data/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/PLASMA_DEDUP/Results_new"
pms = list.files(path = od,pattern = "PM",full.names = T)

#x=pms[1]

minaf = 0.01

for(x in pms){
  message(x)
  abemus.file = list.files(x,pattern = "_F1_",full.names = T)
  f1 = fread(input = abemus.file,header = T,sep='\t',colClasses = list(character=2))
  f1 = data.frame(f1,stringsAsFactors = F)
  f2 = f1[which(f1$af_case >= minaf),]
  name=gsub(basename(abemus.file),pattern = "_F1_",replacement = "_F2_")
  dir.create(path = file.path(wd,basename(x)))
  write.table(x = f2,file = file.path(wd,basename(x),name),quote = F,col.names = T,row.names = F,sep = "\t") 
}