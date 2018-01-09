library(data.table)

# Working directory
wd = "/scratch/sharedCO/Casiraghi/Abemus_data/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/PLASMA_DEDUP/Results_new"

# k factor
sk = 2

pbem_dir = "/scratch/sharedCO/Casiraghi/Abemus_data/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/PLASMA_DEDUP/BaseErrorModel/"
# Import background pbem
tab_bg_pbem = read.delim(file = file.path(pbem_dir,"pbem_background.tsv"),as.is=T,header=T)

# Import matrix for coverages and pbem
load("/elaborazioni/sharedCO/Home_casiraghi/Prog/abemus/data/PBEsim.RData")
tab_cov_pbem = tab.list[[6]]
tab_cov_pbem = tab_cov_pbem*sk

add_class_xbg = function(pmtab,xbg){
  pmtab$CLASS.xbg = NA
  pmtab$CLASS.xbg[which(pmtab$af_control<=xbg & pmtab$bperr<=xbg & pmtab$pbem_allele<=xbg)] = 1
  pmtab$CLASS.xbg[which(pmtab$af_control<=xbg & pmtab$bperr>xbg & pmtab$pbem_allele<=xbg)] = 2
  pmtab$CLASS.xbg[which(pmtab$af_control<=xbg & pmtab$bperr>xbg & pmtab$pbem_allele>xbg)] = 3
  pmtab$CLASS.xbg[which(pmtab$af_control>xbg & pmtab$bperr>xbg & pmtab$pbem_allele>=xbg & pmtab$same_allele == 0)] = 4
  pmtab$CLASS.xbg[which(pmtab$af_control>xbg & pmtab$bperr>xbg & pmtab$pbem_allele>xbg & pmtab$same_allele == 1)] = 5
  return(pmtab)
}

ld = list.dirs(path = wd,full.names = T,recursive = F)

#pm=ld[1]
for(pm in ld){
  setwd(pm)
  message(pm)
  cpmf2.file = list.files(pm,pattern = "pmtab_F2_")
  cpmf2 = fread(input = cpmf2.file,header = T,sep='\t',colClasses = list(character=2))
  cpmf2 = data.frame(cpmf2,stringsAsFactors = F)
  # TABLE 3
  cpmf3 = add_class_xbg(pmtab = cpmf2,xbg = as.numeric(tab_bg_pbem$background_pbem))
  cpmf3$bperr[which(cpmf3$bperr > 0.2)] = 0.2
  
  to.keep = which(sapply(1:nrow(cpmf3), function(k) 
    cpmf3$af_case[k]>=tab_cov_pbem[min(which(covs>=cpmf3$cov_case[k])),
                                   min(which(afs>=cpmf3$bperr[k]))]))
  
  cpmf3$pass.filter.pbem_coverage = 0
  cpmf3$pass.filter.pbem_coverage[to.keep] = 1
  
  f3.name = gsub(basename(cpmf2.file),pattern = "_F2_",replacement = "_F3_")
  f3.name = gsub(f3.name,pattern = "\\.tsv",replacement = paste0(".",sk,".tsv"))
  write.table(cpmf3,file = f3.name,sep = '\t',col.names = T,row.names = F,quote = F)
  
}



