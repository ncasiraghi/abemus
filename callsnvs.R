#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=1){
  message("\nError!\nusage: Rscript callsnvs.R abemus_configure.R\n")
  quit()
}
source(args[1])
abemus_functions = file.path(abemusdir,"functions.R")
source(abemus_functions)

cat(paste("[",Sys.time(),"]\tDetection of somatic SNVs in case samples","\n"))

if(!file.exists(file.path(outdir, "Results"))){
  results_folder_name = "Results"
  dir.create(file.path(outdir,results_folder_name), showWarnings = T)
} else {
  results_folder_name = paste0("Results_",gsub(as.character(Sys.time()),pattern = " ",replacement = "_"))
  dir.create(file.path(outdir,results_folder_name), showWarnings = T)
}

setwd(file.path(outdir,results_folder_name))
timestart <- proc.time()

# Import sample info file
TableSif = read.delim(filesif,as.is=T)

# remove NAs and keep only case samples having matched germline samples
nas = which(is.na(TableSif$plasma.bam) | is.na(TableSif$germline.bam))
if(length(nas)>0){TableSif = TableSif[-nas,]}

# get chromosomes
chromosomes = read.delim(file = file.path(targetbp,"bpcovered.tsv"),as.is=T)
chromosomes = sort(paste0("chr",chromosomes[-nrow(chromosomes),1]))

# Apply filters to define a set of putative SNVs
cat(paste("[",Sys.time(),"]\tApply basic filters"),"\n")
fpam = data.frame(AFbycov = as.character(AFbycov),
                  #spec = as.character(spec),
                  spec = as.character(basename(minaf_corrected_table)),
                  mincov = as.character(mincov),
                  minalt = as.character(minalt),
                  mincovgerm = as.character(mincovgerm),
                  maxafgerm = as.character(maxafgerm),
                  filtering_date = Sys.time(),
                  spec_extended = as.character(minaf_corrected_table),
                  stringsAsFactors = F)
write.table(fpam,file = "filtering_criteria.tsv",col.names = T,row.names = F,sep="\t",quote = F)

# Import background pbem
tab_bg_pbem = read.delim(file = file.path(pbem_dir,"pbem_background.tsv"),as.is=T,header=T)

# Import matrix for coverages and pbem
load("/elaborazioni/sharedCO/Home_casiraghi/Prog/abemus/data/PBEsim.RData")
tab_cov_pbem = tab.list[[6]]

# Check if data tables exists
if(file.exists(file.path(controls_dir,"datathreshold.RData"))){
  cat(paste("[",Sys.time(),"]\tlooking for data table with AF thresholds:",file.path(controls_dir,"datathreshold.RData"),"[ ok ]"),"\n")
  load(file = file.path(controls_dir,"datathreshold.RData"))
} else {
  cat(paste("[",Sys.time(),"]\tlooking for data table with AF thresholds:",file.path(controls_dir,"datathreshold.RData"),"[ NOT found ]"),"\n")
  quit()
}

# Check if corrected AF threshold exists
if(file.exists(file.path(minaf_corrected_table))){
  cat(paste("[",Sys.time(),"]\tlooking for data table with AF thresholds corrected:",file.path(minaf_corrected_table),"[ ok ]"),"\n")
  load(file = file.path(minaf_corrected_table))
} else {
  cat(paste("[",Sys.time(),"]\tlooking for data table with AF thresholds corrected:",file.path(minaf_corrected_table),"[ NOT found ]"),"\n")
  quit()
}

apply_AF_filters <- function(chrpmF1,AFbycov,af.threshold.table,minaf,mybreaks,mc.cores){
  if (AFbycov == FALSE & !is.numeric(AFbycov)){
    chrpmF1[,'af_threshold'] <- minaf
  } else if (is.numeric(AFbycov)){
    chrpmF1[,'af_threshold'] <- AFbycov
  } else if (AFbycov == TRUE & !is.numeric(AFbycov)){
    thresholds = as.numeric(af.threshold.table[,-1])
    af_filter_by_coverage <- function(y,thresholds,chrpmF1){
      this = chrpmF1[y,,drop=F]
      minaf_covth = thresholds[findInterval(this$cov_case,mybreaks)]
      this[,'af_threshold'] <- minaf_covth
      return(this)
    }
    out = mclapply(seq(1,nrow(chrpmF1),1),af_filter_by_coverage,thresholds=thresholds,chrpmF1=chrpmF1,mc.cores = mc.cores)
    chrpmF1 = fromListToDF(out)
  }
  return(chrpmF1)
}

filter = function(i,chromosomes,patient_folder,plasma.folder,germline.folder,out1,out2){
  chrom = chromosomes[i]
  # create chromsome sub-folder
  chromdir = file.path(patient_folder,chrom)
  dir.create(chromdir)
  setwd(chromdir)
  # import files
  plasma_snvs = list.files(file.path(plasma.folder,"snvs"),pattern = paste0("_",chrom,".snvs"),full.names = T)
  n.rows.plasma_snvs = as.numeric(unlist(strsplit(trimws(x = system(paste("wc -l",plasma_snvs),intern = T),which = "left"),split = " "))[[1]])
  if(n.rows.plasma_snvs == 1){
    return()
  }         
  snvs = fread(plasma_snvs,stringsAsFactors = F,showProgress = F,header = F,skip = 1,na.strings = "",colClasses = list(character=3,4,15),verbose = F)
  snvs = unique(snvs)
  snvs = data.frame(snvs)
  names(snvs)=c("chr","pos","ref","alt","A","C","G","T","af","cov","Ars","Crs","Grs","Trs","dbsnp")
  if(nrow(snvs)==0){
    return()
  }
  # F1) Custom basic filters [ in plasma/tumor ]
  out = mclapply(seq(1,nrow(snvs),1),CheckAltReads,snvs=snvs,mc.cores = mc.cores)
  snvs = fromListToDF(out)
  snvs = snvs[which(snvs$af > 0 ),]
  snvs = snvs[which(snvs$cov >= mincov ),]
  snvs = snvs[which(snvs$cov.alt >= minalt ),]
  snvs = unique(snvs)
  if(nrow(snvs)==0){
    return()
  }
  # print filtered positions and grep these pos only from pileup file of germline sample
  cat(unique(snvs$pos),sep = "\n",file = "postogrep.txt",append = F)
  controlfolder_pileup <- list.files(file.path(germline.folder,"pileup"),pattern = paste0("_",chrom,".pileup"),full.names = T)
  cmd = paste("awk -F'\t' '{if (FILENAME == \"postogrep.txt\") { t[$1] = 1; } else { if (t[$2]) { print }}}' postogrep.txt",controlfolder_pileup,"> filtered.germline.pileup.txt")
  system(cmd)
  ctrl.pileup = fread("filtered.germline.pileup.txt",stringsAsFactors = F,showProgress = T,header = F,na.strings = "",colClasses = list(character=10))
  #system("rm postogrep.txt filtered.germline.pileup.txt")
  ctrl.pileup = ctrl.pileup[,1:9]
  ctrl.pileup = unique(ctrl.pileup)
  ctrl.pileup = data.frame(ctrl.pileup)
  #names(ctrl.pileup)=c("chr","pos","ref","A","C","G","T","af","cov","dbsnp")
  names(ctrl.pileup)=c("chr","pos","ref","A","C","G","T","af","cov")
  # F1) Custom basic filters [ in germline ]
  common = merge(x = snvs,y = ctrl.pileup,by = c("chr","pos","ref"),all.x = T,suffixes = c("_case","_control"))
  toremove = which(common$cov_control < mincovgerm | common$af_control > maxafgerm )
  if(length(toremove)>0){
    putsnvs <- common[-toremove,,drop=F]
  } else {
    putsnvs <- common
  }
  if(nrow(putsnvs) > 0){
    # F2) Filters on Variant Allelic Fraction and add pbem [ in plasma/tumor ]
    # import pbem of this chrom
    tabpbem_file = list.files(pbem_dir,pattern = paste0('bperr_',chrom,'.tsv'),full.names = T) 
    tabpbem = data.frame(fread(input = tabpbem_file,stringsAsFactors = F,showProgress = F,header = F,colClasses = list(character=2,character=5)),stringsAsFactors = F)
    colnames(tabpbem) <- c("group","chr","pos","ref","dbsnp","gc","map","uniq","is_rndm","tot_coverage","total.A","total.C","total.G","total.T","n_pos_available",'n_pos_af_lth','n_pos_af_gth','count.A_af_gth','count.C_af_gth','count.G_af_gth','count.T_af_gth',"bperr","tot_reads_supporting_alt")
    # TABLE 1
    chrpmF1 = apply_AF_filters(chrpmF1=putsnvs,AFbycov=AFbycov,mybreaks=covbin,af.threshold.table=minaf_cov_corrected,minaf=minaf,mc.cores=mc.cores)
    chrpmF1 = chrpmF1[,c("chr","pos","ref","alt","A_case","C_case","G_case","T_case","af_case","cov_case","Ars","Crs","Grs","Trs","rev.ref","fwd.ref","cov.alt","rev.alt","fwd.alt","strandbias","A_control","C_control","G_control","T_control","af_control","cov_control","af_threshold")]
    chrpmF1$group <- paste(chrpmF1$chr,chrpmF1$pos,chrpmF1$ref,sep = ":")
    #cpmf1 = merge(x = chrpmF1,y = tabpbem,by = c("group","chr","pos","ref","dbsnp"),all.x = T)
    cpmf1 = merge(x = chrpmF1,y = tabpbem,by = c("group","chr","pos","ref"),all.x = T)
    cpmf1 = cpmf1[,c("group","chr","pos","ref","dbsnp","alt","A_case","C_case","G_case","T_case","af_case","cov_case","Ars","Crs","Grs","Trs","rev.ref","fwd.ref","cov.alt","rev.alt","fwd.alt","strandbias","A_control","C_control","G_control","T_control","af_control","cov_control","af_threshold","gc","map","uniq","is_rndm","tot_coverage","total.A","total.C","total.G","total.T","n_pos_available","n_pos_af_lth","n_pos_af_gth","count.A_af_gth","count.C_af_gth","count.G_af_gth","count.T_af_gth","bperr","tot_reads_supporting_alt")]
    # Add sample/patient IDs 
    cpmf1 = add_names(pm = cpmf1,name.patient = name.patient,name.plasma = name.plasma,name.germline = name.germline)
    # compute pbem allele
    Nids = which(cpmf1$alt=='N')
    if(length(Nids)>0){cpmf1 = cpmf1[-Nids,]}
    cpmf1 = cpmf1[which(!is.na(cpmf1$cov_control)),,drop=F]
    if(nrow(cpmf1)==0){
      return()
    }
    out = mclapply(seq(1,nrow(cpmf1),1),compute_pbem_allele,abemus=cpmf1,mc.cores = mc.cores)
    cpmf1 = fromListToDF(out)
    # add CLASS standard
    cpmf1 = add_class(pmtab = cpmf1)
    # TABLE 2
    cpmf1$af_threshold[which(is.na(cpmf1$af_threshold))] <- -1
    cpmf2 = cpmf1[which(cpmf1$af_case >= cpmf1$af_threshold),,drop=F]
    if(nrow(cpmf2)==0){
      return()
    }
    cpmf2$af_threshold[which(cpmf2$af_threshold == -1)] <- NA
    # TABLE 3
    cpmf3 = add_class_xbg(pmtab = cpmf2,xbg = as.numeric(tab_bg_pbem$background_pbem))
    cpmf3$bperr[which(cpmf3$bperr > 0.2)] = 0.2
    cpmf3$bperr[which(is.na(cpmf3$bperr))] = 0.2 # assign high pbem if it is NA 
    if(nrow(cpmf3)>0){
      pbem_coverage_filter = sapply(1:nrow(cpmf3), function(k) tab_cov_pbem[min(which(covs>=cpmf3$cov_case[k])),min(which(afs>=cpmf3$bperr[k]))])
      cpmf3$filter.pbem_coverage <- pbem_coverage_filter
      cpmf3$pass.filter.pbem_coverage = 0
      cpmf3$pass.filter.pbem_coverage[which(cpmf3$af_case >= cpmf3$filter.pbem_coverage)] = 1
    }
    # Return chromosome tables
    write.table(cpmf1,file = 'chrpm_f1.tsv',sep = '\t',col.names = F,row.names = F,quote = F)
    write.table(cpmf2,file = 'chrpm_f2.tsv',sep = '\t',col.names = F,row.names = F,quote = F)
    write.table(cpmf3,file = 'chrpm_f3.tsv',sep = '\t',col.names = F,row.names = F,quote = F)
    cat(paste(colnames(cpmf1),collapse='\t'),file = file.path(patient_folder,out1),sep = '\n')
    cat(paste(colnames(cpmf2),collapse='\t'),file = file.path(patient_folder,out2),sep = '\n')
    cat(paste(colnames(cpmf3),collapse='\t'),file = file.path(patient_folder,out3),sep = '\n')
  } else {
    return()
  }
}

chrom.in.parallel = 2
mc.cores = 10

# iterate on tumor and matched control
for(id in 1:nrow(TableSif)){
  this = TableSif[id,]
  name.patient = TableSif$patient[id]
  name.plasma = gsub(basename(this$plasma.bam),pattern = ".bam",replacement = "")
  name.germline = gsub(basename(this$germline.bam),pattern = ".bam",replacement = "")
  cat(paste("[",Sys.time(),"]\tPatient:",name.patient,"\tCase:",name.plasma,"\tControl:",name.germline),"\n")
  out1 = paste0("pmtab_F1_",name.plasma,".tsv")
  out2 = paste0("pmtab_F2_",name.plasma,".tsv")
  out3 = paste0("pmtab_F3_",name.plasma,".tsv")
  # create patient sub-folder
  patient_folder = file.path(outdir,results_folder_name,name.patient)
  dir.create(patient_folder)
  germline.folder = list.files(pacbamfolder,pattern = paste0(name.germline,"$"),full.names = T)
  if(length(germline.folder)==0){
    message("[ERROR] Cannot find folder:\t",file.path(pacbamfolder,name.germline))
    quit()
  }
  plasma.folder = list.files(pacbamfolder,pattern = paste0(name.plasma,"$"),full.names = T)
  if(length(plasma.folder)==0){
    message("[ERROR] Cannot find folder:\t",file.path(pacbamfolder,name.plasma))
    quit()
  }
  # run in parallel on chromosomes
  mclapply(seq(1,length(chromosomes),1),filter,chromosomes=chromosomes,patient_folder=patient_folder,plasma.folder=plasma.folder,germline.folder=germline.folder,out1=out1,out2=out2,mc.cores = chrom.in.parallel)
  # collapse all chromosome outs into a single table
  setwd(patient_folder)
  tabs_list = list.files(patient_folder,full.names = T,recursive = T,pattern = 'chrpm_f1.tsv') 
  if(length(tabs_list)>0){
    cmd = paste('cat',paste(tabs_list,collapse = ' '),'>>',out1)
    system(cmd)
  } else {
    cat(file = "f1_table_WARNING.txt","No calls found in chrpm_f1.tsv")
  }
  tabs_list = list.files(patient_folder,full.names = T,recursive = T,pattern = 'chrpm_f2.tsv') 
  if(length(tabs_list)>0){
    cmd = paste('cat',paste(tabs_list,collapse = ' '),'>>',out2)
    system(cmd)
  } else {
    cat(file = "f2_table_WARNING.txt","No calls found in chrpm_f2.tsv")
  }
  tabs_list = list.files(patient_folder,full.names = T,recursive = T,pattern = 'chrpm_f3.tsv') 
  if(length(tabs_list)>0){
    cmd = paste('cat',paste(tabs_list,collapse = ' '),'>>',out3)
    system(cmd)
  } else {
    cat(file = "f3_table_WARNING.txt","No calls found in chrpm_f3.tsv")
  }
}

proc.time()-timestart
