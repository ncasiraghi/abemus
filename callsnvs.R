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
  dir.create(file.path(outdir, "Results"), showWarnings = TRUE)
}
setwd(file.path(outdir, "Results"))
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
                  spec = as.character(spec),
                  mincov = as.character(mincov),
                  minalt = as.character(minalt),
                  mincovgerm = as.character(mincovgerm),
                  maxafgerm = as.character(maxafgerm),
                  filtering_date = Sys.time(),
                  stringsAsFactors = F)
write.table(fpam,file = "filtering_criteria.tsv",col.names = T,row.names = F,sep="\t",quote = F)

#  check if data tables exists
if(!is.numeric(AFbycov)){
  if(file.exists(file.path(controls_dir,"datathreshold.RData"))){
    cat(paste("[",Sys.time(),"]\tlooking for data table with AF thresholds:",file.path(controls_dir,"datathreshold.RData"),"[ ok ]"),"\n")
    load(file = file.path(controls_dir,"datathreshold.RData"))
  } else {
    cat(paste("[",Sys.time(),"]\tlooking for data table with AF thresholds:",file.path(controls_dir,"datathreshold.RData"),"[ NOT found ]"),"\n")
    quit()
  }
}

add_names = function(pm,name.patient,name.plasma,name.germline){
  pm$PatientID = name.patient
  pm$CaseID = name.plasma
  pm$GermlineID = name.germline
  pm = pm[with(pm,order(chr,pos)),]
  return(pm)
}

add_class = function(pmtab){
  pmtab$CLASS = NA
  pmtab$CLASS[which(pmtab$af_control==0 & pmtab$bperr==0 & pmtab$pbem_allele==0)] = 1
  pmtab$CLASS[which(pmtab$af_control==0 & pmtab$bperr>0 & pmtab$pbem_allele==0)] = 2
  pmtab$CLASS[which(pmtab$af_control==0 & pmtab$bperr>0 & pmtab$pbem_allele>0)] = 3
  pmtab$CLASS[which(pmtab$af_control>0 & pmtab$bperr>0 & pmtab$pbem_allele>=0 & pmtab$same_allele == 0)] = 4
  pmtab$CLASS[which(pmtab$af_control>0 & pmtab$bperr>0 & pmtab$pbem_allele>0 & pmtab$same_allele == 1)] = 5
  return(pmtab)
}

apply_AF_filters <- function(chrpmF1,AFbycov,minaf_cov,minaf,mybreaks,mc.cores){
  if (AFbycov == FALSE){
    chrpmF1[,'af_threshold'] <- minaf
    #chrpmF2 = chrpmF1[which(chrpmF1$af_case >= minaf),,drop=F]
  } else if (is.numeric(AFbycov)){
    chrpmF1[,'af_threshold'] <- AFbycov
    #chrpmF2 = chrpmF1[which(chrpmF1$af_case >= AFbycov),,drop=F]
  } else if (AFbycov == TRUE){
    thresholds = as.numeric(minaf_cov[,-1])
    
    af_filter_by_coverage <- function(y,thresholds,chrpmF1){
      this = chrpmF1[y,,drop=F]
      minaf_covth = thresholds[findInterval(this$cov_case,mybreaks)]
      this[,'af_threshold'] <- minaf_covth
      return(this)
    }
    
    out = mclapply(seq(1,nrow(chrpmF1),1),af_filter_by_coverage,thresholds=thresholds,chrpmF1=chrpmF1,mc.cores = mc.cores)
    chrpmF1 = fromListToDF(out)
    #chrpmF2 = chrpmF1[which(chrpmF1$af_case >= chrpmF1$af_threshold),,drop=F]
    
  }
  #return(list(chrpmF1,chrpmF2))
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
  snvs = fread(plasma_snvs,stringsAsFactors = F,showProgress = F,header = F,skip = 1,na.strings = "",colClasses = list(character=15),verbose = F)
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
  snvs = snvs[which(snvs$cov > mincov ),]
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
  ctrl.pileup = unique(ctrl.pileup)
  ctrl.pileup = data.frame(ctrl.pileup)
  names(ctrl.pileup)=c("chr","pos","ref","A","C","G","T","af","cov","dbsnp")   
  # F1) Custom basic filters [ in germline ]
  common = merge(x = snvs,y = ctrl.pileup,by = c("chr","pos","ref","dbsnp"),all.x = T,suffixes = c("_case","_control"))
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
    chrpmF1 = apply_AF_filters(chrpmF1=putsnvs,AFbycov=AFbycov,mybreaks=covbin,minaf_cov=minaf_cov,minaf=minaf,mc.cores=mc.cores)
    chrpmF1$group <- paste(chrpmF1$chr,chrpmF1$pos,chrpmF1$ref,sep = ":")
    cpmf1 = merge(x = chrpmF1,y = tabpbem,by = c("group","chr","pos","ref","dbsnp"),all.x = T)
    # Add sample/patient IDs 
    cpmf1 = add_names(pm = cpmf1,name.patient = name.patient,name.plasma = name.plasma,name.germline = name.germline)
    # compute pbem allele add CLASS 
    Nids = which(cpmf1$alt=='N')
    if(length(Nids)>0){cpmf1 = cpmf1[-Nids,]}
    out = mclapply(seq(1,nrow(cpmf1),1),compute_pbem_allele,abemus=cpmf1,mc.cores = mc.cores)
    cpmf1 = fromListToDF(out)
    cpmf1 = add_class(pmtab = cpmf1)
    # TABLE 2
    cpmf2 = cpmf1[which(cpmf1$af_case >= cpmf1$af_threshold),,drop=F]
    # Return chromosome tables
    write.table(cpmf1,file = 'chrpm_f1.tsv',sep = '\t',col.names = F,row.names = F,quote = F)
    write.table(cpmf2,file = 'chrpm_f2.tsv',sep = '\t',col.names = F,row.names = F,quote = F)
    cat(paste(colnames(cpmf1),collapse='\t'),file = file.path(patient_folder,out1),sep = '\n')
    cat(paste(colnames(cpmf2),collapse='\t'),file = file.path(patient_folder,out2),sep = '\n')
  } else {
    return()
  }
}

chrom.in.parallel = 5
mc.cores = 10

# iterate on tumor and matched control
for(id in 1:nrow(TableSif)){
  this = TableSif[id,]
  name.patient = TableSif$patient[id]
  name.plasma = gsub(basename(this$plasma.bam),pattern = ".bam",replacement = "")
  name.germline = gsub(basename(this$germline.bam),pattern = ".bam",replacement = "")
  cat(paste("\n[",Sys.time(),"]\tPatient:",name.patient,"\tCase:",name.plasma,"\tControl:",name.germline),"\n")
  out1 = paste0("pmtab_F1_",name.plasma,".tsv")
  out2 = paste0("pmtab_F2_",name.plasma,".tsv")
  # create patient sub-folder
  patient_folder = file.path(outdir,'Results',name.patient)
  dir.create(patient_folder)
  germline.folder = list.files(pacbamfolder,pattern = name.germline,full.names = T)
  plasma.folder = list.files(pacbamfolder,pattern = name.plasma,full.names = T)
  # run in parallel on chromosomes
  mclapply(seq(1,length(chromosomes),1),filter,chromosomes=chromosomes,patient_folder=patient_folder,plasma.folder=plasma.folder,germline.folder=germline.folder,out1=out1,out2=out2,mc.cores = chrom.in.parallel)
  # collapse all chromosome outs into a single table
  setwd(patient_folder)
  tabs_list = list.files(patient_folder,full.names = T,recursive = T,pattern = 'chrpm_f1.tsv') 
  cmd = paste('cat',paste(tabs_list,collapse = ' '),'>>',out1)
  system(cmd)
  tabs_list = list.files(patient_folder,full.names = T,recursive = T,pattern = 'chrpm_f2.tsv') 
  cmd = paste('cat',paste(tabs_list,collapse = ' '),'>>',out2)
  system(cmd)
}

proc.time()-timestart
