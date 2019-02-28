setwd("/elaborazioni/sharedCO/Abemus_data_analysis/Analysis_InSilicoDilution_STM2014/")

library(data.table)
library(ggplot2)
library(reshape2)

af_case = 0.0
cov_case = 10
reported_snvs = NULL
Rset <- seq(from = 0.5,to = 1.5,by = 0.1)
pbemAFpar <- TRUE
pars=data.frame(af_case=af_case,cov_case=cov_case,pbemAFpar=pbemAFpar,stringsAsFactors = F)

# Samples
sif = fread("/elaborazioni/sharedCO/Abemus_data_analysis/Analysis_InSilicoDilution_STM2014/SIF_dilution.tsv",data.table=F)

load("/elaborazioni/sharedCO/Home_casiraghi/Prog/abemus/data/PBEsim.RData")
tab = tab.list[[6]]
bed = fread("/elaborazioni/sharedCO/Abemus_data_analysis/Analysis_MultipleSamples_STM2014/IAD44450_30_Designed.bed",data.table=F)
patients = sapply(sif$patient,function(x) strsplit(x,"\\_")[[1]][1])
pacbam=list.files("/scratch/sharedCO/Casiraghi/InSilicoDilution/PaCBAM/",full.names = T,pattern = "\\.pileup$")

# List of reported SNVs in stm2014
stm=read.delim("/elaborazioni/sharedCO/Abemus_data_analysis/Analysis_MultipleSamples_STM2014/tableS5_stm2014_mutations.tsv",header = T,stringsAsFactors = F,as.is = T)
stm$sign=paste(stm$chr,stm$pos,sep = ":")
reported_snvs=unique(stm$sign)

if(pbemAFpar){
  # a. using pbem computed with only AFs < 0.2
  results_folder_abemus <- "/elaborazioni/sharedCO/Abemus_data_analysis/Analysis_InSilicoDilution_STM2014/ABEMUS/Results_pbemAFpar/"
  bperr = fread("/elaborazioni/sharedCO/Abemus_data_analysis/Analysis_MultipleSamples_STM2014/ABEMUS/BaseErrorModel/bperr.tsv",data.table=F)
  setwd("/elaborazioni/sharedCO/Abemus_data_analysis/Analysis_InSilicoDilution_STM2014/comparison_abemus_strelka_ssniper_pbemAFpar")
} else {
  # b. using pbem computed with all AFs
  results_folder_abemus <- "/elaborazioni/sharedCO/Abemus_data_analysis/Analysis_InSilicoDilution_STM2014/ABEMUS/Results_pbem_no_AFpar_cftrlAF_0.2/"
  #bperr = fread("/scratch/sharedCO/Casiraghi/Abemus_data/AmpliSeq_AR_STM2015/BaseErrorModel_113/bperr.tsv",data.table=F)
  bperr = fread("/elaborazioni/sharedCO/Abemus_data_analysis/Analysis_InSilicoDilution_STM2014/ABEMUS/BaseErrorModel_pbem_no_AFpar/bperr.tsv",data.table=F)
  setwd("/elaborazioni/sharedCO/Abemus_data_analysis/Analysis_InSilicoDilution_STM2014/comparison_abemus_strelka_ssniper_pbem_no_AFpar")
}

results_folder_strelka <- "/elaborazioni/sharedCO/Abemus_data_analysis/Analysis_InSilicoDilution_STM2014/Strelka2"
results_folder_ssniper<- "/scratch/sharedCO/Casiraghi/InSilicoDilution/SomaticSniper"

write.table(pars,file="pars_for_filtering.tsv",col.names=T,row.names=F,quote=F,sep="\t")

calls_abemus_strelka <- function(tmp,s,af_case,cov_case,R,reported_snvs=NULL){
  # STRELKA2
  strelka_file = file.path(results_folder_strelka,gsub(".bam","",basename(tmp$plasma.bam[s])),"results/variants/somatic.snvs.vcf")
  strelka = read.delim(file = strelka_file,comment.char = "#",stringsAsFactors = F,header = F,as.is = T)
  strelka = strelka[which(strelka[,7]=="PASS"),,drop=F]
  strelka$sign = paste(strelka[,1],strelka[,2],sep=":")
  a=fread(grep(pacbam,pattern = gsub(".bam",".pileup",basename(tmp$plasma.bam[s])),value = T),data.table = F)
  a$sign=paste(a$chr,a$pos,sep = ":")
  strelka=merge(x = a,y = strelka,by = "sign",all.y = T)
  strelka=strelka[which(strelka$af>=af_case),]
  strelka=strelka[which(strelka$cov>=cov_case),]
  #check=check[which(check$af>=af_case),]
  #check=check[which(check$cov>=cov_case),]
  #strelka=strelka[which(strelka$sign%in%check$sign),,drop=F]
  if(!is.null(reported_snvs)){
    strelka = strelka[which(strelka$sign %in% reported_snvs),]
  }
  # ABEMUS
  abemus_file = file.path(results_folder_abemus,paste(tmp$patient[s],"/pmtab_F3_",gsub(".bam","",basename(tmp$plasma.bam[s])),".tsv",sep=""))
  if(!file.exists(abemus_file)){
    abemus=data.frame(sign=character(),stringsAsFactors = F)
    return(list(abemus,strelka))
  } else {
    abemus = fread(abemus_file,data.table=F)
    abemus = abemus[which(abemus$af_case>=af_case),,drop=F]
    abemus = abemus[which(abemus$cov_case>=cov_case),,drop=F]
    if(nrow(abemus)==0){
      abemus=data.frame(sign=character(),stringsAsFactors = F)
      return(list(abemus,strelka))
    }
    # update scaling factor R
    abemus$filter.pbem_coverage=abemus$filter.pbem_coverage*R
    abemus$pass.filter.pbem_coverage=0
    abemus$pass.filter.pbem_coverage[which(abemus$af_case>=abemus$filter.pbem_coverage)]=1
    abemus = abemus[which(abemus$pass.filter.pbem_coverage==1),]
    if(nrow(abemus)==0){
      abemus=data.frame(sign=character(),stringsAsFactors = F)
      return(list(abemus,strelka))
    }
    abemus$sign = paste(abemus[,2],abemus[,3],sep=":")
    if(!is.null(reported_snvs)){
      abemus = abemus[which(abemus$sign %in% reported_snvs),]
    }
    return(list(abemus,strelka))
  }
  
}

calls_abemus_sniper <- function(tmp,s,af_case,cov_case,R,reported_snvs=NULL){
  # SOMATIC SNIPER
  sniper_file = file.path(results_folder_ssniper,gsub(".bam|.sorted.bam",".pm",basename(tmp$plasma.bam[s])))
  if(!file.exists(sniper_file)){
    sniper=data.frame(sign=character(),stringsAsFactors = F)
  } else {
    sniper = read.delim(file = sniper_file,comment.char = "#",stringsAsFactors = F,header = F,as.is = T)
    sniper$sign = paste(sniper[,1],sniper[,2],sep=":")
    a=fread(grep(pacbam,pattern = gsub(".bam",".pileup",basename(tmp$plasma.bam[s])),value = T),data.table = F)
    a$sign=paste(a$chr,a$pos,sep = ":")
    sniper=merge(x = a,y = sniper,by = "sign",all.y = T)
    sniper=sniper[which(sniper$af>=af_case),]
    sniper=sniper[which(sniper$cov>=cov_case),]
    #check=check[which(check$af>=af_case),]
    #check=check[which(check$cov>=cov_case),]
    #strelka=strelka[which(strelka$sign%in%check$sign),,drop=F]
    if(!is.null(reported_snvs)){
      sniper = sniper[which(sniper$sign %in% reported_snvs),]
    }
  }

  # ABEMUS
  abemus_file = file.path(results_folder_abemus,paste(tmp$patient[s],"/pmtab_F3_",gsub(".bam","",basename(tmp$plasma.bam[s])),".tsv",sep=""))
  if(!file.exists(abemus_file)){
    abemus=data.frame(sign=character(),stringsAsFactors = F)
    return(list(abemus,sniper))
  } else {
    abemus = fread(abemus_file,data.table=F)
    abemus = abemus[which(abemus$af_case>=af_case),,drop=F]
    abemus = abemus[which(abemus$cov_case>=cov_case),,drop=F]
    if(nrow(abemus)==0){
      abemus=data.frame(sign=character(),stringsAsFactors = F)
      return(list(abemus,sniper))
    }
    # update scaling factor R
    abemus$filter.pbem_coverage=abemus$filter.pbem_coverage*R
    abemus$pass.filter.pbem_coverage=0
    abemus$pass.filter.pbem_coverage[which(abemus$af_case>=abemus$filter.pbem_coverage)]=1
    abemus = abemus[which(abemus$pass.filter.pbem_coverage==1),]
    if(nrow(abemus)==0){
      abemus=data.frame(sign=character(),stringsAsFactors = F)
      return(list(abemus,sniper))
    }
    abemus$sign = paste(abemus[,2],abemus[,3],sep=":")
    if(!is.null(reported_snvs)){
      abemus = abemus[which(abemus$sign %in% reported_snvs),]
    }
    return(list(abemus,sniper))
  }
  
}


# Generate count tables for comparison

# ABEMUS vs Strelka
for(R in Rset){

  count_tab = c()
  pbem_only_abemus = list()
  pbem_only_strelka= list()
  pbem_intersect = list()
  af_only_abemus = list()
  af_only_strelka= list()
  af_intersect = list()

  for(p in unique(patients)){
    message(p,"\t[strelka]")
    tmp = sif[grep(sif$patient,pattern = p),]
    tmp = tmp[with(tmp,order(patient,decreasing = F)),]
    abemus.mut = c()
    strelka.mut = c()
    for(s in 1:nrow(tmp)){
      out=calls_abemus_strelka(tmp = tmp,s = s,af_case = af_case,cov_case = cov_case,R = R,reported_snvs = reported_snvs)
      abemus=out[[1]]
      strelka=out[[2]]
      # Start counts
      stm_pm_names=paste(stm$sign[which(stm$Sample==substring(p,4))],sep = "|")
      stm_pm_abemus=rep(0,length(stm_pm_names))
      stm_pm_strelka=rep(0,length(stm_pm_names))
      stm_pm_abemus[which(stm_pm_names%in%abemus$sign)]<-1
      stm_pm_strelka[which(stm_pm_names%in%strelka$sign)]<-1
      count_tab=rbind(count_tab,data.frame(patient=p,
                                           case=tmp$plasma[s],
                                           abemus=length(unique(abemus$sign)),
                                           strelka=length(unique(strelka$sign)),
                                           only_abemus=length(setdiff(unique(abemus$sign),unique(strelka$sign))),
                                           only_strelka=length(setdiff(unique(strelka$sign),unique(abemus$sign))),
                                           intersect_abemus_strelka=length(intersect(unique(abemus$sign),unique(strelka$sign))),
                                           stm_pm=paste(stm_pm_names,collapse= "|"),
                                           stm_pm_abemus=paste(stm_pm_abemus,collapse = "|"),
                                           stm_pm_strelka=paste(stm_pm_strelka,collapse = "|"),
                                           stm_pm_abemus_called=sum(stm_pm_abemus),
                                           stm_pm_strelka_called=sum(stm_pm_strelka),
                                           stringsAsFactors = F))
      pbem_only_abemus[[tmp$plasma[s]]] = abemus$bperr[which(abemus$sign%in%setdiff(unique(abemus$sign),unique(strelka$sign)))]
      pbem_intersect[[tmp$plasma[s]]] = abemus$bperr[which(abemus$sign%in%intersect(unique(abemus$sign),unique(strelka$sign)))]
      strelka_tmp <- strelka
      strelka_tmp$group = paste(strelka_tmp$sign,strelka_tmp$V4,sep = ":")
      strelka_tmp = merge(x = strelka_tmp,y = bperr,by.x = "group",by.y = "V1",all.x = T)
      pbem_only_strelka[[tmp$plasma[s]]]= strelka_tmp$V22[which(strelka_tmp$sign%in%setdiff(unique(strelka_tmp$sign),unique(abemus$sign)))]
      rm(strelka_tmp)
      af_only_abemus[[tmp$plasma[s]]] = abemus$af_case[which(abemus$sign%in%setdiff(unique(abemus$sign),unique(strelka$sign)))]
      af_only_strelka[[tmp$plasma[s]]] = strelka$af[which(strelka$sign%in%setdiff(unique(strelka$sign),unique(abemus$sign)))]
      af_intersect[[tmp$plasma[s]]] = abemus$af_case[which(abemus$sign%in%intersect(unique(abemus$sign),unique(strelka$sign)))]
      # End counts
    }
  }

  save(count_tab,pbem_only_abemus,pbem_only_strelka,pbem_intersect,af_only_abemus,af_only_strelka,af_intersect,file = paste0("count_tab",R,".strelka.RData"))

}

# ABEMUS vs Somatic Sniper
for(R in Rset){
  
  count_tab = c()
  pbem_only_abemus = list()
  pbem_only_sniper= list()
  pbem_intersect = list()
  af_only_abemus = list()
  af_only_sniper= list()
  af_intersect = list()
  
  for(p in unique(patients)){
    message(p,"\t[sniper]")
    tmp = sif[grep(sif$patient,pattern = p),]
    tmp = tmp[with(tmp,order(patient,decreasing = F)),]
    abemus.mut = c()
    sniper.mut = c()
    for(s in 1:nrow(tmp)){
      out=calls_abemus_sniper(tmp = tmp,s = s,af_case = af_case,cov_case = cov_case,R = R,reported_snvs = reported_snvs)
      abemus=out[[1]]
      sniper=out[[2]]
      # Start counts
      stm_pm_names=paste(stm$sign[which(stm$Sample==substring(p,4))],sep = "|")
      stm_pm_abemus=rep(0,length(stm_pm_names))
      stm_pm_sniper=rep(0,length(stm_pm_names))
      stm_pm_abemus[which(stm_pm_names%in%abemus$sign)]<-1
      stm_pm_sniper[which(stm_pm_names%in%sniper$sign)]<-1
      count_tab=rbind(count_tab,data.frame(patient=p,
                                           case=tmp$plasma[s],
                                           abemus=length(unique(abemus$sign)),
                                           sniper=length(unique(sniper$sign)),
                                           only_abemus=length(setdiff(unique(abemus$sign),unique(sniper$sign))),
                                           only_sniper=length(setdiff(unique(sniper$sign),unique(abemus$sign))),
                                           intersect_abemus_sniper=length(intersect(unique(abemus$sign),unique(sniper$sign))),
                                           stm_pm=paste(stm_pm_names,collapse= "|"),
                                           stm_pm_abemus=paste(stm_pm_abemus,collapse = "|"),
                                           stm_pm_sniper=paste(stm_pm_sniper,collapse = "|"),
                                           stm_pm_abemus_called=sum(stm_pm_abemus),
                                           stm_pm_sniper_called=sum(stm_pm_sniper),
                                           stringsAsFactors = F))
      pbem_only_abemus[[tmp$plasma[s]]] = abemus$bperr[which(abemus$sign%in%setdiff(unique(abemus$sign),unique(sniper$sign)))]
      pbem_intersect[[tmp$plasma[s]]] = abemus$bperr[which(abemus$sign%in%intersect(unique(abemus$sign),unique(sniper$sign)))]
      sniper_tmp <- sniper 
      sniper_tmp$group = paste(sniper_tmp$sign,sniper_tmp$V4,sep = ":")
      sniper_tmp = merge(x = sniper_tmp,y = bperr,by.x = "group",by.y = "V1",all.x = T)
      pbem_only_sniper[[tmp$plasma[s]]]= sniper_tmp$V22[which(sniper_tmp$sign%in%setdiff(unique(sniper_tmp$sign),unique(abemus$sign)))]
      rm(sniper_tmp)
      af_only_abemus[[tmp$plasma[s]]] = abemus$af_case[which(abemus$sign%in%setdiff(unique(abemus$sign),unique(sniper$sign)))]
      af_only_sniper[[tmp$plasma[s]]] = sniper$af[which(sniper$sign%in%setdiff(unique(sniper$sign),unique(abemus$sign)))]
      af_intersect[[tmp$plasma[s]]] = abemus$af_case[which(abemus$sign%in%intersect(unique(abemus$sign),unique(sniper$sign)))]
      # End counts
    }
  }
  
  save(count_tab,pbem_only_abemus,pbem_only_sniper,pbem_intersect,af_only_abemus,af_only_sniper,af_intersect,file = paste0("count_tab",R,".sniper.RData"))
  
}
