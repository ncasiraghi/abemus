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
if(T){
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

# #  check if data tables exists
# if(!is.numeric(AFbycov)){
#   if(file.exists(file.path(controls_dir,"datathreshold_high_coverages.RData"))){
#     cat(paste("[",Sys.time(),"]\tlooking for data table with AF thresholds:",file.path(controls_dir,"datathreshold_high_coverages.RData"),"[ ok ]"),"\n")
#     load(file = file.path(controls_dir,"datathreshold_high_coverages.RData"))
#   } else {
#     cat(paste("[",Sys.time(),"]\tlooking for data table with AF thresholds:",file.path(controls_dir,"datathreshold_high_coverages.RData"),"[ NOT found ]"),"\n")
#     quit()
#   }
# }


pmtableF1 = data.frame(stringsAsFactors = F,check.names = F)
pmtableF2 = data.frame(stringsAsFactors = F,check.names = F)
ftabstats = data.frame(stringsAsFactors = F,check.names = F)

for(id in 1:nrow(TableSif)){
  this = TableSif[id,]
  name.patient = TableSif$patient[id]
  name.plasma = gsub(basename(this$plasma.bam),pattern = ".bam",replacement = "")
  name.germline = gsub(basename(this$germline.bam),pattern = ".bam",replacement = "")
  cat(paste("\n[",Sys.time(),"]\tPatient:",name.patient,"\tCase:",name.plasma,"\tControl:",name.germline),"\n")
  pmF1 = data.frame(stringsAsFactors = F,check.names = F)
  pmF2 = data.frame(stringsAsFactors = F,check.names = F)
  germline.folder = list.files(pacbamfolder,pattern = name.germline,full.names = T)
  plasma.folder = list.files(pacbamfolder,pattern = name.plasma,full.names = T)
  for(chrom in chromosomes){
    fts = data.frame(stringsAsFactors = F,check.names = F,row.names = name.patient, id=name.patient,chr=chrom,snvs=0,snvs_f1p=0,snvs_f1g=0,snvs_f2=0) 
    cat(paste("[",Sys.time(),"]\tchromosome:",chrom),"\n")
    # Upload tumor/plasma sample [ snvs ]
    plasma_snvs = list.files(file.path(plasma.folder,"snvs"),pattern = paste0("_",chrom,".snvs"),full.names = T)
    snvs = fread(plasma_snvs,stringsAsFactors = F,showProgress = F,header = F,skip = 1,na.strings = "",colClasses = list(character=15))
    snvs = unique(snvs)
    snvs = data.frame(snvs)
    names(snvs)=c("chr","pos","ref","alt","A","C","G","T","af","cov","Ars","Crs","Grs","Trs","dbsnp")
    cat(paste("[",Sys.time(),"]\tn. of snvs in case :",nrow(snvs)),"\n")
    fts[name.patient,"snvs"] <- nrow(snvs)
    if(nrow(snvs)==0){
      next
    }
    # F1) Custom basic filters [ in plasma/tumor ]
    out = mclapply(seq(1,nrow(snvs),1),CheckAltReads,snvs=snvs,mc.cores = mc.cores)
    snvs = fromListToDF(out)
    snvs = snvs[which(snvs$af > 0 ),]
    snvs = snvs[which(snvs$cov > mincov ),]
    snvs = snvs[which(snvs$cov.alt >= minalt ),]
    snvs = unique(snvs)
    cat(paste("[",Sys.time(),"]\tn. of snvs in case after custom basic filters | F1 applied in case:\t",nrow(snvs)),"\n")
    fts[name.patient,"snvs_f1p"] <- nrow(snvs)
    if(nrow(snvs)==0){
      next
    }
    # print filtered positions and grep these pos only from pileup file of germline sample
    cat(unique(snvs$pos),sep = "\n",file = "postogrep.txt",append = F)
    controlfolder_pileup <- list.files(file.path(germline.folder,"pileup"),pattern = paste0("_",chrom,".pileup"),full.names = T)
    cmd = paste("awk -F'\t' '{if (FILENAME == \"postogrep.txt\") { t[$1] = 1; } else { if (t[$2]) { print }}}' postogrep.txt",controlfolder_pileup,"> filtered.germline.pileup.txt")
    system(cmd)
    ctrl.pileup = fread("filtered.germline.pileup.txt",stringsAsFactors = F,showProgress = T,header = F,na.strings = "",colClasses = list(character=10))
    system("rm postogrep.txt filtered.germline.pileup.txt")
    ctrl.pileup = unique(ctrl.pileup)
    ctrl.pileup = as.data.frame(ctrl.pileup)
    names(ctrl.pileup)=c("chr","pos","ref","A","C","G","T","af","cov","dbsnp")   
    # F1) Custom basic filters [ in germline ]
    common = merge(x = snvs,y = ctrl.pileup,by = c("chr","pos","ref","dbsnp"),all.x = T,suffixes = c("_case","_control"))
    toremove = which(common$cov_control < mincovgerm | common$af_control > maxafgerm )
    if(length(toremove)>0){
      putsnvs <- common[-toremove,,drop=F]
    } else {
      putsnvs <- common
    }
    cat(paste("[",Sys.time(),"]\tn. of snvs in case after custom basic filters | F1 applied in ctrl:\t",nrow(putsnvs)),"\n")
    fts[name.patient,"snvs_f1g"] <- nrow(putsnvs)
    if(nrow(putsnvs) > 0){
      # F2) Filters on Variant Allelic Fraction [ in plasma/tumor ]
      filtafout = apply_AF_filters(chrpmF1=putsnvs,AFbycov=AFbycov,mybreaks=covbin,minaf_cov=minaf_cov,minaf=minaf)
      chrpmF1 = filtafout[[1]]
      chrpmF2 = filtafout[[2]]
      cat(paste("[",Sys.time(),"]\tn. of snvs in case after custom basic filters | F2 applied in case:\t",nrow(chrpmF2)),"\n")
      fts[name.patient,"snvs_f2"] <- nrow(chrpmF2)
      # Update table
      pmF1 <- rbind(pmF1,chrpmF1)
      pmF2 <- rbind(pmF2,chrpmF2)
    }
    ftabstats <- rbind(ftabstats,fts)
  }
  snvs_in_this = ftabstats[which(ftabstats$id==name.patient),]
  total = data.frame(id=name.patient,chr="total",snvs=sum(snvs_in_this$snvs),snvs_f1p=sum(snvs_in_this$snvs_f1p),snvs_f1g=sum(snvs_in_this$snvs_f1g),snvs_f2=sum(snvs_in_this$snvs_f2),stringsAsFactors = F)
  ftabstats = rbind(ftabstats,total)
  # Add IDs
  if(nrow(pmF1)>0){
    pmF1$PatientID = name.patient
    pmF1$CaseID = name.plasma
    pmF1$GermlineID = name.germline
    pmF1 = pmF1[with(pmF1,order(chr,pos)),]
    # write.table(pmF1,file = paste0("pmtab_F1_",name.plasma,".tsv"),quote=F,sep="\t",col.names = T,row.names = F)
    pmtableF1 = rbind(pmtableF1,pmF1)
  }
  if(nrow(pmF2)>0){
    pmF2$PatientID = name.patient
    pmF2$CaseID = name.plasma
    pmF2$GermlineID = name.germline
    pmF2 = pmF2[with(pmF2,order(chr,pos)),]
    # write.table(pmF2,file = paste0("pmtab_F2_",name.plasma,".tsv"),quote=F,sep="\t",col.names = T,row.names = F)
    pmtableF2 = rbind(pmtableF2,pmF2)
  }
}

cat(paste("\n[",Sys.time(),"]\tSorting output tables"),"\n")
pmtableF1 = pmtableF1[with(pmtableF1,order(chr,pos,PatientID)),]
pmtableF2 = pmtableF2[with(pmtableF2,order(chr,pos,PatientID)),]

cat(paste("\n[",Sys.time(),"]\tAdd per-base error measure"),"\n")
tabpbem = data.frame(fread(file.path(pbem_dir,"bperr.tsv"),stringsAsFactors = F,showProgress = F,header = F,colClasses = list(character=2,character=5)),stringsAsFactors = F)  
colnames(tabpbem) <- c("group","chr","pos","ref","dbsnp","gc","map","uniq","is_rndm","tot_coverage","total.A","total.C","total.G","total.T","n_pos_available",'n_pos_af_lth','n_pos_af_gth','count.A_af_gth','count.C_af_gth','count.G_af_gth','count.T_af_gth',"bperr","tot_reads_supporting_alt")

pmtableF1$group <- paste(pmtableF1$chr,pmtableF1$pos,pmtableF1$ref,sep = ":")
pmtableF2$group <- paste(pmtableF2$chr,pmtableF2$pos,pmtableF2$ref,sep = ":")

pmtableF13 = merge(x = pmtableF1,y = tabpbem,by = c("group","chr","pos","ref","dbsnp"),all.x = T)
pmtableF23 = merge(x = pmtableF2,y = tabpbem,by = c("group","chr","pos","ref","dbsnp"),all.x = T)

pmtableF13_final = pmtableF13[with(pmtableF13,order(chr,pos,PatientID)),]
pmtableF23_final = pmtableF23[with(pmtableF23,order(chr,pos,PatientID)),]

cat(paste("[",Sys.time(),"]\tWriting output tables"),"\n")

# By pairs [ table F1 + pbem ]
for(sname in unique(pmtableF13_final$CaseID)){
  tabout_temp = pmtableF13_final[which(pmtableF13_final$CaseID==sname),]
  
  tabout_temp = tabout_temp[which(tabout_temp$alt!='N'),]
  out = mclapply(seq(1,nrow(tabout_temp),1),compute_pbem_allele,abemus=tabout_temp,mc.cores = mc.cores)
  tabout = fromListToDF(out)
  
  write.table(x = tabout,file = paste0('pmtab_F1_pbem_',sname,'.tsv'),quote=F,sep="\t",col.names = T,row.names = F)
}

# By pairs [ table F2 + pbem ]
for(sname in unique(pmtableF23_final$CaseID)){
  tabout = pmtableF23_final[which(pmtableF23_final$CaseID==sname),]
  write.table(x = tabout,file = paste0('pmtab_F2_pbem_',sname,'.tsv'),quote=F,sep="\t",col.names = T,row.names = F)
}

# All samples together
write.table(pmtableF13_final,file = "pmtab_F1_pbem_AllSamples.tsv",quote=F,sep="\t",col.names = T,row.names = F)
write.table(pmtableF23_final,file = "pmtab_F2_pbem_AllSamples.tsv",quote=F,sep="\t",col.names = T,row.names = F)
write.table(ftabstats,file = "ftabstats.txt",quote=F,sep="\t",col.names = T,row.names = F)

if(F){
cat(paste("[",Sys.time(),"]\tConverting output tables into VCF format"),"\n")
# pmtable 1
pmtableF1.vcf = as.data.frame(cbind(pmtableF1$chr,pmtableF1$pos,".",pmtableF1$ref,pmtableF1$alt,".",".","."),stringsAsFactors = F)
names(pmtableF1.vcf)<-c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
pmtableF1.vcf = pmtableF1.vcf[which(pmtableF1.vcf$ALT != "N"),]
pmtableF1.vcf = unique(pmtableF1.vcf)
# pmtable 2
pmtableF2.vcf = as.data.frame(cbind(pmtableF2$chr,pmtableF2$pos,".",pmtableF2$ref,pmtableF2$alt,".",".","."),stringsAsFactors = F)
names(pmtableF2.vcf)<-c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
pmtableF2.vcf = pmtableF2.vcf[which(pmtableF2.vcf$ALT != "N"),]
pmtableF2.vcf = unique(pmtableF2.vcf)
# pmtable 3
pmtableF3.vcf = as.data.frame(cbind(pmtableF3$chr,pmtableF3$pos,".",pmtableF3$ref,pmtableF3$alt,".",".","."),stringsAsFactors = F)
names(pmtableF3.vcf)<-c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
pmtableF3.vcf = pmtableF3.vcf[which(pmtableF3.vcf$ALT != "N"),]
pmtableF3.vcf = unique(pmtableF3.vcf)

cat(paste("[",Sys.time(),"]\tWriting VCF output tables"),"\n")
write.table(pmtableF1.vcf,file = "pmtab_F1_AllSamples.vcf",quote=F,sep="\t",col.names = T,row.names = F)
write.table(pmtableF2.vcf,file = "pmtab_F2_AllSamples.vcf",quote=F,sep="\t",col.names = T,row.names = F)
write.table(pmtableF3.vcf,file = "pmtab_F3_AllSamples.vcf",quote=F,sep="\t",col.names = T,row.names = F)
}

}

# Use Annovar to add functional annotations to retained putative SNVs and save data
if(F){
cat(paste("\n[",Sys.time(),"]\tuse Annovar to add functional annotations to retained putative SNVs and save data"),"\n")  
convert2annovar = paste("/elaborazioni/sharedCO/CO_Shares/Code/annovar/convert2annovar.pl --format vcf4old --outfile pmtab_F2_AllSamples.avinput pmtab_F2_AllSamples.vcf")
system(convert2annovar)
table_annovar = paste("/elaborazioni/sharedCO/CO_Shares/Code/annovar/table_annovar.pl pmtab_F2_AllSamples.avinput /elaborazioni/sharedCO/CO_Shares/Code/annovar/humandb/ -remove -protocol refGene,esp6500siv2_all,exac03,1000g2015aug_all,avsnp144,snp138,snp138NonFlagged,kaviar_20150923,clinvar_20150629,ljb26_all,cosmic70 -operation g,f,f,f,f,f,f,f,f,f,f --buildver hg19 --argument '--splicing_threshold 5 --hgvs --exonicsplicing','','','','','','','','','','' --outfile pmtab_F2_AllSamples.avinput.refGene")
system(table_annovar)
}

# Use SnpEff to add functional annotations to retained putative SNVs and save data
if(F){
cat(paste("\n[",Sys.time(),"]\tuse SnpEff to add functional annotations to retained putative SNVs and save data"),"\n")

snpSift  <- 'java -Xmx16g -jar /elaborazioni/sharedCO/Home_casiraghi/Prog/snpEff_01142015/snpEff/SnpSift.jar annotate -id'
snpEff   <- 'java -Xmx16g -jar /elaborazioni/sharedCO/Home_casiraghi/Prog/snpEff_01142015/snpEff/snpEff.jar -c /elaborazioni/sharedCO/Home_casiraghi/Prog/snpEff_01142015/snpEff/snpEff.config'
database <- 'hg19'
dbsnp    <- '/elaborazioni/sharedCO/CO_Shares/dbSNP/dbSNP144/SNV/Haloplex/All.onlySNV.Header.vcf'  

input <- list.files(path = file.path(outdir, "Results"),pattern = "_F2_AllSamples\\.vcf$")
cat("\nVCF to annotate:",input,"\n")

# Run SnpSift
cat('\nSnpSift annotation: add IDs from dbsnp ...\t')
file.dbsnp.vcf <- gsub(basename(input),pattern = ".vcf",replacement = ".dbsnp.vcf")
cmd = paste(snpSift,dbsnp,input,'>',file.dbsnp.vcf,sep = ' ')
system(cmd)
cat('[ OK ]\n')

# run SnpEff
cat('SnpEff  annotation: add EFF field ...\t')
input  <- file.dbsnp.vcf
output <- gsub(basename(file.dbsnp.vcf),pattern = ".vcf",replacement = ".eff.vcf")
cmd = paste(snpEff,database,input,'>',output,sep = ' ')
system(cmd)
cat('[ OK ]\n')

# SnpSift one variant per line
cat('Cleaning snpEff file output ...\t')
output.clean = gsub(output,pattern = '.vcf',replacement = '.all.tmp')

cmd = paste('cat',output,'| /elaborazioni/sharedCO/Home_casiraghi/Prog/snpEff_01142015/snpEff/scripts/vcfEffOnePerLine.pl |',
            'java -jar /elaborazioni/sharedCO/Home_casiraghi/Prog/snpEff_01142015/snpEff/SnpSift.jar extractFields -s "|" -e "NA" - CHROM POS ID REF ALT FILTER "ANN[*].ALLELE" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].HGVS_P" "ANN[*].EFFECT" "ANN[*].IMPACT"',
            '>',output.clean)
system(cmd)

# Remove possible duplicates from annotated table
all = read.table(file=output.clean,stringsAsFactors = F)
output.all = gsub(output,pattern = '.vcf',replacement = '.all.vcf')
nodup = unique(all)
colnames(nodup) = c("#CHROM","POS","ID","REF","ALT","FILTER","ALLELE","GENE","GENEID","HGVS_P","EFFECT","IMPACT")
write.table(nodup,file=output.all,sep="\t",col.names = T,row.names=F,quote=F)
system("rm snpEff_genes.txt snpEff_summary.html *.dbsnp.vcf *.dbsnp.eff.vcf *.all.tmp")
cat('[ OK ]\n')

}

proc.time()-timestart
