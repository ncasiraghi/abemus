#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=1){
  message("\nError!\nusage: Rscript pbem.R abemus_configure.R\n")
  quit()
}
source(args[1])
abemus_functions = file.path(abemusdir,"functions.R")
source(abemus_functions)

cat(paste("[",Sys.time(),"]\tComputation of per-base error model","\n"))

if(!file.exists(file.path(outdir, "BaseErrorModel"))){
  dir.create(file.path(outdir, "BaseErrorModel"), showWarnings = TRUE)
}
setwd(file.path(outdir, "BaseErrorModel"))
timestart <- proc.time()

# Import sample info file
TableSif = read.delim(filesif,as.is=T)
# remove NAs and keep only unique germline samples
TableSif = TableSif[which(!is.na(TableSif$germline.bam)),]
TableSif = unique(TableSif[,c(1,4,5)])

# get chromosomes
chromosomes = read.delim(file = file.path(targetbp,"bpcovered.tsv"),as.is=T)
chromosomes = sort(paste0("chr",chromosomes[-nrow(chromosomes),1]))

# Import table of bp covered by target
totpos = read.delim(file=file.path(targetbp,"bpcovered.tsv"),header = T,as.is=T,sep = "\t")
rownames(totpos) = paste0("chr",totpos$chromosome)
totpos = totpos[chromosomes,]
totpos$bp_covered_fraction = totpos$bp_covered/sum(totpos$bp_covered)

# parameters about bins of coverage and specificity to use to stratify allelic fractions 
covbin = seq(0,5000,by = coverage_binning)
covbin[length(covbin)] <- Inf
lev = levels(cut(1,breaks=covbin,include.lowest=TRUE))
probs = seq(0.95,1,0.001)

# Compute per-base error on targeted positions and save allelic fractions by bins of coverage

# temp:
bam.with.chr = F

for(chrom in chromosomes){
  cat(paste("[",Sys.time(),"]\tchromosome:",chrom),"\n")
  germlineset = get_germlineset(TableSif)
  tp = list.files(targetbp,pattern = paste0(chrom,"\\.RData$"),full.names = T)
  load(tp,verbose = F)
  chromTargetPositions = as.data.table(chromTargetPositions,keep.rownames = F)
  mytargets = unique(chromTargetPositions)
  targetsFILT <- mytargets
  targetsFILT$randompos <- 1
  # Final targets
  targets = unique(targetsFILT)
  targets = targets[with(targets,order(chr,pos)),]
  targets = as.data.frame(targets)
  if(bam.with.chr){
    targets$chr = paste0('chr',targets$chr)
  }
  cat(paste("[",Sys.time(),"]\tTotal positions to check in this chromosome :",nrow(targets)),"\n")
  step = 5000
  mclapply(seq(1,nrow(targets),step),pos2bperr,
           targets=targets,
           germlineset=germlineset,
           step=step,
           chrom=chrom,
           lev=lev,
           covbin=covbin,
           af_max_to_compute_thresholds=af_max_to_compute_thresholds,
           coverage_min_to_compute_thresholds=coverage_min_to_compute_thresholds,
           af_max_to_compute_pbem=af_max_to_compute_pbem,
           coverage_min_to_compute_pbem=coverage_min_to_compute_pbem,
           n_pos_af_th=n_pos_af_th,
           mc.cores=mc.cores)
  cmda = paste("cat *_pbem.table.txt >",paste0("bperr_",chrom,".tsv"))
  system(cmda)
  system("rm *_pbem.table.txt")
  cmdb = paste("cat *_afgtz.table.txt >",paste0("afgtz_",chrom,".tsv"))
  system(cmdb)
  system("rm *_afgtz.table.txt")
  cmdc = paste("cat *_afz.table.txt >",paste0("afz_",chrom,".tsv"))
  system(cmdc)
  system("rm *_afz.table.txt")
}
cmd.merge = paste("cat afgtz_chr*.tsv > afgtz.tsv") 
system(cmd.merge)
cmd.merge = paste("cat afz_chr*.tsv > afz.tsv") 
system(cmd.merge)

# save counter of afz as RData
afztab = read.delim(file = "afz.tsv",sep="\t",as.is = T,header=F)
afz = apply(afztab,2,sum)
names(afz)=lev
save(afz,file = "afz.RData",compress = T)

# overall statistics on pbems
cmd.merge = paste("cat bperr_chr*.tsv > bperr.tsv") 
system(cmd.merge)

bperr = fread(file.path(outdir, "BaseErrorModel","bperr.tsv"),stringsAsFactors = F,showProgress = F,header = F,colClasses = list(character=2,character=5))  

# summary stats for pbem across the target
bperr_summary = data.frame(as.numeric(summary(bperr$V22)))
bperr_summary = rbind(bperr_summary,sd(x = bperr$V16,na.rm = T))
rownames(bperr_summary) = c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.","std")
write.table(bperr_summary,file = file.path(outdir, "BaseErrorModel","bperr_summary.tsv"),row.names = T,col.names = F,quote = F,sep = "\t")  

# compute background pbem
bgpbem = (sum(as.numeric(bperr$V23)))/(sum(as.numeric(bperr$V10)))
mean_pbem = mean(as.numeric(bperr$V22))
tabstat = data.frame(background_pbem = bgpbem,
                     mean_pbem = mean_pbem,
                     stringsAsFactors = F)
write.table(tabstat,file = file.path(outdir, "BaseErrorModel","pbem_background.tsv"),row.names = F,col.names = T,quote = F,sep = "\t")  

proc.time()-timestart
