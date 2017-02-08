#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=3){
  message("\nERROR:\t3 arguments required")
  message("\nUSAGE:\tRscript dbSNPPER.R <00-All.vcf.gz> <BED> <outdir> \n")
  cat(paste("\t[1] dbsnp\tThe dbsnp file 00-All.vcf.gz",
            "\t[2] bed\t\tTarget regions in BED format",
            "\t[3] outdir\tOutput folder",sep="\n"),"\n\n")  
  quit()
}

ALL = args[1]
BED = args[2]
outfolder = args[3]

setwd(outfolder)
timestart <- proc.time()

message('DBSNP = ',ALL)
message('BED = ',BED)
message('OUTDIR = ',outfolder)

# check if 00-All.vcf.gz.tbi exists and create it if missing
tbi = paste0(ALL,'.tbi')
if(!file.exists(tbi)){
  message("Creating tbi file for ",ALL)
  cmd = paste('tabix -p vcf',ALL)
  system(cmd)
}

b = read.delim(file = BED,as.is = T,header = F,stringsAsFactors = F)
chromosomes = unique(sort(b[,1]))
chromosomes = gsub(chromosomes,pattern = 'chr',replacement = '')

bychromfolder = file.path(outfolder,paste0(basename(ALL),'.bychrom'))
if(!dir.exists(bychromfolder)){
  dir.create(bychromfolder)
}
setwd(bychromfolder)

for(chr in chromosomes){
  message('chrom:\t',chr)
  outname = paste0('dbsnp_chr',chr,'.vcf')
  cmd = paste('tabix',ALL,chr,'>',outname)
  system(cmd)
  # Intersect chrom dbsnp and bed
  out1 = gsub(basename(outname),pattern = "\\.vcf",replacement = ".intersectbed.vcf")
  cmd = paste("bedtools intersect -a",outname,"-b",BED,"-header -wa -u >",out1)
  system(cmd)
  # Keep only SNVs
  out2 = gsub(basename(out1),pattern = "\\.vcf",replacement = ".snv.vcf")
  cmd = paste("grep '\\;VC\\=SNV\\;'",out1,">",out2)
  system(cmd)
  # Keep only SNVs with a single alternative allele
  out3 = gsub(basename(out2),pattern = "\\.vcf",replacement = ".single.vcf")
  cmd = paste("awk '$5 !~ /,/{print($0)}'",out2,">",out3)
  system(cmd)
  # remove intersteps
  cmd = paste('rm',outname,out1,out2)
  system(cmd)
}

# Merge chrom vcfs into a single one 
setwd(outfolder)
vcfs = list.files(path = bychromfolder,pattern = '.intersectbed.snv.single.vcf$',full.names = T)
outfinal = gsub(basename(ALL),pattern = '.vcf.gz',replacement = '.intersectbed.snv.single.vcf')
cmd = paste('cat',paste(vcfs,collapse = ' '),'>',outfinal)
system(cmd)

proc.time()-timestart


