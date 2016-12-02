#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=4){
  message("\nError!\nUSAGE: Invalid number of input parameters\n##### USAGE dbSNPPER.R <00-All.vcf> <clinvar_YYYYMMDD.vcf> <mytarget.bed> <path/to/outfolder>\n\n")
  quit()
}

# 00-All.vcf | the dbSNP input file
ALL=args[1]
# clinvar_YYYYMMDD.vcf | the file reporting clinical rs ids
CLINVAR=args[2]
# BED target regions | BED file with regions covered by sequencing panel
BED=args[3]
# output folder
outfolder=args[4]

setwd(outfolder)
timestart <- proc.time()

message("# 1) Intersect dbSNP and BED of interest")
out1 = gsub(basename(ALL),pattern = "\\.vcf",replacement = "_intersectbed.vcf")
cmd = paste("bedtools intersect -a",ALL,"-b",BED,"-header -wa -u >",out1)
system(cmd)

message("# 2) Keep only SNVs")
out2 = gsub(basename(ALL),pattern = "\\.vcf",replacement = "_intersectbed_snv_single.vcf")
cmd = paste("head -n 500",out1,"| grep '^#' >",out2)
system(cmd)

tmp.snv = gsub(basename(ALL),pattern = "\\.vcf",replacement = "_tmp.vcf")
cmd = paste("grep '\\;VC\\=SNV\\;'",out1,">",tmp.snv)
system(cmd)

message("# 3) Keep only SNVs with a single alternative allele")
tmp.snv.single = gsub(basename(ALL),pattern = "\\.vcf",replacement = "_tmp_single.vcf")
# GetSingleALT = "/home/nicola.casiraghi/Scripts/dbSNP_adjust/get_unique_alt.py"
# cmd = paste("less",tmp.snv,"|",GetSingleALT,">",tmp.snv.single)
cmd = paste("awk '$5 !~ /,/{print($0)}'",tmp.snv,">",tmp.snv.single)
system(cmd)

cmd = paste("cat",tmp.snv.single,">>",out2)
system(cmd)

cmdrm = paste("rm",tmp.snv,tmp.snv.single)
system(cmdrm)

if(FALSE){
message("# 4) Filtering CLINVAR")
clinvar.snv = gsub(basename(CLINVAR),pattern = "\\.vcf",replacement = "_fCLNSIG.vcf")
cmd = paste("head -n 500",CLINVAR,"| grep '^#' >",clinvar.snv)
system(cmd)

tmp.snv = gsub(basename(CLINVAR),pattern = "\\.vcf",replacement = "_tmp.vcf")
cmd = paste("grep '\\;VC\\=SNV\\;'",CLINVAR,">",tmp.snv)
system(cmd)

tmp.snv.single = gsub(basename(CLINVAR),pattern = "\\.vcf",replacement = "_tmp_single.vcf")
cmd = paste("less",tmp.snv,"|",GetSingleALT,">",tmp.snv.single)
system(cmd)

tmp.snv.single.clnsig = gsub(basename(CLINVAR),pattern = "\\.vcf",replacement = "_tmp_single_clnsig.vcf")
filterCLNSIG = "/home/nicola.casiraghi/Scripts/dbSNP_adjust/filter_clinvar_by_CLNSIG.py"
cmd = paste("less",tmp.snv.single,"|",filterCLNSIG,">",tmp.snv.single.clnsig)
system(cmd)

cmd = paste("cat",tmp.snv.single.clnsig,">>",clinvar.snv)
system(cmd)

cmdrm = paste("rm",tmp.snv,tmp.snv.single,tmp.snv.single.clnsig)
system(cmdrm)

message("# 5) Extract CLINVAR somatic")
clinvar.somatic = gsub(basename(CLINVAR),pattern = "\\.vcf",replacement = "_fCLNSIG_somatic.vcf")
cmd = paste("head -n 500",CLINVAR,"| grep '^#' >",clinvar.somatic)
system(cmd)

tmp.somatic1 = gsub(basename(CLINVAR),pattern = "\\.vcf",replacement = "_tmpsomatic.vcf")
cmd = paste("grep '\\;CLNORIGIN\\=2\\;'",clinvar.snv,">",tmp.somatic1)
system(cmd)

cmd = paste("cat",tmp.somatic1,">>",clinvar.somatic)
system(cmd)

cmdrm = paste("rm",tmp.somatic1)
system(cmdrm)

message("# 6) Get dbSNP with and without clinvar SNPs")
nosomatic = gsub(out2,pattern = "\\.vcf",replacement = "_clinvarsomatic_excluded.vcf")
cmd = paste("bedtools intersect -a",out2,"-b",clinvar.somatic,"-wa -v -header >",nosomatic)
system(cmd)

message("# 7) Split dbSNPs by chromosome")
vcf2chrom = "Rscript /home/nicola.casiraghi/Scripts/dbSNP_adjust/VCFbyChromosome.R"
bychr = file.path(outfolder,"ByChromosome")
dir.create(bychr)
cmd = paste(vcf2chrom,nosomatic,bychr)
system(cmd)
}

proc.time()-timestart


