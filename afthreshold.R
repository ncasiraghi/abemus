#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=1){
  message("\nError!\nusage: Rscript afthreshold.R abemus_configure.R\n")
  quit()
}
source(args[1])
abemus_functions = file.path(abemusdir,"functions.R")
source(abemus_functions)

cat(paste("[",Sys.time(),"]\tEstimation of allelic fraction thresold exploiting germline samples only","\n"))

if(!file.exists(file.path(outdir, "Controls"))){
  dir.create(file.path(outdir, "Controls"), showWarnings = TRUE)
}
setwd(file.path(outdir, "Controls"))
timestart <- proc.time()

# get chromosomes
chromosomes = read.delim(file = file.path(targetbp,"bpcovered.tsv"),as.is=T)
chromosomes = sort(paste0("chr",chromosomes[-nrow(chromosomes),1]))

covbin = seq(0,5000,by = coverage_binning)
covbin[length(covbin)] <- Inf
lev = levels(cut(1,breaks=covbin,include.lowest=TRUE))
probs = seq(0.99,1,0.0001)

vafcov_file = file.path(pbem_dir,"afgtz.tsv")
afz_file = file.path(pbem_dir,"afz.RData")

#  check if data tables exists
if(file.exists(file.path(pbem_dir,"afgtz.tsv"))){
  cat(paste("[",Sys.time(),"]\tlooking for data table with AFs > 0 and coverages:",vafcov_file,"[ ok ]"),"\n")  
} else {
  cat(paste("[",Sys.time(),"]\tlooking for data table with AFs > 0 and coverages:",vafcov_file,"[ NOT found ]"),"\n") 
  quit()
}

if(file.exists(file.path(pbem_dir,"afz.RData"))){
  cat(paste("[",Sys.time(),"]\tlooking for data table with AFs = 0 and coverages:",afz_file,"[ ok ]"),"\n")  
} else {
  cat(paste("[",Sys.time(),"]\tlooking for data table with AFs = 0 and coverages:",afz_file,"[ NOT found ]"),"\n") 
  quit()
}

# Load AFs tables and compute thresholds based on distribution quantiles.

cat(paste("[",Sys.time(),"]\tLoading data table with AFs > 0 and coverages:",vafcov_file),"\n")  
vafcov = read.big.matrix(filename = vafcov_file,header = F,sep = "\t",type = "double")
cat(paste("[",Sys.time(),"]\tLoading data table with AFs = 0 and coverages:",afz_file),"\n")
load(afz_file)
cat(paste("[",Sys.time() ,"]\taf.all",dim(vafcov)[1] + sum(afz),"af.gtz",dim(vafcov)[1],"af.etz",sum(afz)),"\n")

cat(paste("[",Sys.time(),"]\tcompute AF quantiles ( not stratified by coverage )"),"\n")
nz = sum(afz)
th_results = quantile.zaf(x = vafcov[,1],probs = probs,nz = nz)
minaf = as.numeric(th_results[as.character(spec)])

cat(paste("[",Sys.time(),"]\tcompute AF quantiles ( stratified by coverage )"),"\n")
a <- cut(vafcov[,2],breaks = covbin,include.lowest = T)
th_results_bin = data.frame(specificity = probs)
datacount_bin = array(data = 0,dim = length(lev),dimnames = list(lev))
for(l in lev){
  cat(paste("[",Sys.time(),"]\tprocessing coverage bin\t",l),"\n")
  idx = which(a==l)
  datacount_bin[l] <- length(idx)
  if(length(idx)>0){
    b = vafcov[idx,1]
    nz = as.integer(afz[l])
    kk = quantile.zaf(x = b,probs = probs,nz = nz)
    th_results_bin[,l] <- kk
  } else {
    th_results_bin[,l] <- NA
  }
}

minaf_cov = th_results_bin[which(round(th_results_bin$specificity,4)==spec),]

cat(paste("[",Sys.time(),"]\tSaving threshold data\n"))
save(minaf,minaf_cov,th_results,th_results_bin,covbin,lev,datacount_bin,file = file.path(outdir,"Controls",paste0("datathreshold.RData")),compress = T)

# Plot HeatMap data

cat(paste("[",Sys.time(),"]\tPlotting results"),"\n")
load(file = file.path(outdir,"Controls",paste0("datathreshold.RData")))
pdf(file = file.path(outdir,"Controls",paste0("datathreshold_heatmap.pdf")), h = 219, w = 297, paper='A4r')
layout(mat = matrix(c(1,1,1,
                      2,2,2,
                      2,2,2,
                      2,2,2,
                      2,2,2,
                      3,3,3,
                      3,3,3),nrow = 7,byrow = T))
par(xaxs = "i") 
th_results_bin_matrix = as.matrix(th_results_bin[,-1])
keepcols = colSums(is.na(th_results_bin_matrix))<nrow(th_results_bin_matrix)
th_results_bin_matrix <- th_results_bin_matrix[,keepcols]
datacount_bin_filt <- datacount_bin[names(which(keepcols))]
datacount_bin_filt <- datacount_bin_filt+afz[names(which(keepcols))]
maxaf = max(th_results_bin_matrix,na.rm = T)
# plot barplot 
par(mar=c(0,6,0,0))
x_coordinates <- barplot(datacount_bin_filt,space=0,col = "grey70",border = "white",xaxt="n",las=2,log="y")
text(x = x_coordinates,y = 1e+04,labels = datacount_bin_filt,srt=90,adj = 0)
# plot heat map
N = 1000
my_palette <- colorRampPalette(rev(c("#d7191c","#fdae61","#ffffbf","#abdda4","#2b83ba")))(n = N)
par(mar=c(0,6,1,0))
image(x = x_coordinates,y=probs,z = t(th_results_bin_matrix),breaks = seq(0.0,maxaf, length.out = N+1),col=my_palette,axes=F,ann=F)
axis(side = 2,at=probs[seq(1,length(probs),5)],labels = probs[seq(1,length(probs),5)],las=2)
# boxplot
par(mar=c(6,6,0,0))
boxplot(at = x_coordinates,th_results_bin_matrix,pch=20,outcol="#636363",las=2,ylim=c(0,maxaf+sd(th_results_bin_matrix,na.rm = T)),frame.plot=FALSE)
dev.off()
  
proc.time()-timestart

