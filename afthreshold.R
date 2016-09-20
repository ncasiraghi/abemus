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

vafcov_file = file.path(outdir,"BaseErrorModel","afgtz.tsv")
afz_file = file.path(outdir,"BaseErrorModel","afz.RData")

#  check if data tables exists
if(file.exists(file.path(outdir,"BaseErrorModel","afgtz.tsv"))){
  cat(paste("[",Sys.time(),"]\tlooking for data table with AFs > 0 and coverages:",vafcov_file,"[ ok ]"),"\n")  
} else {
  cat(paste("[",Sys.time(),"]\tlooking for data table with AFs > 0 and coverages:",vafcov_file,"[ NOT found ]"),"\n") 
  quit()
}

if(file.exists(file.path(outdir,"BaseErrorModel","afz.tsv"))){
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

minaf_cov = th_results_bin[which(th_results_bin$specificity==spec),]

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
  
# Estimate AF threshold for not available coverage bins.
cat(paste("\n[",Sys.time(),"]\tEstimation of AF thresholds for high coverages\n"))  

# find out most represented bin of coverage
cat(paste("[",Sys.time(),"]\tSelecting most represented coverage bin\n"))  
datacount_bin_complete <- datacount_bin+afz
datacount_bin_complete <- datacount_bin_complete[grep(pattern = "Inf",names(datacount_bin_complete),invert = T,value = T)]
maxbin = names(which.max(datacount_bin_complete))
maxbin_count_gtz = as.numeric(datacount_bin[maxbin])
maxbin_count_etz = as.numeric(afz[maxbin])
maxbin_count_total = maxbin_count_gtz+maxbin_count_etz
cat(paste("[",Sys.time(),"]\tMost represented coverage bin:\t",maxbin),"\n")  
cat(paste("[",Sys.time(),"]\tData points in coverage bin",maxbin,":\t",maxbin_count_total),"\n")
idx = which(a==maxbin)
b = vafcov[idx,1]

# subsampling coverage bin and compute quantiles
cat(paste("[",Sys.time(),"]\tSubsampling max coverage bin and compute quantiles\n"))  

sampling = data.frame(bin = maxbin,
                      tot.pos = maxbin_count_total,
                      etz.pos = maxbin_count_etz,
                      etz.frc = maxbin_count_etz/maxbin_count_total,
                      gtz.pos = maxbin_count_gtz,
                      gtz.frc = maxbin_count_gtz/maxbin_count_total,
                      sampid = 1,
                      stringsAsFactors = F)

tot.pos = maxbin_count_total
threhsold_sampling = data.frame(spec=probs,row.names = probs,stringsAsFactors = F)
threhsold_sampling = cbind(threhsold_sampling,th.original=quantile.zaf(x = b,probs = probs,nz = maxbin_count_etz))

thsamp = list(threhsold_sampling)
times = 1000 
keepsubf = 1/2

sampid = 2
tot.pos = as.integer(tot.pos*keepsubf)
while(tot.pos>10){
  cat(paste("[",Sys.time(),"]\tSubsampling",tot.pos,"positions from total max bin\n"))
  npos_etz_tosample = as.integer(tot.pos*sampling$etz.frc[1])
  npos_gtz_tosample = as.integer(tot.pos*sampling$gtz.frc[1])
  this_sampling = data.frame(bin = maxbin,
                             tot.pos = tot.pos,
                             etz.pos = npos_etz_tosample,
                             etz.frc = npos_etz_tosample/tot.pos,
                             gtz.pos = npos_gtz_tosample,
                             gtz.frc = npos_gtz_tosample/tot.pos,
                             sampid = sampid,
                             stringsAsFactors = F)
  sampling = rbind(sampling,this_sampling)
  this_data = data.frame(spec=probs,row.names = probs,stringsAsFactors = F)
  cat(paste("[",Sys.time(),"]\tComputing",times,"replicas sampling","\n"))  
  for(n in 1:times){
    # pick up AFs indicated in sampling from AF tables in correct coverage bin
    bsamp = sample(x = b,size = npos_gtz_tosample,replace = F)
    # compute AF thresholds
    replica = quantile.zaf(x = bsamp,probs = probs,nz = npos_etz_tosample)
    this_data[,paste0("replica.",n)] = replica
  }
  this_data$mean = apply(this_data[,-1],MARGIN = 1,FUN = mean,na.rm=T)
  this_data$median = apply(this_data[,-1],MARGIN = 1,FUN = median,na.rm=T)
  this_data$std = apply(this_data[,-1],MARGIN = 1,FUN = sd,na.rm=T)
  this_data$coeffvar = (this_data$std)/(this_data$mean)
  thsamp[[sampid]] <- this_data 
  sampid = sampid+1
  tot.pos = as.integer(tot.pos*keepsubf)
}

tabsummary.mean = data.frame(spec=probs,row.names = probs,stringsAsFactors = F)
tabsummary.std = data.frame(spec=probs,row.names = probs,stringsAsFactors = F)
for(sampid in sampling$sampid){
  tt = thsamp[[sampid]]
  tabsummary.mean[,as.character(sampid)]<-tt$mean
  tabsummary.std[,as.character(sampid)]<-tt$std
}
tabsummary.coeffvar = tabsummary.std/tabsummary.mean
tabsummary.coeffvar$spec <- rownames(tabsummary.coeffvar)
  
cat(paste("[",Sys.time(),"]\tSaving threshold data from random sampling\n"))
save(sampling,thsamp,tabsummary.mean,tabsummary.std,tabsummary.coeffvar,datacount_bin_complete,file = file.path(outdir,"Controls",paste0("datathreshold_random_sampling.RData")),compress = T)

pdf(file = file.path(outdir,"Controls","datathreshold_random_sampling.pdf"))
par(oma=c(5,5,2,5))
for(currentspec in 1:length(probs)){
  layout(mat = matrix(c(1,1,1,2,2,2),nrow = 2,byrow = T))
  ylim = range(as.matrix(tabsummary.mean[,-1]))
  original = thsamp[[1]]
  original = original[currentspec,2]
  dataspec = as.numeric(tabsummary.mean[currentspec,-1])
  dataspec.std = as.numeric(tabsummary.std[currentspec,-1])
  dataspec.cfv = as.numeric(tabsummary.coeffvar[currentspec,-1])
  xxx = barplot(height = c(0,dataspec.cfv),plot = F,space=0)
  par(mar=c(1,5,3,3))
  plot(x = xxx,c(original,dataspec),las=2,xaxt="n",ann = F,ylim = ylim)
  errbar(x = xxx,y = c(original,dataspec),yplus = c(original,dataspec)+c(0,dataspec.std),yminus = c(original,dataspec)-c(0,dataspec.std),add = T,col = "grey40",type = 'n')
  abline(v = xxx,lwd=0.5,lty=3,col="grey30")
  abline(h = original,lwd=0.7,col="grey40")
  title(ylab="Allelic Fraction threshold",line=4,cex.lab = 1.2)
  title(main = paste(maxbin,probs[currentspec],original,sep="   "),line = 1)
  par(mar=c(8,5,1,3))
  plot(x = xxx,c(0,dataspec.cfv),xaxt="n",ann=F,las=2,pch=16,col="grey40",ylim=c(0,1))
  axis(side = 1,at = xxx,labels = sampling$tot.pos,las=2,cex.axis = 1.2)
  abline(v = xxx,lwd=0.5,lty=3,col="grey30")
  title(ylab="Coefficient of variation",line=4,cex.lab = 1.2)
  title(xlab="Number of data points sampled",line=7,cex.lab = 1.2)
  text(x = xxx,labels = c(0,round(dataspec.cfv,5)),y = 0.5,srt = 90,adj = 0)
}

dev.off()

#  Select the min number of required positions to get a stable AF threshold
cat(paste("[",Sys.time(),"]\tSelection of minimum number of positions to have a stable AF threshold\n"))

mincoeffvar = 0.01

npos = as.numeric(tabsummary.coeffvar[as.character(spec),-1])
names(npos) = sampling$tot.pos[2:nrow(sampling)]
min_npos = min(as.numeric(names(npos[which(npos<=mincoeffvar)])))
cat(paste("[",Sys.time(),"]\tMin number of positions:\t",min_npos,"\n"))

th_results_bin_minpos = th_results_bin[,names(which(datacount_bin_complete>=min_npos)),drop=F]

lastbin = names(th_results_bin_minpos)[ncol(th_results_bin_minpos)] 
lastbin = paste0(unlist(strsplit(lastbin,split = ","))[1],",Inf]")
colnames(th_results_bin_minpos)[ncol(th_results_bin_minpos)] <- lastbin

cat(paste("[",Sys.time(),"]\tSaving threshold data\n"))
th_results_bin_minpos = cbind(th_results_bin$specificity,th_results_bin_minpos)
colnames(th_results_bin_minpos)[1] = "specificity"
minaf_cov_minpos = th_results_bin_minpos[which(th_results_bin_minpos$specificity == spec),]

covbin_minpos = covbin[1:length(minaf_cov_minpos[,-1])]
covbin_minpos = c(covbin_minpos,Inf)
save(min_npos,covbin_minpos,minaf_cov_minpos,th_results_bin_minpos,file = file.path(outdir,"Controls",paste0("datathreshold_high_coverages.RData")),compress = T)

proc.time()-timestart

