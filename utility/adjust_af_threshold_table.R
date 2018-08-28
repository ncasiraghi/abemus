#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=4){
  message("\nERROR:\t2 argument required")
  message("\nUSAGE:\tRscript adjust_af_threshold_table.R <abemus_configure.R> <replicas> <n.cores> <detection.spec>\n")
  cat(paste("\t[1] abemus_configure.R","\n"))
  cat(paste("\t[2] replicas","\n"))
  cat(paste("\t[3] replicas in parallel","\n"))
  cat(paste("\t[4] detection spec","\n\n"))
  quit()
}

# input data
config = args[1]
#config="/scratch/sharedCO/Casiraghi/Abemus_data/InSilicoData/HALO_t1/abemus_configure.R"
replicas = as.numeric(args[2])
#replicas=10
replicas.in.parallel = as.numeric(args[3])
#replicas.in.parallel = 2
this.detection.specificity = as.numeric(args[4]) 
#this.detection.specificity = 0.995

source(file = config)
setwd(controls_dir)

abemus_functions = file.path(abemusdir,"functions.R")
source(abemus_functions)

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

# Load data from controls folder
load(file = file.path(controls_dir,"datathreshold.RData"))

a <- cut(vafcov[,2],breaks = covbin,include.lowest = T)

# Find out most represented bin of coverage

datacount_bin_complete = datacount_bin+afz
datacount_bin_complete = datacount_bin_complete[grep(pattern = "Inf",names(datacount_bin_complete),invert = T,value = T)]
datacount_bin_complete = sort(datacount_bin_complete,decreasing = T)
datacount_bin_complete = datacount_bin_complete[which(datacount_bin_complete>0)]

stop = length(datacount_bin_complete)-1
last.stable.card = 0
tab = c()

get_afth <- function(num,idx,next.bin.AFgtz,vafcov,this.detection.specificity){
  if(length(idx) < next.bin.AFgtz){
    b = vafcov[sample(x = idx,size = next.bin.AFgtz,replace = T),1]
  } else {
    b = vafcov[sample(x = idx,size = next.bin.AFgtz,replace = F),1]
  }
  afth = quantile.zaf(x = b,probs = this.detection.specificity,nz = next.bin.AFetz)
  return(afth)
}

for(i in 1:stop){
  current.bin = names(datacount_bin_complete)[i]
  current.bin.card = as.numeric(datacount_bin_complete[i])
  message("current bin: ",current.bin)
  if(current.bin.card >= last.stable.card){
    for(j in (i+1):(stop+1)){
      next.bin = names(datacount_bin_complete)[j]
      next.bin.card = as.numeric(datacount_bin_complete[j])
      message("sampling according to: ",next.bin)
      next.bin.AFetz = as.numeric(afz[next.bin])
      next.bin.AFgtz = as.numeric(datacount_bin[next.bin])
      ReplicasTable = c()
      idx = which(a==current.bin)
      out = mclapply(1:replicas,get_afth,idx=idx,next.bin.AFgtz=next.bin.AFgtz,vafcov=vafcov,this.detection.specificity=this.detection.specificity,mc.cores = replicas.in.parallel)
      ReplicasTable = matrix(unlist(out),ncol = replicas)
      coeffvar = apply(ReplicasTable,MARGIN = 1,FUN = sd,na.rm=T)/apply(ReplicasTable,MARGIN = 1,FUN = mean,na.rm=T)
      if(is.na(coeffvar)){
        coeffvar <- 0
      }
      x = data.frame(current.bin=current.bin,next.bin=next.bin,coeffvar=as.numeric(coeffvar),median.afth=median(ReplicasTable,na.rm = T),last.stable.card=last.stable.card,stringsAsFactors = F)
      tab=rbind(tab,x)
      if(coeffvar > 0.01){
        last.stable.card = max(last.stable.card,next.bin.card)
        break
      }
    }
  }
}

message("ROSPO")

name.out = paste0("tab",replicas,"r_",this.detection.specificity,".RData")
save(tab,file = file.path(controls_dir,name.out),compress = T)

#load(file = replica.rdata)
message("GRIGIO")

# correct the threhsold table

minaf_cov_corrected <- th_results_bin[which(th_results_bin$specificity == this.detection.specificity),]
C = tab$last.stable.card[nrow(tab)]

last.afth.used = 1
start = 2
end = length(minaf_cov_corrected)-1

for(i in start:end){
  #message(i)
  bin.name = names(minaf_cov_corrected)[i]
  bin.card = as.numeric(datacount_bin_complete[which(names(datacount_bin_complete)==bin.name)])
  if(!identical(bin.card,numeric(0))){
    if(bin.card < C & last.afth.used == 1){
      minaf_cov_corrected[i] <- last.afth.used
    }
    if(bin.card >= C){
      last.afth.used <- minaf_cov_corrected[i]
    }
    if(bin.card < C & last.afth.used != 1){
      minaf_cov_corrected[i] <- last.afth.used
    }
  }
}

minaf_cov_corrected[which(is.na(minaf_cov_corrected))] <- last.afth.used
minaf_cov_corrected[which(minaf_cov_corrected==1)] <- NA

# correct the AF threhsold in bin where there is Inf as limit
N = length(minaf_cov_corrected)
minaf_cov_corrected[N] <- minaf_cov_corrected[N-1]

name.out = paste0("minaf_cov_corrected_",replicas,"r_",this.detection.specificity,".RData")
save(minaf_cov_corrected,file = file.path(controls_dir,name.out))

message("IN CATTIVITA")

pdf.out = paste0("minaf_cov_corrected_",replicas,"r_",this.detection.specificity,".pdf")
pdf(pdf.out, h = 219, w = 297, paper='A4r')

par(mfrow=c(3,1),oma=c(2,2,0,0))
a = as.numeric(th_results_bin[which(th_results_bin$specificity==this.detection.specificity),2:length(th_results_bin)])
b = as.numeric(minaf_cov_corrected[2:length(minaf_cov_corrected)])

count = rbind(as.numeric(datacount_bin),as.numeric(afz))
total.pos = datacount_bin+afz

col = rep("grey",ncol(count))
col[which(total.pos>=C)] <- "forestgreen"

xcord=barplot(count,las=2,main=this.detection.specificity,col=c("grey60","grey90"),border = c("grey60","grey90"),names.arg = names(total.pos),ylab = "n. of positions")
mtext(controls_dir,cex = 0.6)
text(x = xcord[which(total.pos>=C)],y = total.pos[which(total.pos>=C)],pos = 1,labels = "*",cex = 1.5,col = "forestgreen")
abline(h = C,lwd=0.5,lty=2)

barplot(a,las=2,ylab="Automatic AF threshold",col=col,border = col,names.arg = round(a,4))
barplot(b,las=2,ylab="Adjusted AF threshold",ylim=par("yaxp")[1:2],col=col,border=col,names.arg = round(b,4))

dev.off()
