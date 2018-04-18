options(scipen=99999,digits = 10)
library(parallel)
library(data.table)
library(bigmemory)
library(Hmisc)

quantile.zaf <- function (x, probs = seq(0, 1, 0.25),nz){
  N <- length(x) + round(nz)
  x <- sort(x)
  qnt <- 1 + round((probs * (N - 1))) # !
  ql <- c()
  for(j in qnt){
    if(j-nz <=0){
      ql = c(ql,0)
    } else {
      ql = c(ql,x[j-nz])
    }
  }
  names(ql) <- probs
  return(ql)
}

strandbias = function(fwd.ref,fwd.alt,rev.ref,rev.alt){
  if(fwd.alt+rev.alt==0 | fwd.ref+rev.ref==0){
    return(NA)
  } else {
    #sb = abs((fwd.alt/(fwd.ref+fwd.alt))-(rev.alt/(rev.ref+rev.alt)))/((fwd.alt+rev.alt)/sum(fwd.ref,fwd.alt,rev.ref,rev.alt))
    sb = abs((fwd.alt/(fwd.ref+fwd.alt))-(rev.alt/(rev.ref+rev.alt)))
    sb = round(sb,5)
    return(sb)
  }
}

strandbias_fisher = function(fwd.ref,fwd.alt,rev.ref,rev.alt){
  fisher = fisher.test(cbind(c(fwd.ref,fwd.alt),c(rev.ref,rev.alt)))$p.value
  return(fisher)
}

strandbias_OR = function(fwd.ref,fwd.alt,rev.ref,rev.alt){
  R = (fwd.ref*rev.alt)/(fwd.alt*rev.ref)
  or = sum(R,1/R)  
  return(or)
}

strandbias_perc = function(fwd.ref,fwd.alt,rev.ref,rev.alt){
  dperc = abs(fwd.ref/rev.ref-fwd.alt/rev.alt)
  return(dperc)
}

pos2bperr = function(id,targets,germlineset,step,chrom,lev,covbin,af_max_to_compute_thresholds,coverage_min_to_compute_thresholds,af_max_to_compute_pbem,coverage_min_to_compute_pbem,n_pos_af_th){
  upto = id+step-1
  if(upto>nrow(targets)){upto <- nrow(targets)}
  this = targets[id:upto,,drop=F]
  outfile = paste(this$chr[1],this$pos[1],this$pos[nrow(this)],"postogrep.txt",sep="_")
  filtpileup = paste(this$chr[1],this$pos[1],this$pos[nrow(this)],"filtered.pileup.txt",sep="_")
  taboutchrom = paste(this$chr[1],this$pos[1],this$pos[nrow(this)],"pbem.table.txt",sep="_")
  mytabafchrom = paste(this$chr[1],this$pos[1],this$pos[nrow(this)],"afgtz.table.txt",sep="_")
  afzchrom = paste(this$chr[1],this$pos[1],this$pos[nrow(this)],"afz.table.txt",sep="_")
  afz = array(data = 0,dim = length(lev),dimnames = list(lev))
  cat(unique(this$pos),sep = "\n",file = outfile,append = F)
  cmd = paste0("awk -F'\t' '{if (FILENAME == \"",outfile,"\") { t[$1] = 1; } else { if (t[$2]) { print }}}' ",outfile," ",germlineset," > ",filtpileup)
  system(cmd)
  if(file.info(filtpileup)$size == 0){
    cat(paste("[",Sys.time(),"]\tpositions in ",outfile,"not found in any pileups."),"\n")
  } else {
    completetab_all = fread(filtpileup,stringsAsFactors = F,showProgress = T,header = F,na.strings = "",colClasses=list(character=10))
    names(completetab_all)=c("chr","pos","ref","A","C","G","T","af","RD","dbsnp")
    # exclude annotated and private SNPs [ to compute AF threshold]
    completetab = completetab_all[which(is.na(completetab_all$dbsnp)),,drop=F]
    completetab = completetab[which(completetab$af <= af_max_to_compute_thresholds & completetab$RD >= coverage_min_to_compute_thresholds),,drop=F]
    # Compute pbem also in positions that are SNPs
    completetab_dbsnp = completetab_all[which(completetab_all$af <= af_max_to_compute_pbem & completetab_all$RD >= coverage_min_to_compute_pbem),,drop=F]
    if(nrow(completetab)>0){
      # save allelic fractions by bins of coverage
      mytabaf = completetab[which(completetab$af > 0),,drop=F]
      mytabz =  completetab[which(completetab$af == 0),,drop=F]
      afz = afz + as.array(table(cut(mytabz$RD,breaks = covbin,include.lowest = T)))
      write.table(t(afz),file = afzchrom,sep="\t",col.names = F,row.names = F,quote = F,append = F)
      write.table(mytabaf[,8:9,with=F],file = mytabafchrom,append = F,sep = "\t",quote = F,row.names = F,col.names = F)
      # compute pbem
      completetab_dbsnp$group = paste(completetab_dbsnp$chr,completetab_dbsnp$pos,completetab_dbsnp$ref,sep=":")
      dat <- data.table(completetab_dbsnp, key="group")
      ans <- dat[,{list(tot_coverage=sum(RD,na.rm = T),
                        total.A=sum(A,na.rm = T),
                        total.C=sum(C,na.rm = T),
                        total.G=sum(G,na.rm = T),
                        total.T=sum(T,na.rm = T),
                        n_pos_available = length(which(RD > 0)),
                        n_pos_af_lth=length(which(af < n_pos_af_th)),
                        n_pos_af_gth=length(which(af >= n_pos_af_th)),
                        count.A_af_gth=length(A[which(af >= n_pos_af_th & A>0)]),
                        count.C_af_gth=length(C[which(af >= n_pos_af_th & C>0)]),
                        count.G_af_gth=length(G[which(af >= n_pos_af_th & G>0)]),
                        count.T_af_gth=length(T[which(af >= n_pos_af_th & T>0)]))},by="group"]
      this = as.data.table(this)
      this$group = paste(this$chr,this$pos,this$ref,sep=":")
      this_ans = merge(x = this,y = ans,by = "group",all.y = T)
      this_ans = as.data.frame(this_ans)
      tabstats = fromListToDF(mclapply(seq(1,nrow(this_ans),1),perbaserror_stats,tgf=this_ans,mc.cores = 2 ))
      tabstats = tabstats[with(tabstats, order(pos)), ]
      cat(paste("[",Sys.time(),"]\tWriting output for positions in: ",filtpileup),"\n")
      write.table(tabstats,file = taboutchrom,append = F,quote = F,row.names = F,col.names = F,sep="\t")
    }
  }
  system(paste("rm",filtpileup,outfile))
}  

get_germlineset <- function(TableSif){
  germlineset = c()
  for(id in 1:nrow(TableSif)){
    thisSample=TableSif[id,]
    name = gsub(basename(thisSample$germline.bam),pattern = ".bam",replacement = "")
    thisPileup.file = list.files(file.path(pacbamfolder,name,"pileup"),full.names = T,pattern = paste0(chrom,"\\.pileup$"))
    germlineset = c(germlineset,thisPileup.file)
  }  
  germlineset = paste(unique(germlineset),collapse = " ")
  return(germlineset)
}

perbaserror_stats = function(id,tgf){
  this = tgf[id,,drop=F]
  alts = setdiff(c("total.A","total.C","total.G","total.T"),paste0("total.",this$ref))
  rd.alts = as.numeric(sum(this[alts],na.rm = T))
  this$bperr = sum(rd.alts)/this$tot_coverage
  this$tot_reads_supporting_alt = rd.alts
  return(this)
}

CheckAltReads <- function(i,snvs){
  this = snvs[i,,drop=F]
  altbase = this$alt
  refbase = this$ref
  cov.ref = this[,refbase]
  rev.ref = this[,paste0(refbase,"rs")]
  fwd.ref = cov.ref - rev.ref
  this$rev.ref = rev.ref
  this$fwd.ref = fwd.ref
  if(altbase=="N"){
    altbases = setdiff(c("A","C","G","T"),refbase) 
    cov.alt = as.numeric(max(this[,altbases]))
    this$cov.alt <- cov.alt
    this$af <- round(cov.alt/sum(cov.alt,this[,refbase]),4)
    this$rev.alt = NA
    this$fwd.alt = NA
    this$strandbias = NA
  } else {
    cov.alt = this[,altbase]
    this$cov.alt <- cov.alt
    rev.alt = this[,paste0(altbase,"rs")]
    fwd.alt = cov.alt - rev.alt
    this$rev.alt = rev.alt
    this$fwd.alt = fwd.alt
    this$strandbias = strandbias(fwd.ref=fwd.ref,fwd.alt=fwd.alt,rev.ref=rev.ref,rev.alt=rev.alt)
  }
  return(this)
}

fromListToDF <- function(inputList){
  if (is.null(inputList)){return(NULL)}
  #check if some is null and remove
  nullPositions <- which(sapply(inputList,is.null))
  if (length(nullPositions) > 0){
    inputList <- inputList[-nullPositions]    
  }
  firstEl <- inputList[[1]]  
  inputList <- lapply(inputList, function(x){ matrix(unlist(x), ncol=ncol(x))} )
  outDF <- as.data.frame(do.call(rbind, inputList),stringsAsFactors=F)
  colnames(outDF) <-names(firstEl)
  
  for(idx in c(1:ncol(outDF))){
    if (class(firstEl[[idx]]) == "logical"){       
      if (is.na(firstEl[[idx]])){
        class(outDF[[idx]]) <- "numeric"
      }else if (outDF[1,idx] == "TRUE" ||outDF[1,idx] == "FALSE"  ){
        outDF[,idx] <- as.logical(outDF[,idx]) * 1      
      }
      class(outDF[[idx]]) <- "numeric"
    }else{
      class(outDF[[idx]]) <- class(firstEl[[idx]])
    }
  }
  return(outDF)
}

compute_proptest <- function(id,dataf){
  locus=dataf[id,,drop=F]
  if(locus$bperr>0){
    locus$proptest = prop.test(x = c(locus$cov.alt,locus$tot_reads_supporting_alt),n = c(locus$cov_case,locus$tot_coverage))$p.value
  }
  return(locus)
}

compute_pbem_allele <- function(id,abemus){
  locus = abemus[id,,drop=F]
  locus$pbem_allele = locus[,paste0('total.',locus$alt)]/locus$tot_coverage
  # check if same alt allele is observed in case and germline (when af germline > 0)
  alt_case = as.numeric(locus[,paste0(locus$alt,'_case')])
  alt_ctrl = as.numeric(locus[,paste0(locus$alt,'_control')])
  if(alt_ctrl>0){
    locus$same_allele = 1
  } else {
    locus$same_allele = 0
  }
  return(locus)
}

getAltCountCovBased = function(id,target){
  locus = target[id,]
  out = locus[,c('total.A','total.C','total.G','total.T')]
  out[,paste0('total.',locus$ref)] = 0
  names(out) = c('cov.total.A','cov.total.C','cov.total.G','cov.total.T')
  out.count <- out
  out.count[which(out.count>0)] = 1
  names(out.count) = c('n.total.A','n.total.C','n.total.G','n.total.T')
  locusout = cbind(out,out.count)
  return(locusout)
}

getTransTrasvCount = function(id,target){
  locus = target[id,]
  ref = locus$ref
  out = locus[,c('total.A','total.C','total.G','total.T')]
  out[,paste0('total.',locus$ref)] = 0
  alts=gsub(names(out[,which(out>0)]),pattern = 'total.',replacement = '')
  if(length(alts)>0){
    out = paste(ref,alts,sep = '>')
    return(out)
  } else {
    return(NA)
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

add_class_xbg = function(pmtab,xbg){
  if(nrow(pmtab)>0){
    pmtab$CLASS.xbg = NA
    pmtab$CLASS.xbg[which(pmtab$af_control<=xbg & pmtab$bperr<=xbg & pmtab$pbem_allele<=xbg)] = 1
    pmtab$CLASS.xbg[which(pmtab$af_control<=xbg & pmtab$bperr>xbg & pmtab$pbem_allele<=xbg)] = 2
    pmtab$CLASS.xbg[which(pmtab$af_control<=xbg & pmtab$bperr>xbg & pmtab$pbem_allele>xbg)] = 3
    pmtab$CLASS.xbg[which(pmtab$af_control>xbg & pmtab$bperr>xbg & pmtab$pbem_allele>=xbg & pmtab$same_allele == 0)] = 4
    pmtab$CLASS.xbg[which(pmtab$af_control>xbg & pmtab$bperr>xbg & pmtab$pbem_allele>xbg & pmtab$same_allele == 1)] = 5
    return(pmtab)
  } else {
    return(pmtab)
  }
}
