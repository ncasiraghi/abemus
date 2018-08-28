library(data.table)
library(parallel)

mc.cores = 30

args <- commandArgs(trailingOnly = TRUE)

# pileups to scale
#folder.pacbam = "/scratch/sharedCO/Casiraghi/InSilicoData_Tumors/HALO_byChrom"
#list.of.samples = list.files(folder.pacbam,full.names = T,recursive = F)
#list.of.samples = readLines("/scratch/sharedCO/Casiraghi/InSilicoData_Tumors/HALO_byChrom_c10_model/list_tumors.txt")
list.of.samples = readLines(args[1])
list.of.pileups = list.files(file.path(list.of.samples,"pileup"),full.names = T,recursive = F)

# Base Error model on real data
pbem.model.bychrom = list.files("/scratch/sharedCO/Casiraghi/Abemus_data/IPM_tissues_Haloplex/BaseErrorModel",pattern = "bperr_chr",full.names = T)

# Out folder where put scaled pileups
out.folder = "/scratch/sharedCO/Casiraghi/InSilicoData_Tumors/HALO_byChrom_c10_model/"

# Functions
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

scaling <- function(i,p,cov.scaled,pbem.chr){
  pos = p[i,,drop=F]
  if(!is.na(pos$bperr)){
    # apply scaling factor on coverage
    pos$cov = round(pos$cov*cov.scaled) 
    # how many reads with error
    n.of.success = rbinom(size = pos$cov,n = 1,prob = pos$bperr)
    # adjust coverage of the reference
    pos[,pos$ref] = pos$cov - n.of.success
    alts = setdiff(c("A","C","G","T"),pos$ref)
    pos[,alts] = 0
    # how many errors for each not-reference allele
    if(n.of.success > 0){
      outs = sample(x = alts,size = n.of.success,replace = T,prob = pos[,paste0("total.",alts)]/pos$tot_coverage)
      pos[,alts[1]] = length(which(outs==alts[1]))
      pos[,alts[2]] = length(which(outs==alts[2]))
      pos[,alts[3]] = length(which(outs==alts[3]))
    }
    # compute af
    pos$af = sum(pos[,alts])/pos$cov
  } else {
    # apply scaling factor on coverage
    pos$cov = round(pos$cov*cov.scaled) 
    # how many reads with error
    n.of.success = rbinom(size = pos$cov,n = 1,prob = 0.000181845800784122)
    # adjust coverage of the reference
    pos[,pos$ref] = pos$cov - n.of.success
    alts = setdiff(c("A","C","G","T"),pos$ref)
    pos[,alts] = 0
    # how many errors for each not-reference allele
    if(n.of.success > 0){
      outs = sample(x = alts,size = n.of.success,replace = T)
      pos[,alts[1]] = length(which(outs==alts[1]))
      pos[,alts[2]] = length(which(outs==alts[2]))
      pos[,alts[3]] = length(which(outs==alts[3]))
    }
    # compute af
    pos$af = sum(pos[,alts])/pos$cov
  }
  return(pos[,c("chr","pos","ref","A","C","G","T","af","cov","dbsnp")])
} 

# Compute scaling
for(id in list.of.samples){
  cat(paste("[",Sys.time(),"]\t >>",basename(id),"\n"))
  # create pileup and snvs folders
  dir.create(path = file.path(out.folder,basename(id)))
  dir.create(path = file.path(out.folder,basename(id),"pileup"))
  dir.create(path = file.path(out.folder,basename(id),"snvs"))
  # scale pileup
  ps = grep(x = list.of.pileups,pattern = basename(id),value = T)
  for(x in ps){
    cat(paste("[",Sys.time(),"]\t",basename(x),"\n"))
    # split pileup file
    split.folder = file.path(out.folder,basename(id),"pileup","split.folder")
    dir.create(path = split.folder)
    setwd(split.folder)
    cmd = paste("tail -n +2",x,"| split - -l 100000")
    system(cmd)
    # load chr pbem 
    chr=gsub(unlist(strsplit(basename(x),split = "_"))[3],pattern = "\\.pileup",replacement = "")
    input = grep(pbem.model.bychrom,pattern = paste0(chr,".tsv$"),value = T)
    pbem.chr = fread(input = input,stringsAsFactors = F,showProgress = F,header = F,data.table = T,select = c(1,10:14,22))
    colnames(pbem.chr) <- c("group","tot_coverage","total.A","total.C","total.G","total.T","bperr")
    splt = list.files(split.folder,full.names = T)
    for(y in splt){
      # load pileup 
      p = fread(input = y,stringsAsFactors = F,showProgress = F,header = F,na.strings = "",colClasses = list(character=10),data.table = F)
      colnames(p) = c("chr","pos","ref","A","C","G","T","af","cov","dbsnp")
      p = p[which(p$cov > 10 ),,drop=F]
      p$group = paste(p$chr,p$pos,p$ref,sep = ":")
      p = data.table(p)
      # merge
      p = merge(x = p,y = pbem.chr,all.x = T,sort = T,by = "group")
      p = data.frame(p,stringsAsFactors = F)
      # apply scaling
      out = mclapply(seq(1,nrow(p),1),scaling,p=p,cov.scaled = 0.10,mc.cores = mc.cores,mc.preschedule = T)
      p.scaled = fromListToDF(out)
      # write pileup file
      write.table(x = p.scaled,file = paste0(basename(y),".scaled"),col.names = F,row.names = F,sep = "\t",quote = F,na = "")
    }
    # cat together pileups
    setwd(file.path(out.folder,basename(id),"pileup"))
    cmd = paste("head -n 1",x,">",basename(x))
    system(cmd)
    cmd = paste("cat split.folder/*.scaled >>",basename(x))
    system(cmd)
    # remove split folder
    system("rm -r split.folder")
    # generate and write snvs file
    setwd(file.path(out.folder,basename(id),"snvs"))
    p.scaled = fread(input = file.path(out.folder,basename(id),"pileup",basename(x)),stringsAsFactors = F,showProgress = F,header = T,na.strings = "",colClasses = list(character=10),data.table = F)
    s.scaled = p.scaled[which(p.scaled$af > 0),,drop=F]
    s.scaled$Ars = 0 
    s.scaled$Crs = 0
    s.scaled$Grs = 0
    s.scaled$Trs = 0
    s.scaled = s.scaled[,c("chr","pos","ref","A","C","G","T","af","cov","Ars","Crs","Grs","Trs","dbsnp")]
    write.table(x = s.scaled,file = gsub(basename(x),pattern = "\\.pileup$",replacement = ".snvs"),col.names = T,row.names = F,sep = "\t",quote = F,na = "")
  }
}

