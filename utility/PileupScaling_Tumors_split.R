setwd("/scratch/sharedCO/Casiraghi/InSilicoData_Tumors/HALO_byChrom_c10_model")
library(data.table)

admixtures = grep(basename(list.dirs("/scratch/sharedCO/Casiraghi/InSilicoData_Tumors/HALO_byChrom",recursive = F,full.names = T)),pattern = "adm",value = T)
mutations = read.delim("/scratch/sharedCO/Casiraghi/Abemus_data/InSilicoData_Tumors/HALO_c10_model/Scripts/HALO_byChrom_200genes.tsv",header = F,stringsAsFactors = F)
mutations = mutations[,1:2]
colnames(mutations)=c("chr","pos")
mutations <- within(mutations,  group <- paste(chr, pos, sep=":"))

sf = 0.10

scaling <- function(x,sf){
  for(i in 1:nrow(x)){
    x[i,"A"] = round(x[i,"A"]*sf)
    x[i,"C"] = round(x[i,"C"]*sf)
    x[i,"G"] = round(x[i,"G"]*sf)
    x[i,"T"] = round(x[i,"T"]*sf)
    x[i,"cov"] = sum(x[i,c("A","C","G","T")])
    alt.cov = sum(x[i,setdiff(c("A","C","G","T"),x$ref[i])])
    x[i,"af"] = alt.cov/x[i,"cov"]
  }
  return(x)
}

# Scaling in pileup
admixtures=setdiff(admixtures,"adm20_0")

for(x in admixtures){
  message(x)
  ps = list.files(path = file.path("/scratch/sharedCO/Casiraghi/InSilicoData_Tumors/HALO_byChrom/",x,"pileup"),full.names = T)
  #p=ps[15]
  for(p in ps){
    message(p)
    # upload original 100% coverage file
    a=fread(input = p,stringsAsFactors = F,data.table = F,showProgress = F)
    ga=paste(a$chr, a$pos, sep=":")
    idx_a = which(ga %in% mutations$group)
    # scaling
    to_scale=a[idx_a,,drop=F]
    if(nrow(to_scale)>0){
      xsc = scaling(x = to_scale,sf = sf)
      # upload the scaled file
      file.to.correct = gsub(p,pattern = "HALO_byChrom",replacement = "HALO_byChrom_c10_model")
      b=fread(input = file.to.correct,stringsAsFactors = F,data.table = F,showProgress = F)
      gb=paste(b$chr, b$pos, sep=":")
      idx_b = which(gb %in% mutations$group)
      b=b[-idx_b,]
      # replace
      b=rbind(b,xsc)
      b=b[with(b,order(pos,decreasing = F)),]
      # write
      write.table(x = b,file = file.to.correct,col.names = T,row.names = F,sep = "\t",quote = F,na = "")
    } else {
      next()
    }
  }
}
