library(data.table)

args <- commandArgs(trailingOnly = TRUE)

if(length(args)!=4){
  message("\nERROR! Required arguments:\n\n1. out.folder\n2. mutations_file\n3. full_coverage_folder\n4. coverage scaling factor\n")
  quit()
}

wd = args[1]
mutation_file = args[2]
full_coverage_folder = args[3]
sf = as.numeric(args[4])

#wd = "/scratch/sharedCO/Casiraghi/InSilicoData_Tumors/HALO_t0.1_byChrom_c10_model"
#mutation_file="/scratch/sharedCO/Casiraghi/InSilicoData_Tumors/HALO_t0.1_byChrom_100genes.tsv"
#full_coverage_folder="/scratch/sharedCO/Casiraghi/InSilicoData_Tumors/HALO_t0.1_byChrom"
#sf = 0.10

setwd(wd)
admixtures = grep(basename(list.dirs(full_coverage_folder,recursive = F,full.names = T)),pattern = "adm",value = T)
mutations = read.delim(mutation_file,header = F,stringsAsFactors = F)
mutations = mutations[,1:2]
colnames(mutations)=c("chr","pos")
mutations <- within(mutations,  group <- paste(chr, pos, sep=":"))

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

# Scaling pileup
for(x in admixtures){
  message(x)
  # folders output
  ps = list.files(path = file.path(full_coverage_folder,x,"pileup"),full.names = T)
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
      file.to.correct = gsub(p,pattern = basename(full_coverage_folder),replacement = basename(wd))
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
