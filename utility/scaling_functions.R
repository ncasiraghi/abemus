# based on AF

scaling <- function(i,p,cov.scaled,pbem.chr){
  message(i)
  pos = p[i,,drop=F]
  group = paste(pos$chr,pos$pos,pos$ref,sep = ":")
  # af = 0
  if(pos$af == 0){
    # apply scaling factor on coverage
    pos$cov = round(pos$cov*cov.scaled) 
    pos[,pos$ref] = pos$cov
  }
  # af > 0
  if(pos$af > 0){
    pos.pbem = pbem.chr[which(pbem.chr$group==group),,drop=F]
    if(nrow(pos.pbem)>0){
      # apply scaling factor on coverage
      pos$cov = round(pos$cov*cov.scaled) 
      # how many reads with error
      n.of.success = rbinom(size = pos$cov,n = 1,prob = pos.pbem$bperr)
      # adjust coverage of the reference
      pos[,pos$ref] = pos$cov - n.of.success
      alts = setdiff(c("A","C","G","T"),pos$ref)
      pos[,alts] = 0
      # how many errors for each not-reference allele
      if(n.of.success > 0){
        outs = sample(x = alts,size = n.of.success,replace = T,prob = pos.pbem[,paste0("total.",alts)]/pos.pbem$tot_coverage)
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
  }
  return(pos)
} 


# all positions

scaling <- function(i,p,cov.scaled,pbem.chr){
  message(i)
  pos = p[i,,drop=F]
  group = paste(pos$chr,pos$pos,pos$ref,sep = ":")
  pos.pbem = pbem.chr[which(pbem.chr$group==group),,drop=F]
  if(nrow(pos.pbem)>0){
    # apply scaling factor on coverage
    pos$cov = round(pos$cov*cov.scaled) 
    # how many reads with error
    n.of.success = rbinom(size = pos$cov,n = 1,prob = pos.pbem$bperr)
    # adjust coverage of the reference
    pos[,pos$ref] = pos$cov - n.of.success
    alts = setdiff(c("A","C","G","T"),pos$ref)
    pos[,alts] = 0
    # how many errors for each not-reference allele
    if(n.of.success > 0){
      outs = sample(x = alts,size = n.of.success,replace = T,prob = pos.pbem[,paste0("total.",alts)]/pos.pbem$tot_coverage)
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
  return(pos)
} 

# only for coverages lower than 10
