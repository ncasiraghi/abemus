#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=1){
  message("\nERROR:\t1 argument required")
  message("\nUSAGE:\tfolder with count_tabs.RData\n")
  quit()
}

library(data.table)
library(ggplot2)
library(reshape2)
library(gtable)
library(grid)
library(gridExtra)
library(plyr)

wd <- args[1]
setwd(wd)
rdata = list.files(path = wd,pattern = "count_tab",full.names = T)

# collect data for Strelka
rdata_strelka = grep(rdata,pattern = "\\.strelka.RData$",value = T)

fulldata_strelka=c()
df_abemus_strelka=c()
for(x in rdata_strelka){
  message(x)
  load(file = x)
  count_tab$R = as.numeric(gsub(basename(x),pattern = "count_tab|\\.strelka.RData",replacement = ""))
  fulldata_strelka=rbind(fulldata_strelka,count_tab)
  for(s in count_tab$case){
    
    if(length(as.numeric(pbem_only_strelka[[s]]))>0){
      df.strelka=data.frame(class="strelka",
                            sample=s,
                            pbem=as.numeric(pbem_only_strelka[[s]]),
                            af=as.numeric(af_only_strelka[[s]]),
                            stringsAsFactors = F)
    } else {
      df.strelka=data.frame(class="strelka",
                            sample=s,
                            pbem=NA,
                            af=NA,
                            stringsAsFactors = F)
    }
    
    if(length(as.numeric(pbem_only_abemus[[s]]))>0){
      df.abemus=data.frame(class="abemus",
                           sample=s,
                           pbem=as.numeric(pbem_only_abemus[[s]]),
                           af=as.numeric(af_only_abemus[[s]]),
                           stringsAsFactors = F)
    } else {
      df.abemus=data.frame(class="abemus",
                           sample=s,
                           pbem=NA,
                           af=NA,
                           stringsAsFactors = F)
    }
    
    if(length(as.numeric(pbem_intersect[[s]]))>0){
      df.intersect=data.frame(class="abemus AND strelka",
                              sample=s,
                              pbem=as.numeric(pbem_intersect[[s]]),
                              af=as.numeric(af_intersect[[s]]),
                              stringsAsFactors = F)
    } else {
      df.intersect=data.frame(class="abemus AND strelka",
                              sample=s,
                              pbem=NA,
                              af=NA,
                              stringsAsFactors = F)
    }
    df.tmp=rbind(df.strelka,df.abemus,df.intersect)
    df.tmp$R = as.numeric(gsub(basename(x),pattern = "count_tab|\\.strelka.RData",replacement = ""))
    df_abemus_strelka=rbind(df_abemus_strelka,df.tmp)
    rm(df.strelka,df.abemus,df.intersect)
  }
}


# collect data for Somatic Sniper
rdata_sniper = grep(rdata,pattern = "\\.sniper.RData$",value = T)

fulldata_sniper=c()
df_abemus_sniper=c()
for(x in rdata_sniper){
  message(x)
  load(file = x)
  count_tab$R = as.numeric(gsub(basename(x),pattern = "count_tab|\\.sniper.RData",replacement = ""))
  fulldata_sniper=rbind(fulldata_sniper,count_tab)
  for(s in count_tab$case){
    
    if(length(as.numeric(pbem_only_sniper[[s]]))>0){
      df.sniper=data.frame(class="sniper",
                           sample=s,
                           pbem=as.numeric(pbem_only_sniper[[s]]),
                           af=as.numeric(af_only_sniper[[s]]),
                           stringsAsFactors = F)
    } else {
      df.sniper=data.frame(class="sniper",
                           sample=s,
                           pbem=NA,
                           af=NA,
                           stringsAsFactors = F)
    }
    
    if(length(as.numeric(pbem_only_abemus[[s]]))>0){
      df.abemus=data.frame(class="abemus",
                           sample=s,
                           pbem=as.numeric(pbem_only_abemus[[s]]),
                           af=as.numeric(af_only_abemus[[s]]),
                           stringsAsFactors = F)
    } else {
      df.abemus=data.frame(class="abemus",
                           sample=s,
                           pbem=NA,
                           af=NA,
                           stringsAsFactors = F)
    }
    
    if(length(as.numeric(pbem_intersect[[s]]))>0){
      df.intersect=data.frame(class="abemus AND sniper",
                              sample=s,
                              pbem=as.numeric(pbem_intersect[[s]]),
                              af=as.numeric(af_intersect[[s]]),
                              stringsAsFactors = F)
    } else {
      df.intersect=data.frame(class="abemus AND sniper",
                              sample=s,
                              pbem=NA,
                              af=NA,
                              stringsAsFactors = F)
    }
    df.tmp=rbind(df.sniper,df.abemus,df.intersect)
    df.tmp$R = as.numeric(gsub(basename(x),pattern = "count_tab|\\.sniper.RData",replacement = ""))
    df_abemus_sniper=rbind(df_abemus_sniper,df.tmp)
    rm(df.sniper,df.abemus,df.intersect)
  }
}


# Merge data
fulldata = merge(x = fulldata_strelka,y = fulldata_sniper,by = c("patient","case","R","stm_pm","abemus","stm_pm_abemus","stm_pm_abemus_called"))
df = rbind(df_abemus_strelka,df_abemus_sniper)

# Add admixture levels to fulldata 
adms=c()
for(i in 1:nrow(fulldata)){
  a=unlist(strsplit(fulldata$case[i],split = "_adm"))
  if(length(a)==1){
    adms=c(adms,0)
  }
  else {
    a=gsub(a[2],pattern = "\\.sorted",replacement = "")
    a=as.numeric(gsub(a,pattern = "_",replacement = "."))
    adms=c(adms,a)
  }
}
fulldata$admixture<-adms

#------------------------------------------------------------------------
# Studio6 [ Fig.4 ]
#------------------------------------------------------------------------
pdf(file = "studio6.pdf",h = 219, w = 297, paper='A4r')

tmp=sapply(X = fulldata$stm_pm,FUN = strsplit,split = "\\|")
fulldata$pms_available=as.numeric(lapply(tmp,length))
pmtab_aggregate = aggregate(pms_available ~ R + admixture,data=fulldata, sum)
nsample_count = count(fulldata,c("R","admixture"))

tab=c()
for(R in unique(fulldata$R)){
  message(R)
  for(a in unique(fulldata$admixture)){
    message(a)
    this=fulldata[which(fulldata$R==R & fulldata$admixture==a),]
    out_abemus=data.frame(class=paste0("abemus_R",R),
                          level=a,
                          sum_overall=sum(this$abemus),
                          sum_stm=sum(this$stm_pm_abemus_called),
                          stringsAsFactors = F)
    tab=rbind(tab,out_abemus)
  }
}

for(a in unique(fulldata$admixture)){
  message(a)
  this=fulldata[which(fulldata$R==1 & fulldata$admixture==a),]
  out_strelka=data.frame(class="strelka",
                         level=a,
                         sum_overall=sum(this$strelka),
                         sum_stm=sum(this$stm_pm_strelka_called),
                         stringsAsFactors = F)
  tab=rbind(tab,out_strelka)
}

for(a in unique(fulldata$admixture)){
  message(a)
  this=fulldata[which(fulldata$R==1 & fulldata$admixture==a),]
  out_sniper=data.frame(class="sniper",
                        level=a,
                        sum_overall=sum(this$sniper),
                        sum_stm=sum(this$stm_pm_sniper_called),
                        stringsAsFactors = F)
  tab=rbind(tab,out_sniper)
}


pmstot = data.frame(level=unique(fulldata$admixture),npms=NA,nsamples=NA,stringsAsFactors = F)
for(i in 1:nrow(pmstot)){
  pmstot$npms[i]<-unique(pmtab_aggregate$pms_available[which(pmtab_aggregate$admixture==pmstot$level[i])])
  pmstot$nsamples[i]<-unique(nsample_count$freq[which(nsample_count$admixture==pmstot$level[i])])
}

bp1<-ggplot(pmstot, aes(factor(level),nsamples,label = nsamples)) + geom_bar(stat="identity") + geom_label() + theme(aspect.ratio=1)
bp2<-ggplot(pmstot, aes(factor(level),npms,label = npms)) + geom_bar(stat="identity") + geom_label() + theme(aspect.ratio=1)
grid.arrange(bp1,bp2,nrow=1)
grid.newpage()

tab=merge(x = tab,y = pmstot,by = "level",all.x = T)

tab$class = gsub(tab$class,pattern = "abemus_",replacement = "")

tab$score = tab$sum_stm/tab$npms

bpcovered=read.delim("/scratch/sharedCO/Casiraghi/Abemus_data/AmpliSeq_AR_STM2015/TargetPositions/bpcovered.tsv",as.is = T)
tot.bpcovered = as.numeric(bpcovered$bp_covered[which(bpcovered$chromosome=="TOTAL")])
tab$sum_overall_norm = tab$sum_overall/(tab$nsamples*tot.bpcovered)
tab$sum_stm_norm = tab$sum_stm/tab$sum_overall
tab$double_ratio = (tab$sum_stm_norm)*(tab$score)

tab = tab[which(!tab$level%in%c(0,20,60)),]

# p1<-ggplot(tab, aes(x=class, y=sum_overall)) +
#   geom_bar(stat="identity",fill="grey20") +
#   facet_grid(.~level) +
#   theme(axis.text.x=element_text(angle=90,hjust = 1,vjust = 0.5,size=4),strip.text.y = element_text(size = 4))

# q1<-ggplot(tab, aes(x=class, y=sum_overall_norm)) +
#   geom_bar(stat="identity",fill="grey20") +
#   facet_grid(.~level) +
#   theme(axis.text.x=element_text(angle=90,hjust = 1,vjust = 0.5,size=4),strip.text.y = element_text(size = 4))

colfunc <- colorRampPalette(c("#9ecae1", "#08519c"))
col_values=c(colfunc(11),"#ff7f00","#b2df8a")

p2<-ggplot(tab, aes(x=class, y=double_ratio,fill=class)) + scale_fill_manual(values = col_values) + 
  geom_bar(stat="identity") +
  facet_grid(.~level) +
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust = 0.5,size=4),strip.text.y = element_text(size = 4),legend.position="none")

q2<-ggplot(tab, aes(x=class, y=sum_stm_norm, fill=class)) + ylim(0,1) + scale_fill_manual(values = col_values) +
  geom_bar(stat="identity") +
  facet_grid(.~level) +
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust = 0.5,size=4),strip.text.y = element_text(size = 4),legend.position="none")

p3<-ggplot(tab, aes(x=class, y=score, fill=class)) + ylim(0,1) + scale_fill_manual(values = col_values) +
  geom_bar(stat="identity") +
  facet_grid(.~level) + 
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust = 0.5,size=4),strip.text.y = element_text(size = 4),legend.position="none")

g1 <- ggplotGrob(q2)
g2 <- ggplotGrob(p3)
g3 <- ggplotGrob(p2)
g <- rbind(g1, g2, g3, size = "first")
#g <- rbind(g1, g2, size = "first")
g$widths <- unit.pmax(g1$widths, g2$widths, g3$widths)
#g$widths <- unit.pmax(g1$widths, g2$widths)
grid.draw(g)

# p1<-ggplot(tab, aes(x=factor(level), y=sum_stm, group=class, color=class)) + 
#   geom_point(position=position_dodge(0.7),shape=19,size=4)+
#   scale_color_brewer(palette="Paired")+theme_minimal()
# 
# p2<-ggplot(tab, aes(x=factor(level), y=sum_overall, group=class, color=class)) + 
#   geom_point(position=position_dodge(0.7),shape=19,size=4)+
#   scale_color_brewer(palette="Paired")+theme_minimal()

dev.off()



# Studio 1
pdf(file = "studio1.pdf",h = 219, w = 297, paper='A4r')
col_values=c("#636363","#636363","#3182bd","#3182bd","#feb24c","#feb24c")
for(pt in unique(fulldata$patient)){
  message(pt)
  m=fulldata[which(fulldata$patient==pt),]
  h=rbind(data.frame(class="only abemus (vs strelka)",sample=m$case,count=m$only_abemus.x,R=m$R,stringsAsFactors = F),
          data.frame(class="only strelka (vs abemus)",sample=m$case,count=m$only_strelka,R=m$R,stringsAsFactors = F),
          data.frame(class="abemus AND strelka",sample=m$case,count=m$intersect_abemus_strelka,R=m$R,stringsAsFactors = F),
          data.frame(class="only abemus (vs sniper)",sample=m$case,count=m$only_abemus.y,R=m$R,stringsAsFactors = F),
          data.frame(class="only sniper (vs abemus)",sample=m$case,count=m$only_sniper,R=m$R,stringsAsFactors = F),
          data.frame(class="abemus AND sniper",sample=m$case,count=m$intersect_abemus_sniper,R=m$R,stringsAsFactors = F))
  gb=ggplot(h, aes(x=class, y=count, fill=class)) + scale_fill_manual(values = col_values) + geom_bar(stat="identity") + facet_grid(sample~R) + theme(legend.position="none") + ggtitle(pt) + theme(axis.text.x=element_text(angle=90,hjust = 1,vjust = 0.5),strip.text.y = element_text(size = 4))
  plot(gb)
}
dev.off()

# Studio 2a
# pdf(file = "studio2a.pdf",h = 219, w = 297, paper='A4r')
# # to be checked the pbem in sniper
# for(pt in unique(count_tab$patient)){
#   message(pt)
#   m=df[grep(df$sample,pattern = pt),]
#   gb=ggplot(m, aes(x=class, y=pbem, fill=R)) + geom_boxplot(fill="white", color="black",outlier.shape = NA) + facet_grid(sample~R) + theme(legend.position="none") + ggtitle(pt) + geom_jitter(shape=16, position=position_jitter(0.2)) + theme(axis.text.x=element_text(angle=90,hjust = 1,vjust = 0.5),strip.text.y = element_text(size = 4))
#   plot(gb)
# }
# dev.off()
# 
# # Studio 2b
# pdf(file = "studio2b.pdf",h = 219, w = 297, paper='A4r')
# for(pt in unique(count_tab$patient)){
#   message(pt)
#   m=df[grep(df$sample,pattern = pt),]
#   gb=ggplot(m, aes(x=class, y=af, fill=R)) + geom_boxplot(fill="white", color="black",outlier.shape = NA) + facet_grid(sample~R) + theme(legend.position="none") + ggtitle(pt) + geom_jitter(shape=16, position=position_jitter(0.2)) + theme(axis.text.x=element_text(angle=90,hjust = 1,vjust = 0.5),strip.text.y = element_text(size = 4))
#   plot(gb)
# }
# dev.off()

# Studio 3
pdf(file = "studio3.pdf",h = 219, w = 297, paper='A4r')
for(pt in unique(fulldata$patient)){
  message(pt)
  m=fulldata[which(fulldata$patient==pt),]
  m=m[with(m,order(case,decreasing = F)),]
  pms=unlist(strsplit(unique(m$stm_pm),split = "\\|"))
  for(j in 1:length(pms)){
    message(pms[j])
    x=data.frame(sample=unique(m$case),stringsAsFactors = F)
    for(i in sort(unique(m$R))){
      k=strsplit(m$stm_pm_abemus[which(m$R==i)],split = "\\|")
      abemus=as.numeric(unlist(lapply(k, `[[`, j)))
      z=strsplit(m$stm_pm_strelka[which(m$R==i)],split = "\\|")
      strelka=as.numeric(unlist(lapply(z, `[[`, j)))
      strelka[which(strelka==1)]<-3
      w=strsplit(m$stm_pm_sniper[which(m$R==i)],split = "\\|")
      sniper=as.numeric(unlist(lapply(w, `[[`, j)))
      sniper[which(sniper==1)]<-5
      x[,as.character(i)]=abemus+strelka+sniper
    }
    df_hm=melt(x, id.vars = "sample")
    names(df_hm)[2:3] <- c("R", "call")
    #colors <- c("0"= "white","1"="gold", "3"="red", "5"="blue","4"="orange","6"="darkolivegreen2","8"="darkmagenta")
    df_hm$call <- factor(df_hm$call,levels=c(0,1,3,5,4,6,8))
    gpl=ggplot(df_hm, aes(R, sample )) + geom_tile(aes(fill = call), color = "black") + coord_fixed(ratio=1) + ggtitle(pms[j]) + 
        scale_fill_manual(values = c("white","yellow","red","blue","orange","darkolivegreen2","darkmagenta"),
                          labels=c("none","only abemus","only strelka","only sniper","abemus AND strelka","abemus AND sniper","strelka AND sniper"),
                          drop=FALSE,name="who")
    plot(gpl)
  }
}
dev.off()







if(FALSE){
# Studio 4
#pdf(file = "studio4.pdf",h = 219, w = 297, paper='A4r')
tab=c()
for(R in unique(fulldata$R)){
  message(R)
  for(a in unique(fulldata$admixture)){
    message(a)
    this=fulldata[which(fulldata$R==R & fulldata$admixture==a),]
    out_abemus=data.frame(class=paste0("abemus_R",R),
                          level=a,
                          mean=mean(this$abemus),
                          median=median(this$abemus),
                          sd=sd(this$abemus),
                          mean_stm_called=mean(this$stm_pm_abemus_called),
                          median_stm_called=median(this$stm_pm_abemus_called),
                          stringsAsFactors = F)
    tab=rbind(tab,out_abemus)
  }
}

for(a in unique(fulldata$admixture)){
  message(a)
  this=fulldata[which(fulldata$R==1 & fulldata$admixture==a),]
  out_strelka=data.frame(class="strelka",
                         level=a,
                         mean=mean(this$strelka),
                         median=median(this$strelka),
                         sd=sd(this$strelka),
                         mean_stm_called=mean(this$stm_pm_strelka_called),
                         median_stm_called=median(this$stm_pm_strelka_called),
                         stringsAsFactors = F)
  tab=rbind(tab,out_strelka)
}

ggplot(tab, aes(x=factor(level), y=mean, group=class, color=class)) + 
  geom_point(position=position_dodge(0.7),shape=18,size=4)+
  scale_color_brewer(palette="Paired")+theme_minimal()

ggplot(tab, aes(x=factor(level), y=median, group=class, color=class)) + 
  geom_point(position=position_dodge(0.7),shape=18,size=4)+
  scale_color_brewer(palette="Paired")+theme_minimal()

ggplot(tab, aes(x=factor(level), y=median_stm_called/median, group=class, color=class)) + 
  geom_line()+geom_point()+
  scale_color_brewer(palette="Paired")+theme_minimal()

ggplot(tab, aes(x=factor(level), y=mean_stm_called/mean, group=class, color=class)) + 
  geom_line()+geom_point()+
  scale_color_brewer(palette="Paired")+theme_minimal()

# Studio 5
#pdf(file = "studio5.pdf",h = 219, w = 297, paper='A4r')
tab=c()
for(i in 1:nrow(fulldata)){
  this=fulldata[i,]
  out_abemus=data.frame(patient=this$patient,
                        case=this$case,
                        class=paste0("abemus_R",this$R),
                        level=this$admixture,
                        score=this$stm_pm_abemus_called/this$abemus,
                        stringsAsFactors = F)
  tab=rbind(tab,out_abemus)
}

reduced=fulldata[which(fulldata$R==1),]
for(i in 1:nrow(reduced)){
  this=reduced[i,]
  out_strelka=data.frame(patient=this$patient,
                         case=this$case,
                         class="strelka",
                         level=this$admixture,
                         score=this$stm_pm_strelka_called/this$strelka,
                         stringsAsFactors = F)
  tab=rbind(tab,out_strelka)
}

for(pt in unique(tab$patient)){
  message(pt)
  h=tab[which(tab$patient==pt),]
  gg<-ggplot(h, aes(x=factor(level), y=score, group=class, color=class))+geom_line()+geom_point()+scale_color_brewer(palette="Paired")+theme_minimal()+ggtitle(pt)
  #gg<-ggplot(h, aes(x=factor(level), y=score, group=class, color=class))+geom_point(position=position_dodge(0.7),shape=16,size=3)+scale_color_brewer(palette="Paired")+theme_minimal()+ggtitle(pt)
  plot(gg)
}

dev.off()

}
