library(here)
library(glue)
library(data.table)
library(ggplot2)
library(scales)
rot45 <- theme(axis.text.x = element_text(angle = 45, hjust = 1))
gh <- function(x)glue(here(x))

cstrata <- c('DECENT_EQUIPMENT',
             'DECENT_PERSONNEL',
             'DECENT_SUPPLY',
             'DECENT_TRAINING')

## load data
all <- list()
for(cst in cstrata){
  fn <- gh('graphs/allscout_iphbased.{cst}.Rdata')
  load(fn)
  allscout[,stratum:=gsub('DECENT_','',cst)]
  all[[cst]] <- allscout
}
all <- rbindlist(all)

## drop presumed,presented
nmz <- names(all)
drop <- grep('pres',nmz,value=TRUE)
all[,(drop):=NULL]
allm <- melt(all,id=c('iso3','stratum'))
allm[,variable:=gsub('perATT\\.','',variable)]
allm[,variable:=gsub('cost\\.','',variable)]
allm[,variable:=gsub('\\.mid','',variable)]
allm[,c('clinic level','stage','arm'):=tstrsplit(variable,split='\\.')]
allm[,`stage/level`:=paste(stage,`clinic level`,sep='/')]
allm[,ARM:='SOC']
allm[arm=='idh',ARM:='DH-focused']
allm[arm=='iph',ARM:='PHC-focused']
allm[,ARM:=factor(ARM,levels = c('SOC','DH-focused','PHC-focused'),ordered = TRUE)]


## load data by activity from CSV
all2 <- list()
for(cst in cstrata){
  fn <- gh('graphs/allout_iphbased.{cst}.csv')
  tmp <- fread(fn)
  tmp[,stratum:=gsub('DECENT_','',cst)]
  all2[[cst]] <- tmp
}
all2 <- rbindlist(all2)

## drop presumed,presented
all2 <- all2[,.(iso3,stratum,costperATT.soc.mid,costperATT.iph.mid, costperATT.idh.mid)]
allm2 <- melt(all2,id=c('iso3','stratum'))
allm2[,variable:=gsub('costperATT\\.','',variable)]
allm2[,variable:=gsub('\\.mid','',variable)]
allm2[,ARM:='SOC']
allm2[variable=='idh',ARM:='DH-focused']
allm2[variable=='iph',ARM:='PHC-focused']
allm2[,ARM:=factor(ARM,levels = c('SOC','DH-focused','PHC-focused'),ordered = TRUE)]


## make graphs
## first by stage/level/arm
gd1 <- allm[,.(value=sum(value)),by=.(iso3,`stage/level`,ARM)]
gd1[,ARM:=gsub('ssed','sed',ARM)]
gd1$ARM <- factor(gd1$ARM,levels=unique(gd1$ARM),ordered = TRUE)
gd1 <- gd1[iso3 != 'ZMB']
gd1[,lbl:=paste0('US$',round(value))]
gdt <- gd1[,.(value=sum(value)),by=.(iso3,ARM)]
gdt[,lbl:=paste0('US$',round(value))]
gdt[,`stage/level`:=NA]
gd1[,tot:=sum(value),by=.(iso3,ARM)]
gd1[,pclbl:=paste0(round(1e2*value/tot,1),'%')]
gd1[,pc:=value/tot]

## absolute
GP <- ggplot(gd1,aes(iso3,value,fill=`stage/level`,label=lbl))+
  geom_bar(stat = 'identity')+
  facet_wrap(~ARM)+
  xlab('Country')+ylab('Cost per child treated (US$)')+
  geom_text(position = position_stack(vjust = 0.5))+
  geom_text(data=gdt,vjust = -0.2,fontface='bold')+
  ggthemes::scale_fill_few()+
  theme_linedraw()
GP

ggsave(GP,file=gh('graphs/pp_stagecosts_full.png'),w=15,h=10)


## proportional
GP <- ggplot(gd1,aes(iso3,value,fill=`stage/level`))+
  geom_bar(stat = 'identity',position='fill')+
  scale_y_continuous(label=percent)+
  facet_wrap(~ARM)+
  ggthemes::scale_fill_few()+
  geom_text(aes(x=iso3,y=pc,label=pclbl),position = position_stack(vjust = 0.5))+
  xlab('Country')+ylab('Proportion of cost per child treated (US$)')+
  theme_linedraw()
GP

ggsave(GP,file=gh('graphs/pp_stagecosts_pc.png'),w=15,h=10)


## versions that work by stage
## gd2 <- allm[,.(value=sum(value)),by=.(iso3,ARM,activity=tolower(stratum))]
gd2 <- allm2
gd2[,`cost category`:=tolower(stratum)]
gd2[,ARM:=gsub('ssed','sed',ARM)]
gd2$ARM <- factor(gd2$ARM,levels=unique(gd2$ARM),ordered = TRUE)
gd2 <- gd2[iso3 != 'ZMB']
gd2[,lbl:=paste0('US$',round(value))]
gd2t <- gd2[,.(value=sum(value)),by=.(iso3,ARM)]
gd2t[,lbl:=paste0('US$',round(value))]
gd2t[,`cost category`:=NA]
gd2[,tot:=sum(value),by=.(iso3,ARM)]
gd2[,pclbl:=paste0(round(1e2*value/tot,1),'%')]
gd2[,pc:=value/tot]


## absolute
GP <- ggplot(gd2,aes(iso3,value,fill=`cost category`,label=lbl))+
  geom_bar(stat = 'identity')+
  facet_wrap(~ARM)+
  xlab('Country')+ylab('Cost per child treated (US$)')+
  geom_text(position = position_stack(vjust = 0.5))+
  geom_text(data=gd2t,vjust = -0.2,fontface='bold')+
  ggthemes::scale_fill_few()+
  theme_linedraw()
GP

ggsave(GP,file=gh('graphs/pp_actcosts_full.png'),w=15,h=10)


## proportional
GP <- ggplot(gd2,aes(iso3,value,fill=`cost category`))+
  geom_bar(stat = 'identity',position='fill')+
  scale_y_continuous(label=percent)+
  facet_wrap(~ARM)+
  xlab('Country')+ylab('Proportion of cost per child treated (US$)')+
  ggthemes::scale_fill_few()+
  geom_text(aes(x=iso3,y=pc,label=pclbl),position = position_stack(vjust = 0.5))+
  xlab('Country')+ylab('Proportion of cost per child treated (US$)')+
  theme_linedraw()
GP

ggsave(GP,file=gh('graphs/pp_actcosts_pc.png'),w=15,h=10)


## check
gd1[,sum(value),by=.(iso3,ARM)]
gd2[,sum(value),by=.(iso3,ARM)]


## new cascadegraph
load(file=gh('graphs/CT_iphbased.Rdata'))

CT <- CT[,.(value=sum(value)),by=.(id,arm,stage,age)]
CTT <- CT[,.(value=sum(value)),by=.(id,arm,stage)]
CTT[,age:='0-14']
CT <- rbind(CT,CTT)
CT$stage <- factor(CT$stage,levels=c('presented','screened','presumed','treated'),ordered = TRUE)
CT$age <- factor(CT$age,levels=c('0-14','0-4','5-14'),ordered = TRUE)
CT[,stagel:='Presented at health facility']
CT[stage=='screened',stagel:='Screened for tuberculosis']
CT[stage=='presumed',stagel:='Assessed as presumptive tuberculosis']
CT[stage=='treated',stagel:='Treated for tuberculosis']
CT$stagel <- factor(CT$stagel,levels=c('Presented at health facility',
                                       'Screened for tuberculosis',
                                       'Assessed as presumptive tuberculosis',
                                       'Treated for tuberculosis'),ordered = TRUE)
CT[,agel:=paste0(age,' years')]
CT$agel <- factor(CT$agel,levels=c('0-14 years','0-4 years','5-14 years'),ordered = TRUE)
CT[,arml:='Standard of care']
CT[arm=='idh',arml:='DH-focused intervention']
CT[arm=='iph',arml:='PH-focused intervention']
CT$arml <- factor(CT$arml,levels=c('Standard of care',
                                   'DH-focused intervention',
                                   'PH-focused intervention'),ordered = TRUE)

## text summaries of step-downs
CTS0 <- CT[,.(value=mean(value)),by=.(arml,stage,agel)]
CTT <- dcast(CTS0,arml+agel~stage,value.var = 'value')
CTT[,treated:=format(round(1e2*treated/presumed,1),nsmall = 1)]
CTT[,presumed:=format(round(1e2*presumed/screened,1),nsmall = 1)]
CTT[,screened:=format(round(1e2*screened/presented,1),nsmall = 1)]
CTT[,presented:='']
CTT <- melt(CTT,id=c('arml','agel'))
CTT[,stagel:='Presented at health facility']
CTT[variable=='screened',stagel:='Screened for tuberculosis']
CTT[variable=='presumed',stagel:='Assessed as presumptive tuberculosis']
CTT[variable=='treated',stagel:='Treated for tuberculosis']

## final dataset
CTS <- CT[,.(value=mean(value)),by=.(arml,stagel,agel)]
CTS <- merge(CTS,CTT[,.(arml,stagel,agel,txt=paste0(value,'%'))],by=c('arml','stagel','agel'))
CTS[txt=='NA%',txt:='']
CTS$stagel <- factor(CTS$stagel,levels=c('Presented at health facility',
                                       'Screened for tuberculosis',
                                       'Assessed as presumptive tuberculosis',
                                       'Treated for tuberculosis'),ordered = TRUE)
CTS$arml <- factor(CTS$arml,levels=c('Standard of care',
                                     'DH-focused intervention',
                                     'PH-focused intervention'),ordered = TRUE)
CTS$agel <- factor(CTS$agel,levels=c('0-14 years','0-4 years','5-14 years'),ordered = TRUE)

w <- 0.8
dg <- position_dodge(w)
GP <- ggplot(CTS,aes(stagel,value,group=arml,fill=arml,label=txt))+
  facet_wrap(~agel)+
  geom_bar(stat='identity',position=dg,width = w)+
  ggthemes::scale_fill_few()+
  scale_y_sqrt(label=percent)+
  geom_text(position=dg,vjust=-0.25,size=2.5)+
  theme_linedraw()+
  ylab('Proportion of all children presenting (square root scale)')+
  xlab('Cascade stage')+
  theme(legend.position = 'top') + labs(fill=NULL)+rot45
GP

fn <- gh('graphs/model_cascade.png')
ggsave(GP,file=fn,h=7,w=14)
fn <- gh('graphs/model_cascade.pdf')
ggsave(GP,file=fn,h=7,w=14)

chk <- CT[,.(value=mean(value)),by=.(arm,stage,age)]
chk <- chk[age=='0-14']
chk <- dcast(chk,arm~stage,value.var = 'value')
chk[arm=='iph',treated]/chk[arm=='soc',treated]
chk[arm=='idh',treated]/chk[arm=='soc',treated]
chk[arm=='idh',treated]/chk[arm=='iph',treated]

## simpler version
w <- 0.8
dg <- position_dodge(w)
GP <- ggplot(CTS[agel=='0-14 years'],
             aes(stagel,value,group=arml,fill=arml,label=txt))+
  geom_bar(stat='identity',position=dg,width = w)+
  ggthemes::scale_fill_few()+
  scale_y_sqrt(label=percent)+
  geom_text(position=dg,vjust=-0.25,size=2.5)+
  theme_linedraw()+
  ylab('Proportion of all children presenting (square root scale)')+
  xlab('Cascade stage')+
  theme(legend.position = 'top') + labs(fill=NULL)+rot45
GP

fn <- gh('graphs/model_cascade_simple.png')
ggsave(GP,file=fn,h=7,w=7)
fn <- gh('graphs/model_cascade_simple.pdf')
ggsave(GP,file=fn,h=7,w=7)



## grabbing sensitivity analysis results

## prev = c(500,200,50)
## disr = c(0,0.03,0.06)
## DH, PH vs iso3

prev <- c('0.005','0.002','5e-04')
disr <- c('0.00','0.03','0.06')
ans <- list()
for(pv in prev){
  for(dr in disr){
    fn <- gh('graphs/{pv}{dr}tab1_iphbased.DECENT.csv')
    tmp <- fread(fn)
    tmp <- tmp[variable %in% c('ICER.idh','ICER.iph')]
    tmp[,prev:=pv]
    tmp[,discount.rate:=dr]
    ans[[paste(pv,dr)]] <- tmp
  }
}
ans <- rbindlist(ans)
ans <- ans[order(variable,prev,discount.rate)]
fwrite(ans,file=gh('graphs/sensitivity_tab1_iphbased.DECENT.csv'))
