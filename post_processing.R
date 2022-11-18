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
  fn <- gh('graphs/allscout_gm.{cst}.Rdata')
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
allm[arm=='idh',ARM:='DH-focussed']
allm[arm=='iph',ARM:='PHC-focussed']
allm[,ARM:=factor(ARM,levels = c('SOC','DH-focussed','PHC-focussed'),ordered = TRUE)]


## load data by activity from CSV
all2 <- list()
for(cst in cstrata){
  fn <- gh('graphs/allout_gm.{cst}.csv')
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
allm2[variable=='idh',ARM:='DH-focussed']
allm2[variable=='iph',ARM:='PHC-focussed']
allm2[,ARM:=factor(ARM,levels = c('SOC','DH-focussed','PHC-focussed'),ordered = TRUE)]


## make graphs
## first by stage/level/arm
gd1 <- allm[,.(value=sum(value)),by=.(iso3,`stage/level`,ARM)]

## absolute
GP <- ggplot(gd1,aes(iso3,value,fill=`stage/level`))+
  geom_bar(stat = 'identity')+
  facet_wrap(~ARM)+
  xlab('Country')+ylab('Cost per child treated (USD)')+
  rot45

ggsave(GP,file=gh('graphs/pp_stagecosts_full.png'),w=7,h=5)


## proportional
GP <- ggplot(gd1,aes(iso3,value,fill=`stage/level`))+
  geom_bar(stat = 'identity',position='fill')+
  scale_y_continuous(label=percent)+
  facet_wrap(~ARM)+
  xlab('Country')+ylab('Proportion of cost per child treated (USD)')+
  rot45

ggsave(GP,file=gh('graphs/pp_stagecosts_pc.png'),w=7,h=5)


## versions that work by stage
## gd2 <- allm[,.(value=sum(value)),by=.(iso3,ARM,activity=tolower(stratum))]
gd2 <- allm2
gd2[,activity:=tolower(stratum)]

## absolute
GP <- ggplot(gd2,aes(iso3,value,fill=activity))+
  geom_bar(stat = 'identity')+
  facet_wrap(~ARM)+
  xlab('Country')+ylab('Cost per child treated (USD)')+
  rot45

ggsave(GP,file=gh('graphs/pp_actcosts_full.png'),w=7,h=5)


## proportional
GP <- ggplot(gd2,aes(iso3,value,fill=activity))+
  geom_bar(stat = 'identity',position='fill')+
  scale_y_continuous(label=percent)+
  facet_wrap(~ARM)+
  xlab('Country')+ylab('Proportion of cost per child treated (USD)')+
  rot45

ggsave(GP,file=gh('graphs/pp_actcosts_pc.png'),w=7,h=5)


## check
gd1[,sum(value),by=.(iso3,ARM)]
gd2[,sum(value),by=.(iso3,ARM)]
