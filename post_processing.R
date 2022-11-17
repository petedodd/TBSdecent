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

## make graphs

## absolute
GP <- ggplot(allm,aes(iso3,value,fill=`stage/level`))+
  geom_bar(stat = 'identity')+
  facet_wrap(~ARM)+
  xlab('Country')+ylab('Cost per child treated (USD)')+
  rot45

ggsave(GP,file=gh('graphs/pp_stagecosts_full.png'),w=7,h=5)


## proportional
GP <- ggplot(allm,aes(iso3,value,fill=`stage/level`))+
  geom_bar(stat = 'identity',position='fill')+
  scale_y_continuous(label=percent)+
  facet_wrap(~ARM)+
  xlab('Country')+ylab('Proportion of cost per child treated (USD)')+
  rot45

ggsave(GP,file=gh('graphs/pp_stagecosts_pc.png'),w=7,h=5)
