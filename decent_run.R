## top level handler for DECENTRALIZATION
rm(list=ls())
library(here)

## load other scripts
source(here('R/decent_tree.R'))           #tree structure and namings: also tree functions & libraries
source(here('R/decent_functions.R'))      #functions for tree parameters

DD <- fread(here('indata/DD.csv')) #cascade data for plots
DD$location <- toupper(DD$location)

## number of reps
nreps <- 1e3

## attributes to use
tblevels <- c('TB+','TB-','noTB') #bac confirmable TB, bac unconfirmable TB, not TB
hivlevels <- c(0,1)
artlevels <- c(0,1)
agelevels <- c('0-4','5-14')

## prior parameters
PD0 <- read.csv(here('indata/DecentParms - distributions.csv')) #read in

## temporary work with new parameters from data
PD1 <- read.csv(here('indata/PMZ.csv')) #read in
both <- intersect(PD0$NAME,PD1$NAME)
PD0 <- PD0[!PD0$NAME %in% both,] #supercede with PD1 version


## combine with above
P <- parse.parmtable(PD0)             #convert into parameter object
P2 <- parse.parmtable(PD1)             #convert into parameter object

P2$d.TBprev.ICS.o5 <- P2$d.TBprev.ICS.o5*1/5 #TODO expt
P2$d.TBprev.ICS.u5 <- P2$d.TBprev.ICS.u5*1/5
## P2$d.OR.dh.if.TB <- 8

P <- c(P,P2)


## ## version making test plots
## P <- parse.parmtable(PD0,
##                      testdir = here('graphs/test'),
##                      outfile = here('graphs/test/00out.csv'))

names(P)

## make base PSA dataset
set.seed(1234) #random number seed
D <- makePSA(nreps,P,dbls = list(c('cfrhivor','cfrartor')))

## ## NOTE temporary introduction of noise:
## for(nm in setdiff(PD1$NAME,'d.OR.dh.if.TB')){ #loop over probs
##   D[[nm]] <- ilogit(logit(D[[nm]]) + rnorm(nreps)/5)
## }
## D[['d.OR.dh.if.TB']] <- exp(log(D[['d.OR.dh.if.TB']]) + rnorm(nreps)/5)


## use these parameters to construct intput data by attribute
D <- makeAttributes(D)
D[,sum(value),by=id] #CHECK
D[tb!='noTB',sum(value),by=id] #CHECK
D[,sum(value),by=.(id,age)] #CHECK

## compute other parameters (adds by side-effect)
MakeTreeParms(D)

## check for leaks
## head(IPD.F$checkfun(D)) #IPD arm NOTE omitted from current country work
head(IPH.F$checkfun(D)) #IPH arm
head(IDH.F$checkfun(D)) #IDH arm
head(SOC.F$checkfun(D)) #SOC arm

## add cost data
csts <- fread(here('indata/testcosts.csv'))         #read cost data
C <- MakeCostData(csts,nreps)               #make cost PSA
D <- merge(D,C,by='id',all.x = TRUE)        #merge into PSA

names(D)
## now split tree into intervention and not
## D <- runallfuns(D,arm='all')                      #appends anwers

notIPD <- c('SOC','IDH','IPH')
D <- runallfuns(D,arm=notIPD)                      #appends anwers

## TODO need a function to compute weights by id


## --- cascade variables checking
nmz <- names(D)
rnmz0 <- grep("DH\\.|PHC\\.",nmz,value=TRUE)
rnmz <- c('id','age','tb','value',rnmz0)
A <- D[,..rnmz]
A[,sum(value),by=id] #CHECK
A[,pop:=value]       #rename for melting
A[,value:=NULL]
A[,TB:=ifelse(tb!='noTB','TB','not TB')] #simpler version of TB indicator
A[,tb:=NULL]                             #drop

nrow(A) #240K (attributes x 1000)
A[,sum(pop),by=id] #CHECK
A[TB=='TB',1e2*sum(pop),by=id] #CHECK
A[,sum(pop),by=.(id,age)] #CHECK


## population scaling
A[,c(rnmz0):=lapply(.SD,function(x) x*pop),
  .SDcols=rnmz0] #multiply variables by population
A[,pop:=NULL]                                             #can drop now

## melt
AM <- melt(A,id=c('id','age','TB'))
AM[,c('location','stage','arm'):=tstrsplit(variable,split='\\.')]
## sum over other attributes
AM <- AM[,.(value=sum(value)),by=.(id,arm,age,location,stage,TB)]
nrow(AM)
## aggregate/average
## TB version
AS <- AM[,.(mid=mean(value)),by=.(arm,age,TB,location,stage)]

## no TB version
AS2 <- AM[,.(mid=sum(value)),by=.(id,arm,age,location,stage)]
AS2 <- AS2[,.(mid=mean(mid)),by=.(arm,age,location,stage)]

## noTB graph version
tpl <- AS2[stage=='presented']
AS2 <- merge(AS2,tpl[,.(arm,age,location,pmid=mid)],
             by=c('arm','age','location'),all.x=TRUE)
AS2[,vpl:=1e5*mid/pmid]
lvls <- c('presented','screened','presumed','treated')
AS2$stage <- factor(AS2$stage,levels=lvls,ordered=TRUE)
AS2[,txt:=round(vpl)]
AS2[stage=='presented',txt:=NA]


## cascade plot
GP2 <- ggplot(AS2,aes(stage,vpl))+
  geom_bar(stat='identity')+
  facet_grid(location+age ~ arm)+
  scale_y_log10(label=comma)+
  ylab('Number per 100K attending')+
  xlab('Stage')+
  theme_light() + rot45+
  geom_text(aes(label=txt),col='red',vjust=0.1) +
  geom_point(data=DD,col='cyan',size=2)+
  geom_text(data=DD,aes(label=txt),col='cyan',vjust=2)
GP2

ggsave(GP2,file=here('graphs/cascade_plt.png'),w=15,h=15)

## TB graph version
AS <- merge(AS,tpl[,.(arm,age,location,pmid=mid)],
             by=c('arm','age','location'),all.x=TRUE)
AS[,vpl:=1e5*mid/pmid]
AS$stage <- factor(AS$stage,levels=lvls,ordered=TRUE)


GP3 <- ggplot(AS,aes(stage,vpl,fill=TB))+
  geom_bar(stat='identity',position='fill')+
  facet_grid(location+age ~ arm)+
  scale_y_sqrt(label=percent)+
  ylab('Number per 100K attending')+
  xlab('Stage')+
  theme_light() + rot45
## GP3

ggsave(GP3,file=here('graphs/cascade_plt_TB.png'),w=15,h=15)

## Treated true TB per 100K presented by arm
TTB <- AS[stage=='treated' & TB=='TB',.(TTBpl=1e5*sum(mid)),by=arm]
fwrite(TTB,file=here('graphs/TTB.csv'))

LYSdone <- TRUE
if(!LYSdone){
  ## make discounted life-years if they haven't been done
  isoz <- c('CIV','CMR','UGA','MOZ','ZMB','KHM') #relevant countries?
  LYK <- GetLifeYears(isolist=isoz,discount.rate=0.03,yearfrom=2021)
  save(LYK,file=here('indata/LYK.Rdata'))
} else {load(file=here('indata/LYK.Rdata'))}
LYK <- LYK[,.(LYS=mean(LYS)),by=age] #averaged life-years 4 generic tests

## generate some CEA outputs in graphs/ & outdata/
## NOTE these folders need to be created
## NOTE need ggpubr, BCEA installed
MakeCEAoutputs(D, #PSA dataset
               LYK, #discounted expected life-years by age
               file.id='test', #string to identify output files 
               Kmax=5e3,wtp=5e3)
