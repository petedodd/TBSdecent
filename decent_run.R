## top level handler for DECENTRALIZATION
rm(list=ls())
library(here)

## load other scripts
source(here('R/decent_tree.R'))           #tree structure and namings: also tree functions & libraries
source(here('R/decent_functions.R'))      #functions for tree parameters

## cascade data
DD <- fread(here('indata/DD.csv')) #cascade data for plots
DD$location <- toupper(DD$location)
DDW <- dcast(DD,age~location+stage+arm,value.var='vpl')
ICS <- fread(here('indata/ICS.csv')) #initial care-seeking props

## number of reps
nreps <- 1e3
set.seed(1234)

## attributes to use
tblevels <- c('TB+','TB-','noTB') #bac confirmable TB, bac unconfirmable TB, not TB
hivlevels <- c(0,1)
artlevels <- c(0,1)
agelevels <- c('0-4','5-14')
isoz <- c('KHM','CMR','CIV','MOZ','SLE','UGA','ZMB') #relevant countries


## --- life years and other outputs NOTE needs to be set FALSE on first run thru
LYSdone <- TRUE
if(!LYSdone){
  ## make discounted life-years if they haven't been done
  LYKc <- GetLifeYears(isolist=isoz,discount.rate=0.03,yearfrom=2021)
  LYK <- LYKc[,.(LYS=mean(LYS)),by=.(age)] #averaged life-years 4 generic tests
  save(LYKc,file=here('indata/LYKc.Rdata'))
  save(LYK,file=here('indata/LYK.Rdata'))
} else {
  load(file=here('indata/LYKc.Rdata'))
  load(file=here('indata/LYK.Rdata'))
}


## read and make cost data
csts <- fread(here('indata/testcosts.csv'))         #read cost data
rcsts <- fread(here('indata/TBS.DECENT.costs.csv'),skip = 1)    #read cost data
## check
setdiff(unique(rcsts$NAME),
        unique(csts$cost))
setdiff(unique(csts$cost),
        unique(rcsts$NAME))
allcosts <- reformatCosts(rcsts)
allcosts[cost.sd==0,cost.sd:=cost.m/40]        #SD such that 95% UI ~ 10% of mean
C <- MakeCostData(allcosts[iso3=='CIV'],nreps)               #make cost PSA NOTE using CIV cost data

## prior parameters
PD0 <- read.csv(here('indata/DecentParms - distributions.csv')) #read in
## parameters to be determined from cascade data
PD1 <- PD0[PD0$DISTRIBUTION=="",]
## the rest
PD0 <- PD0[PD0$DISTRIBUTION!="",]

## this computes and saves out the average accuracy of dx cascades
if(!file.exists(here('graphs'))) dir.create(here('graphs'))
DxA <- computeDxAccuracy(PD0,PD1,C,nreps)

## this computes and saves model parameters derived from cascade data
prevapproach <- 'gm'
PD1 <- computeCascadeParameters(DD,ICS,DxA,prevapproach)

## combine different parameter types
P1 <- parse.parmtable(PD0)             #convert into parameter object
P2 <- parse.parmtable(PD1)             #convert into parameter object
P <- c(P1,P2)
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
MakeTreeParms(D,P)

## check for leaks
## head(IPD.F$checkfun(D)) #IPD arm NOTE omitted from current country work
head(IPH.F$checkfun(D)) #IPH arm
head(IDH.F$checkfun(D)) #IDH arm
head(SOC.F$checkfun(D)) #SOC arm
names(SOC.F)

## add cost data
D <- merge(D,C,by='id',all.x = TRUE)        #merge into PSA

## run model
notIPD <- c('SOC','IDH','IPH')
D <- runallfuns(D,arm=notIPD)                      #appends anwers

## --- cascade variables checking
CO <- computeCascadeData(D)
AS2 <- CO$woTB
AS <- CO$wTB

## cascade plot
GP2 <- ggplot(AS2,aes(stage,vpl))+
  geom_bar(stat='identity')+
  facet_grid(location+age ~ arm)+
  scale_y_log10(label=comma)+
  ylab('Number per 100K attending')+
  xlab('Stage')+
  theme_light() + rot45+
  geom_text(aes(label=txt),col='red',vjust=0.1,hjust=-0.5) +
  geom_point(data=DD,col='cyan',size=2)+
  geom_text(data=DD,aes(label=txt),col='cyan',vjust=2,hjust=1)
GP2

ggsave(GP2,file=gh('graphs/cascade_plt_{prevapproach}.png'),w=15,h=15)

## Treated true TB per 100K presented by arm
TTB <- AS[stage=='treated' & TB=='TB',.(TTBpl=1e5*sum(mid)),by=arm]
fwrite(TTB,file=gh('graphs/TTB_{prevapproach}.csv'))



## --- run over different countries
cnmz <- names(C)
cnmz <- cnmz[cnmz!='id']
toget <- c('id','cost.soc','cost.iph','cost.idh',
           'att.soc','att.idh','att.iph',
           'deaths.soc','deaths.idh','deaths.iph',
           'LYS','value'
           )
notwt <- c('id','LYS','value') #variables not to weight against value
lyarm <- c('LYL.soc','LYL.idh','LYL.iph')
tosum <- c(setdiff(toget,notwt),lyarm)
## heuristic to scale top value for thresholds:
heur <- c('id','value','deaths.iph','deaths.soc')
out <- D[,..heur]
out <- out[,lapply(.SD,function(x) sum(x*value)),.SDcols=c('deaths.iph','deaths.soc'),by=id] #sum against popn
## topl <- 0.25/out[,mean(deaths.soc-deaths.iph)]
topl <- 300
lz <- seq(from = 0,to=topl,length.out = 1000) #threshold vector for CEACs



## containers & loop
allout <- allpout <- list() #tabular outputs
ceacl <- list()             #CEAC outputs
## cn <- isoz[1]
for(cn in isoz){
  cat('running model for:',cn,'\n')
  ## --- costs
  ## drop previous costs
  D[,c(cnmz):=NULL]
  ## add cost data
  C <- MakeCostData(allcosts[iso3==cn],nreps) #make cost PSA
  D <- merge(D,C,by='id',all.x = TRUE)        #merge into PSA
  ## --- DALYs
  ## drop any that are there
  if('LYS' %in% names(D)) D[,LYS:=NULL]
  D <- merge(D,LYKc[iso3==cn,.(age,LYS)],by='age',all.x = TRUE)        #merge into PSA
  ## --- run model (quietly)
  invisible(capture.output(D <- runallfuns(D,arm=notIPD)))
  ## --- grather outcomes
  out <- D[,..toget]
  out[,c(lyarm):=.(LYS*deaths.soc,LYS*deaths.idh,LYS*deaths.iph)] #LYL per pop by arm
  ## out[,sum(value),by=id]                                       #CHECK
  out <- out[,lapply(.SD,function(x) sum(x*value)),.SDcols=tosum,by=id] #sum against popn
  ## increments wrt SOC (per child presenting at either DH/PHC)
  out[,Dcost.iph:=cost.iph-cost.soc]; out[,Dcost.idh:=cost.idh-cost.soc] #inc costs
  out[,Datt.iph:=att.iph-att.soc]; out[,Datt.idh:=att.idh-att.soc] #inc atts
  out[,Ddeaths.iph:=deaths.iph-deaths.soc]; out[,Ddeaths.idh:=deaths.idh-deaths.soc] #inc deaths
  out[,DLYL.iph:=LYL.iph-LYL.soc]; out[,DLYL.idh:=LYL.idh-LYL.soc] #inc LYLs
  ## per whatever
  out[,DcostperATT.iph:=Dcost.iph/Datt.iph];out[,DcostperATT.idh:=Dcost.idh/Datt.idh];out[,DcostperATT.soc:=cost.soc/att.soc]
  out[,Dcostperdeaths.iph:=-Dcost.iph/Ddeaths.iph];out[,Dcostperdeaths.idh:=-Dcost.idh/Ddeaths.idh]
  out[,DcostperLYS.iph:=-Dcost.iph/DLYL.iph];out[,DcostperLYS.idh:=-Dcost.idh/DLYL.idh]
  ## summarize
  smy <- outsummary(out)
  outs <- smy$outs; pouts <- smy$pouts;
  outs[,iso3:=cn]; pouts[,iso3:=cn]
  ## capture tabular
  allout[[cn]] <- outs; allpout[[cn]] <- pouts
  ## ceac data
  ceacl[[cn]] <- data.table(iso3=cn,
                            iph=make.ceac(out[,.(Q=-DLYL.iph,P=Dcost.iph)],lz),
                            idh=make.ceac(out[,.(Q=-DLYL.idh,P=Dcost.idh)],lz),
                            threshold=lz)
}
allout <- rbindlist(allout)
allpout <- rbindlist(allpout)
ceacl <- rbindlist(ceacl)

fwrite(allout,file=gh('graphs/allout_{prevapproach}.csv'))
fwrite(allpout,file=gh('graphs/allpout_{prevapproach}.csv'))
save(ceacl,file=gh('graphs/ceacl_{prevapproach}.Rdata'))


## CEAC plot
cbPalette <- c("#999999", "#E69F00", "#56B4E9","#009E73",
               "#F0E442", "#0072B2","#D55E00", "#CC79A7")
ceaclm <- melt(ceacl,id=c('iso3','threshold'))
ceaclm[,Intervention:=ifelse(variable=='iph','PHC-focussed','DH-focussed')]
## name key
ckey <- data.table(iso3=c('KHM','CMR','CIV','MOZ','SLE','UGA','ZMB'),
                   country=c('Cambodia','Cameroon',"C\u00F4te d'Ivoire",
                             'Mozambique','Sierra Leone','Uganda','Zambia'))

ceaclm <- merge(ceaclm,ckey,by='iso3',all.x=TRUE)

## plot
GP <- ggplot(ceaclm,aes(threshold,value,
                        col=country,lty=Intervention)) +
  geom_line() +
  theme_classic() +
  theme(legend.position = 'top',legend.title = element_blank())+
  ggpubr::grids()+
  ylab('Probability cost-effective')+
  xlab('Cost-effectiveness threshold (USD/DALY)')+
  scale_colour_manual(values=cbPalette)
GP

ggsave(GP,file=gh('graphs/CEAC_{prevapproach}.png'),w=7,h=5)



## plot
GP <- ggplot(ceaclm[variable=='idh'],aes(threshold,value,
                        col=country)) +
  geom_line() +
  theme_classic() +
  theme(legend.position = 'top',legend.title = element_blank())+
  ggpubr::grids()+
  ylab('Probability cost-effective')+
  xlab('Cost-effectiveness threshold (USD/DALY)')+
  scale_colour_manual(values=cbPalette) ## + xlim(x=c(0,1500))
GP

ggsave(GP,file=gh('graphs/CEAC1_{prevapproach}.png'),w=7,h=5)


## ------ no ZMB versions of graphs ------

## plot
GP <- ggplot(ceaclm[iso3 !='ZMB'],
             aes(threshold,value,
                 col=country,lty=Intervention)) +
  geom_line() +
  theme_classic() +
  theme(legend.position = 'top',legend.title = element_blank())+
  ggpubr::grids()+
  ylab('Probability cost-effective')+
  xlab('Cost-effectiveness threshold (USD/DALY)')+
  scale_colour_manual(values=cbPalette)
GP

ggsave(GP,file=gh('graphs/CEAC_noZMB_{prevapproach}.png'),w=7,h=5)



## plot
GP <- ggplot(ceaclm[variable=='idh' & iso3 !='ZMB'],
             aes(threshold,value,
                 col=country)) +
  geom_line() +
  theme_classic() +
  theme(legend.position = 'top',legend.title = element_blank())+
  ggpubr::grids()+
  ylab('Probability cost-effective')+
  xlab('Cost-effectiveness threshold (USD/DALY)')+
  scale_colour_manual(values=cbPalette) ## + xlim(x=c(0,1500))
GP

ggsave(GP,file=gh('graphs/CEAC1_noZMB_{prevapproach}.png'),w=7,h=5)





## ## generate some CEA outputs in graphs/ & outdata/
## ## NOTE these folders need to be created
## ## NOTE need ggpubr, BCEA installed
## MakeCEAoutputs(D, #PSA dataset
##                LYK, #discounted expected life-years by age
##                file.id='test', #string to identify output files 
##                Kmax=5e3,wtp=5e3)
