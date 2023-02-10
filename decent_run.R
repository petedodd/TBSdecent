## top level handler for DECENTRALIZATION
## rm(list=ls())

## argument handling:
## see NOTE comments for how to set by hand for interactive use of script below
args <- commandArgs(trailingOnly = TRUE)

nargs <- length(args)
variant <- as.character(args[1])
if(variant=='BIA'){
  cat('Running BIA! \n')
  bia <- 'BIA' #NOTE set variable by hand if BIA
} else if(variant=='CEA') {
  cat('Running CEA! \n')
  bia <- '' #NOTE set variable by hand if CEA
} else {
  stop('First argument must either by BIA or CEA!')
}

if(nargs==1){
  postpend <- 'DECENT' #NOTE set by hand for base full-cost analysis
  cat('...using full costs \n')
} else{
  postpend <- as.character(args[2]) #NOTE set variable to 'DECENT_EQUIPMENT' etc if partial cost analysis
  cat('...using costs for ',postpend,'\n')
}

## some by-hand setting
## bia <- ''
## postpend <- 'DECENT'
## postpend <- 'DECENT_PERSONNEL'
## postpend <- 'DECENT_SUPPLY'

## arguments defined: begin script
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
  LYKc0 <- GetLifeYears(isolist=isoz,discount.rate=0.00,yearfrom=2021)
  LYKc <- merge(LYKc,LYKc0[,.(iso3,age,LYS0=LYS)],by=c('iso3','age'))
  LYK <- LYKc[,.(LYS=mean(LYS),LYS0=mean(LYS0)),by=.(age)] #averaged life-years 4 generic tests
  save(LYKc,file=here('indata/LYKc.Rdata'))
  save(LYK,file=here('indata/LYK.Rdata'))
} else {
  load(file=here('indata/LYKc.Rdata'))
  load(file=here('indata/LYK.Rdata'))
}


## read and make cost data
csts <- fread(here('indata/testcosts.csv'))         #read cost data
fn <- gh('indata/TBS.{postpend}.costs{bia}.csv')
rcsts <- fread(fn,skip = 1)    #read cost data
cat('!!! Using costs from ',fn,' !!!\n')

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


## calculate TB-independent cascade parameters
tbicp <- TBinptCascadeParms(CDD,ICS) #calculate from data
save(tbicp,file=here('graphs/cascades/data/tbicp.Rdata')) #NOTE for reporting

## make into parmtable
PDx <- rbind(data.frame(NAME=tbicp$prps$nm,
                         DISTRIBUTION=paste0('B(',tbicp$prps$a,',',tbicp$prps$b,')')),
              data.frame(NAME=tbicp$Fu5$nm,
                         DISTRIBUTION=paste0('B(',tbicp$Fu5$a,',',tbicp$Fu5$b,')'))
              )
PDx1 <- data.frame(NAME=tbicp$pstv$nm,DISTRIBUTION=paste0('G(',tbicp$pstv$gmk,',',tbicp$pstv$gmsc,')'))
PDx <- rbind(PDx,PDx1)

## ps - look to compute dxa for each rep of PSA
## build PSA
PD1[,2] <- "0.5" #not relevant but needed to generate answers (tbc)
P1 <- parse.parmtable(PD0)             #convert into parameter object
P2 <- parse.parmtable(PD1)             #convert into parameter object
Px <- parse.parmtable(PDx)             #the TB-indepdent cascade parms
notIPD <- c('SOC','IDH','IPH')
P2$d.TBprev.ICS.o5 <- P2$d.TBprev.ICS.u5 <- 0.9999 ## NOTE for computing sense
P <- c(P1,P2,Px)
D <- makePSA(nreps,P,dbls = list(c('cfrhivor','cfrartor')))

## rename other parms based on Px
D[,d.F.u5:=Fu5] #age mix
D[,c("d.idh.pphc.u5","d.iph.pphc.u5","d.ipd.pphc.u5",
  "d.soc.pphc.u5","d.idh.pphc.o5","d.iph.pphc.o5",
  "d.ipd.pphc.o5","d.soc.pphc.o5"):=F_ICS] #ICS
D[,c("d.soc.dh.assess.u5","d.soc.dh.assess.o5",
     "d.soc.phc.assess.u5","d.soc.phc.assess.o5"):=F_SOC_BoA_modelled] #SOC screen
D[,c("d.idh.dh.assess.u5","d.idh.dh.assess.o5",
  "d.iph.dh.assess.u5","d.iph.dh.assess.o5",
  "d.idh.phc.assess.u5","d.idh.phc.assess.o5",
  "d.iph.phc.assess.u5","d.iph.phc.assess.o5"):=F_alpha] #INT screen
## & the following are *over* written
D[,c("d.idh.rltfu","d.iph.rltfu","d.ipd.rltfu","d.soc.rltfu"):=F_refu] #TODO check
D[,c("d.iph.phc.test7.referDH","d.ipd.phc.test.referDH",
     "d.ipd.phc.notest.referDH","d.soc.phc.test.referDH",
     "d.soc.phc.notest.referDH"):=F_refs] #TODO check
## D[,f:=F_omega_flat] NOTE only used in TB cascade calx


## spec and prev priors
prior.tbprev <- list(mean=odds(0.1),sd=0.03) #normal in odds
prior.phi <- list(mean=(0.1),sd=0.025) #normal: phi = 1-spec clinical
D[,tbprev:=iodds(rnorm(nreps,mean=prior.tbprev$mean,sd=prior.tbprev$sd))]
D[,phi:=rnorm(nreps,mean=prior.phi$mean,sd=prior.phi$sd)]
D$spec.clin <- D$spec.clinCXR.soc <- pmin(1,1-D$phi)

## DI <- copy(D)
## AM[arm=='idh' & location=='PHC']
## ## d.soc.phc.ptbxsp
## D[,.(d.idh.pphc,d.idh.phc.assess,d.phc.prsmptv,d.idh.phc.presumed)]


## add in relevant parameters:
DxAZ <- computeDxAccuracy4PSA(D)
addon <- D[,.(tbprev,phi,
              BoA=F_alpha,CoB=F_presume,
              `DoC_5-14`,`DoC_0-4`,
              F_ICS)] #corresonding parameters
addon[,id:=1:nreps]
DxAZ <- merge(DxAZ,addon,by=c('id'),all.x=TRUE)
## DxAZ[(arm=='idh' & location=='PHC')] #NOTE check
DxAZ <- DxAZ[!(arm=='idh' & location=='PHC')] #NOTE check

## linear approx for algorithm specificity as fn of phi
## 1-spec = a + b*phi
## 1-spec1 = a + b*0
## 1-spec = a + b*phi
## a = 1-spec1; b = (1-spec-a)/phi
DxAZ[,a:=1-spec1]
DxAZ[,b:=(1-spec-a)/phi]

## age-specific DoC
DxAZ[,DoC:=ifelse(age=='0-4',`DoC_0-4`,`DoC_5-14`)]


## ## NOTE test
## DxAZ[,DoC:=0.1] #test value for final cascade step now

## ## sample conditional on D/C cascade value
## BayesSpecPrev(DxAZ) #acts by side-effect

## ## checks
## DxAZ[,summary(tbprev)]
## DxAZ[,summary(phi)]
## DxAZ[,summary(sense*tbprev+(a+b*phi)*(1-tbprev))] #at level

## ## TODO
## ## BoA is alpha
## ## CoB is from F_presume
## ## DoC is TPS.csv

## ## TODO real DoC and next steps
## DDW <- dcast(DD,arm+location+age~stage,value.var='vpl')
## DxAZ <- merge(DxAZ,DDW[,.(arm,location,age,
##                           DoC=treated/presumed,CoB=presumed/screened,BoA=screened/1e5)],
##               by=c('arm','location','age'),all.x=TRUE)

BayesSpecPrev(DxAZ) #acts by side-effect

## ## checks
## DxAZ[,summary(tbprev)]
## DxAZ[,summary(phi)]
## DxAZ[,summary(sense*tbprev+(a+b*phi)*(1-tbprev)-DoC)] #at level

## ## looking at results TODO save some out
## DxAZ[,.(clinspec.m=mean(1-phi),clinspec.sd=sd(1-phi),
##         prev.m=mean(tbprev),prev.s=mean(tbprev)),by=.(age)]

## DxAZ[,.(clinspec.m=mean(1-phi),clinspec.sd=sd(1-phi),
##         prev.m=mean(tbprev),prev.s=mean(tbprev)),by=.(arm)]

## DxAZ[,.(clinspec.m=mean(1-phi),clinspec.sd=sd(1-phi),
##         prev.m=mean(tbprev),prev.s=mean(tbprev)),by=.(location)]


## tmp <- DxAZ[,.(clinspec.m=mean(1-phi),clinspec.sd=sd(1-phi),
##                prev.m=mean(tbprev),prev.s=mean(tbprev)),by=.(arm,location,age)]

## tmp[order(age,location,arm)][arm!='idh']
## tmp[order(age,arm,location)][arm!='idh']

## ## version with SOC and IPH sharing prevalence
## Dx2arm <- dcast(DxAZ[arm!='idh'],age+location+id ~ arm,value.var = c('sense','a','b','DoC'))
## setnames(Dx2arm,old=c("sense_iph","sense_soc",
##                       "a_iph","a_soc",
##                       "b_iph","b_soc",
##                       "DoC_iph","DoC_soc"),
##          new=c('sense2','sense1','a2','a1','b2','b1','DoC2','DoC1'))

## BayesSpecPrev2(Dx2arm)

## Dx2arm


DxAZ





## TODO choose value for X-arm comparisons

## this computes and saves out the average accuracy of dx cascades
if(!file.exists(here('graphs'))) dir.create(here('graphs'))
## DxA <- computeDxAccuracy(PD0,PD1,C,nreps)

## ## this computes and saves model parameters derived from cascade data
## prevapproach <- 'gm'
## ## PD1 <- computeCascadeParameters(DD,ICS,DxA,prevapproach)
## PD1 <- read.csv(here('indata/calcparms_gm_in.csv'))


## ## combine different parameter types
## P1 <- parse.parmtable(PD0)             #convert into parameter object
## P2 <- parse.parmtable(PD1)             #convert into parameter object
## P <- c(P1,P2)
## names(P)

## ## make base PSA dataset

## D <- makePSA(nreps,P,dbls = list(c('cfrhivor','cfrartor')))

## ## NOTE temporary introduction of noise:
## for(nm in setdiff(PD1$NAME,'d.OR.dh.if.TB')){ #loop over probs
##   D[[nm]] <- ilogit(logit(D[[nm]]) + rnorm(nreps)/5)
## }
## D[['d.OR.dh.if.TB']] <- exp(log(D[['d.OR.dh.if.TB']]) + rnorm(nreps)/5)

Dx2arm
## jj
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

ggsave(GP2,file=gh('graphs/{bia}cascade_plt_{prevapproach}.{postpend}.png'),w=15,h=15)

## Treated true TB per 100K presented by arm
TTB <- AS[stage=='treated' & TB=='TB',.(TTBpl=1e5*sum(mid)),by=arm]
fwrite(TTB,file=gh('graphs/{bia}TTB_{prevapproach}.{postpend}.csv'))



## --- run over different countries
cnmz <- names(C)
cnmz <- cnmz[cnmz!='id']
costsbystg <- c('DH.presented.cost.soc','DH.presented.cost.iph','DH.presented.cost.idh',
                'DH.screened.cost.soc','DH.screened.cost.iph','DH.screened.cost.idh',
                'DH.presumed.cost.soc','DH.presumed.cost.iph','DH.presumed.cost.idh',
                'DH.treated.cost.soc','DH.treated.cost.iph','DH.treated.cost.idh',
                'PHC.presented.cost.soc','PHC.presented.cost.iph','PHC.presented.cost.idh',
                'PHC.screened.cost.soc','PHC.screened.cost.iph','PHC.screened.cost.idh',
                'PHC.presumed.cost.soc','PHC.presumed.cost.iph','PHC.presumed.cost.idh',
                'PHC.treated.cost.soc','PHC.treated.cost.iph','PHC.treated.cost.idh')
toget <- c('id',
           'cost.soc','cost.iph','cost.idh',
           costsbystg,
           'att.soc','att.idh','att.iph',
           'deaths.soc','deaths.idh','deaths.iph',
           'LYS','LYS0','value'
           )
notwt <- c('id','LYS','LYS0','value') #variables not to weight against value
lyarm <- c('LYL.soc','LYL.idh','LYL.iph')
lyarm <- c(lyarm,gsub('\\.','0\\.',lyarm)) #include undiscounted
tosum <- c(setdiff(toget,notwt),lyarm)
## heuristic to scale top value for thresholds:
heur <- c('id','value','deaths.iph','deaths.soc')
out <- D[,..heur]
out <- out[,lapply(.SD,function(x) sum(x*value)),.SDcols=c('deaths.iph','deaths.soc'),by=id] #sum against popn
## topl <- 0.25/out[,mean(deaths.soc-deaths.iph)]
topl <- 300 #100
lz <- seq(from = 0,to=topl,length.out = 1000) #threshold vector for CEACs
## staged costs by arm
soc.sc <- grep('soc',costsbystg,value=TRUE); psoc.sc <- paste0('perATT.',soc.sc)
iph.sc <- grep('iph',costsbystg,value=TRUE); piph.sc <- paste0('perATT.',iph.sc)
idh.sc <- grep('idh',costsbystg,value=TRUE); pidh.sc <- paste0('perATT.',idh.sc)


## containers & loop
allout <- allpout <- allscout <- list() #tabular outputs
ceacl <- NMB <- list()             #CEAC outputs etc
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
  if('LYS' %in% names(D)) D[,c('LYS','LYS0'):=NULL]
  D <- merge(D,LYKc[iso3==cn,.(age,LYS,LYS0)],by='age',all.x = TRUE)        #merge into PSA
  ## --- run model (quietly)
  invisible(capture.output(D <- runallfuns(D,arm=notIPD)))
  ## --- grather outcomes
  out <- D[,..toget]
  out[,c(lyarm):=.(LYS*deaths.soc,LYS*deaths.idh,LYS*deaths.iph,
                   LYS0*deaths.soc,LYS0*deaths.idh,LYS0*deaths.iph)] #LYL per pop by arm
  ## out[,sum(value),by=id]                                       #CHECK
  out <- out[,lapply(.SD,function(x) sum(x*value)),.SDcols=tosum,by=id] #sum against popn
  ## non-incremental cost per ATT
  out[,costperATT.soc:=cost.soc/att.soc];
  out[,(psoc.sc):=lapply(.SD,function(x) x/att.soc),.SDcols=soc.sc]
  out[,costperATT.iph:=cost.iph/att.iph];
  out[,(piph.sc):=lapply(.SD,function(x) x/att.iph),.SDcols=iph.sc]
  out[,costperATT.idh:=cost.idh/att.idh];
  out[,(pidh.sc):=lapply(.SD,function(x) x/att.idh),.SDcols=idh.sc]
  ## increments wrt SOC (per child presenting at either DH/PHC)
  out[,Dcost.iph:=cost.iph-cost.soc]; out[,Dcost.idh:=cost.idh-cost.soc] #inc costs
  out[,Datt.iph:=att.iph-att.soc]; out[,Datt.idh:=att.idh-att.soc] #inc atts
  out[,attPC.iph:=1e2*att.iph/att.soc]; out[,attPC.idh:=1e2*att.idh/att.soc] #rel inc atts
  out[,Ddeaths.iph:=deaths.iph-deaths.soc]; out[,Ddeaths.idh:=deaths.idh-deaths.soc] #inc deaths
  out[,DLYL0.iph:=LYL0.iph-LYL0.soc]; out[,DLYL0.idh:=LYL0.idh-LYL0.soc] #inc LYLs w/o discount
  out[,DLYL.iph:=LYL.iph-LYL.soc]; out[,DLYL.idh:=LYL.idh-LYL.soc] #inc LYLs
  ## per whatever
  out[,DcostperATT.iph:=cost.iph/att.iph-cost.soc/att.soc];
  out[,DcostperATT.idh:=cost.idh/att.idh-cost.soc/att.soc];
  out[,Dcostperdeaths.iph:=-cost.iph/deaths.iph+cost.soc/deaths.soc];
  out[,Dcostperdeaths.idh:=-cost.idh/deaths.idh+cost.soc/deaths.soc]
  out[,DcostperLYS0.iph:=-cost.iph/LYL0.iph+cost.soc/LYL0.soc];
  out[,DcostperLYS0.idh:=-cost.idh/LYL0.idh+cost.soc/LYL0.soc]
  out[,DcostperLYS.iph:=-cost.iph/LYL.iph+cost.soc/LYL.soc];
  out[,DcostperLYS.idh:=-cost.idh/LYL.idh+cost.soc/LYL.soc]
  ## D/D
  out[,DcostperDATT.iph:=Dcost.iph/Datt.iph];out[,DcostperDATT.idh:=Dcost.idh/Datt.idh];
  out[,DcostperDdeaths.iph:=-Dcost.iph/Ddeaths.iph];out[,DcostperDdeaths.idh:=-Dcost.idh/Ddeaths.idh]
  out[,DcostperDLYS0.iph:=-Dcost.iph/DLYL0.iph];out[,DcostperDLYS0.idh:=-Dcost.idh/DLYL0.idh]
  out[,DcostperDLYS.iph:=-Dcost.iph/DLYL.iph];out[,DcostperDLYS.idh:=-Dcost.idh/DLYL.idh]
  ## summarize
  smy <- outsummary(out)
  outs <- smy$outs; pouts <- smy$pouts; scouts <- smy$scouts
  outs[,iso3:=cn]; pouts[,iso3:=cn]; scouts[,iso3:=cn]
  ## capture tabular
  allout[[cn]] <- outs; allpout[[cn]] <- pouts; allscout[[cn]] <- scouts
  ## capture data for NMB
  NMB[[cn]] <- out[,.(iso3=cn,DLYL.iph,Dcost.iph,DLYL.idh,Dcost.idh)]
  ## ceac data
  ceacl[[cn]] <- data.table(iso3=cn,
                            iph=make.ceac(out[,.(Q=-DLYL.iph,P=Dcost.iph)],lz),
                            idh=make.ceac(out[,.(Q=-DLYL.idh,P=Dcost.idh)],lz),
                            threshold=lz)
}
allout <- rbindlist(allout)
allpout <- rbindlist(allpout)
allscout <- rbindlist(allscout)
ceacl <- rbindlist(ceacl)
NMB <- rbindlist(NMB)

## checks
out[,.(att.iph/att.soc,att.idh/att.soc)]
out[,.(Datt.iph/att.soc,Datt.idh/att.soc)]
out[,.(att.soc,att.iph,att.idh)]
out[,.(Ddeaths.iph/Datt.iph,Ddeaths.idh/Datt.idh)] #OK
out[,.(DLYL.iph/Ddeaths.iph,DLYL.idh/Ddeaths.idh)] #OK
## check
allout[,.(costperATT.iph.mid-costperATT.soc.mid,DcostperATT.iph.mid)]



fwrite(allout,file=gh('graphs/{bia}allout_{prevapproach}.{postpend}.csv'))
fwrite(allpout,file=gh('graphs/{bia}allpout_{prevapproach}.{postpend}.csv'))
save(ceacl,file=gh('graphs/{bia}ceacl_{prevapproach}.{postpend}.Rdata'))
save(NMB,file=gh('graphs/{bia}NMB_{prevapproach}.{postpend}.Rdata'))
save(allscout,file=gh('graphs/{bia}allscout_{prevapproach}.{postpend}.Rdata'))




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



## plot: IPH only
GP <- ggplot(ceaclm[variable=='iph' &
                    iso3 %in% c('KHM', 'CIV', 'CMR', 'MOZ', 'SLE', 'UGA')],
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

ggsave(GP,file=gh('graphs/{bia}CEAC_IPHonly_{prevapproach}.{postpend}.png'),w=7,h=5)


## plot: IDH only
GP <- ggplot(ceaclm[variable=='idh'  &
                    iso3 %in% c('KHM', 'CIV', 'CMR')],aes(threshold,value,
                        col=country,lty=Intervention)) +
  geom_line() +
  theme_classic() +
  theme(legend.position = 'top',legend.title = element_blank())+
  ggpubr::grids()+
  ylab('Probability cost-effective')+
  xlab('Cost-effectiveness threshold (USD/DALY)')+
  scale_colour_manual(values=cbPalette)
GP

ggsave(GP,file=gh('graphs/{bia}CEAC_IDHonly_{prevapproach}.{postpend}.png'),w=7,h=5)


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

ggsave(GP,file=gh('graphs/{bia}CEAC_{prevapproach}.{postpend}.png'),w=7,h=5)


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

ggsave(GP,file=gh('graphs/{bia}CEAC1_{prevapproach}.{postpend}.png'),w=7,h=5)


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

ggsave(GP,file=gh('graphs/{bia}CEAC_noZMB_{prevapproach}.{postpend}.png'),w=7,h=5)



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

ggsave(GP,file=gh('graphs/{bia}CEAC1_noZMB_{prevapproach}.{postpend}.png'),w=7,h=5)





## ## generate some CEA outputs in graphs/ & outdata/
## ## NOTE these folders need to be created
## ## NOTE need ggpubr, BCEA installed
## MakeCEAoutputs(D, #PSA dataset
##                LYK, #discounted expected life-years by age
##                file.id='test', #string to identify output files 
##                Kmax=5e3,wtp=5e3)
