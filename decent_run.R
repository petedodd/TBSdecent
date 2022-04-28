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

## read and make cost data
csts <- fread(here('indata/testcosts.csv'))         #read cost data
rcsts <- fread(here('indata/TBS.DECENT.costs.csv'),skip = 1)    #read cost data
## check
setdiff(unique(rcsts$NAME),
        unique(csts$cost))
setdiff(unique(csts$cost),
        unique(rcsts$NAME))
allcosts <- reformatCosts(rcsts)
C <- MakeCostData(allcosts[iso3=='CIV'],nreps)               #make cost PSA

## prior parameters
PD0 <- read.csv(here('indata/DecentParms - distributions.csv')) #read in
## parameters to be determined from cascade data
PD1 <- PD0[PD0$DISTRIBUTION=="",]
## the rest
PD0 <- PD0[PD0$DISTRIBUTION!="",]

## this computes and saves out the average accuracy of dx cascades
DxA <- computeDxAccuracy(PD0,PD1,C,nreps)

## this computes and saves model parameters derived from cascade data
PD1 <- computeCascadeParameters(DD,ICS,DxA)

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

ggsave(GP2,file=here('graphs/cascade_plt.png'),w=15,h=15)

## Treated true TB per 100K presented by arm
TTB <- AS[stage=='treated' & TB=='TB',.(TTBpl=1e5*sum(mid)),by=arm]
fwrite(TTB,file=here('graphs/TTB.csv'))

## --- run over different countries
cnmz <- names(C)
cnmz <- cnmz[cnmz!='id']
toget <- c('id','cost.soc','cost.iph','cost.idh','att.soc','att.idh','att.iph')

allout <- list()
## cn <- isoz[1]
for(cn in isoz){
  cat('running model for:',cn,'\n')
  ## drop previous costs
  D[,c(cnmz):=NULL]
  ## add cost data
  C <- MakeCostData(allcosts[iso3==cn],nreps) #make cost PSA
  D <- merge(D,C,by='id',all.x = TRUE)        #merge into PSA
  ## run model (quietly)
  invisible(capture.output(D <- runallfuns(D,arm=notIPD)))
  ## grather outcomes
  out <- D[,..toget]
  out <- out[,lapply(.SD,sum),.SDcols=toget[-1],by=id]
  out <- out[,.(DcostperATT.iph=(cost.iph-cost.soc)/(att.iph-att.soc),
                DcostperATT.idh=(cost.idh-cost.soc)/(att.idh-att.soc)),
             by=id]
  out <- out[,.(DcostperATT.iph.mid=mean(DcostperATT.iph),
                DcostperATT.iph.lo=lo(DcostperATT.iph),
                DcostperATT.iph.hi=hi(DcostperATT.iph),
                DcostperATT.idh.mid=mean(DcostperATT.idh),
                DcostperATT.idh.lo=lo(DcostperATT.idh),
                DcostperATT.idh.hi=hi(DcostperATT.idh))]
  out[,iso3:=cn]
  ## capture
  allout[[cn]] <- out
}
allout <- rbindlist(allout)
allout[,pty.iph:=paste0(round(DcostperATT.iph.mid,0),' (',
                        round(DcostperATT.iph.lo,0),' - ',
                        round(DcostperATT.iph.hi,0),')')]

allout[,pty.idh:=paste0(round(DcostperATT.idh.mid,0),' (',
                        round(DcostperATT.idh.lo,0),' - ',
                        round(DcostperATT.idh.hi,0),')')]

fwrite(allout,file=here('graphs/allout.csv'))

## TODO check total population
## TODO DALY outputs etc


## ## --- life years and other outputs
## LYSdone <- TRUE
## if(!LYSdone){
##   ## make discounted life-years if they haven't been done
##   
##   LYK <- GetLifeYears(isolist=isoz,discount.rate=0.03,yearfrom=2021)
##   save(LYK,file=here('indata/LYK.Rdata'))
## } else {load(file=here('indata/LYK.Rdata'))}
## LYK <- LYK[,.(LYS=mean(LYS)),by=age] #averaged life-years 4 generic tests

## ## generate some CEA outputs in graphs/ & outdata/
## ## NOTE these folders need to be created
## ## NOTE need ggpubr, BCEA installed
## MakeCEAoutputs(D, #PSA dataset
##                LYK, #discounted expected life-years by age
##                file.id='test', #string to identify output files 
##                Kmax=5e3,wtp=5e3)
