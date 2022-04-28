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

## read and make cost data
csts <- fread(here('indata/testcosts.csv'))         #read cost data
C <- MakeCostData(csts,nreps)               #make cost PSA

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

## --- life years and other outputs
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
