## top level handler for DECENTRALIZATION
rm(list=ls())
library(here)

## load other scripts
source(here('R/decent_tree.R'))           #tree structure and namings: also tree functions & libraries
source(here('R/decent_functions.R'))      #functions for tree parameters

## number of reps
nreps <- 1e4

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
P <- c(P,P2)


## ## version making test plots
## P <- parse.parmtable(PD0,
##                      testdir = here('graphs/test'),
##                      outfile = here('graphs/test/00out.csv'))

## make base PSA dataset
set.seed(1234) #random number seed
D <- makePSA(nreps,P,dbls = list(c('cfrhivor','cfrartor')))

## use these parameters to construct intput data by attribute
D <- makeAttributes(D)
D[,sum(value),by=id] #CHECK

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

## --- cascade variables checking
nmz <- names(D)
rnmz <- grep("DH\\.|PHC\\.",nmz,value=TRUE)
rnmz <- c('id','age','tb',rnmz)
A <- D[,..rnmz]
A[,TB:=ifelse(tb!='noTB','TB','not TB')] #simpler version of TB indicator
A[,tb:=NULL]
A <- melt(A,id=c('id','age','TB'))
A[,c('location','stage','arm'):=tstrsplit(variable,split='\\.')]

## summary version
AS <- A[,.(mid=mean(value),lo=lo(value),hi=hi(value)),by=.(age,TB,arm,location,stage)]
lvls <- c('presented','screened','presumed','treated')
AS$stage <- factor(AS$stage,levels=lvls,ordered = TRUE)
tpl <- AS[stage=='presented',.(top=sum(mid)),by=.(arm,location,age)]
AS <- merge(AS,tpl[,.(arm,location,age,top)],by=c('arm','location','age'),all.x=TRUE)
AS[,vpl:=1e5*mid/top] #value per lakh
AS[,txt:=round(vpl)]
AS[stage=='presented',txt:=NA]


## graphical output
GP <- ggplot(AS,aes(stage,vpl,fill=TB))+
  geom_bar(stat='identity')+
  facet_grid(location+age ~ arm)+
  scale_y_continuous(label=comma)+
  ylab('Number per 100K attending')+
  xlab('Stage')+
  theme_light() + rot45
## GP
ggsave(GP,file=here('graphs/cascade_plt.png'),w=15,h=15)



## TODO
## errors!
## check prev
## check presume <> ATT


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
