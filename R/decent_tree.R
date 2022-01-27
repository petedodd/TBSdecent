## looking at tree 7
## rm(list=ls())
library(here)
library(HEdtree)
library(data.tree)
library(data.table)
library(glue)


## === outcomes subtree ===
notbdxo <- txt2tree(here('indata/tbnotx.txt')) #no tx
tbtxb <- txt2tree(here('indata/tbdxb.txt')) #bac+
tbtxc <- txt2tree(here('indata/tbdxc.txt')) #clin

## default prob/cost:
notbdxo$Set(p=1)
tbtxb$Set(p=1)
tbtxc$Set(p=1)
notbdxo$Set(cost=0)
tbtxb$Set(cost=0)
tbtxc$Set(cost=0)

## set probabilities
## NOTE these namings actually get overwritten by CSV read-ins below
## -- clinical:
tbtxc$`No TB treatment`$p <- 'p.ptltfu'
tbtxc$`No TB treatment`$Dies$p <- 'p.cfr.notx'
tbtxc$`No TB treatment`$Survives$p <- '1-p.cfr.notx'
tbtxc$`RifS-TB treatment`$p <- '1-p.ptltfu'
tbtxc$`RifS-TB treatment`$Dies$p <- 'p.cfr.tx'
tbtxc$`RifS-TB treatment`$Survives$p <- '1-p.cfr.tx'

## -- bac:
## tx
tbtxb$`TB treatment`$p <- '1-p.ptltfu'
## RS
tbtxb$`TB treatment`$`RifS-TB treatment`$p <- '1-p.rr'
tbtxb$`TB treatment`$`RifS-TB treatment`$Survives$p <- '1-p.cfr.tx'
tbtxb$`TB treatment`$`RifS-TB treatment`$Dies$p <- 'p.cfr.tx'
## RR
tbtxb$`TB treatment`$`RifR-TB treatment`$p <- 'p.rr'
tbtxb$`TB treatment`$`RifR-TB treatment`$Survives$p <- '1-p.cfr.tx'
tbtxb$`TB treatment`$`RifR-TB treatment`$Dies$p <- 'p.cfr.tx'
## untreated
tbtxb$`No TB treatment`$p <- 'p.ptltfu'
tbtxb$`No TB treatment`$Survives$p <- '1-p.cfr.notx'
tbtxb$`No TB treatment`$Dies$p <- 'p.cfr.notx'

## -- no tx:
notbdxo$`No TB treatment`$Dies$p <- 'p.cfr.notx'
notbdxo$`No TB treatment`$Survives$p <- '1-p.cfr.notx'

## restrict to no TB tx (rather than dx)
notbtxo <- top(notbdxo) #remove top


## ====== function to add outcomes & counters
AddOutcomes <- function(D){
  ## === cost and probs (defaults)
  D$Set(p=1)
  D$Set(cost=0)

  ## === merge to create final tree ===
  MergeByName(D,notbtxo,'No TB diagnosed',leavesonly = TRUE) #NOTE need to restrict to leaves, although not necessary
  MergeByName(D,notbtxo,'All other patients',leavesonly = TRUE) #NOTE need to restrict to leaves, although not necessary
  MergeByName(D,notbtxo,'Unidentified presumptive TB',leavesonly = TRUE) #NOTE need to restrict to leaves, although not necessary
  MergeByName(D,notbtxo,'Does not reach DH',leavesonly = TRUE) #NOTE need to restrict to leaves, although not necessary
  MergeByName(D,tbtxb,'TB diagnosed (bacteriological)')
  MergeByName(D,tbtxc,'TB diagnosed (clinical)')

  ## ===========  other counters
  ## check
  D$Set(check=1)
  D$Set(check=0,filterFun=function(x) length(x$children)>0)

  ## deaths
  D$Set(deaths=0)
  D$Set(deaths=1,filterFun=function(x) (x$name=='Dies'))

  ## lives
  D$Set(lives=0)
  D$Set(lives=1,filterFun=function(x) (x$name=='Survives'))

  ## dx clinical
  D$Set(dxc=0)
  D$Set(dxc=1,filterFun=function(x) x$name=='TB diagnosed (clinical)')

  ## dx bac
  D$Set(dxb=0)
  D$Set(dxb=1,
        filterFun=function(x)x$name=='TB diagnosed (bacteriological)')
  ## ATT
  D$Set(att=0)
  D$Set(att=1,
        filterFun=function(x)x$name %in% c('RifS-TB treatment','RifR-TB treatment'))

  ## referrals
  D$Set(refers=0)
  D$Set(refers=1,filterFun=function(x) grepl('Refer',x$name))

  return(D)
}


print(SOC,'p','cost','deaths','lives','refers','dxc','dxb','att','check')

## === SOC
SOC <- MSorg2tree(here('indata/dSOC.txt'))
SOC <- top(SOC)
print(SOC)
## merge in extras, write out
SOC <- AddOutcomes(SOC)

tree2file(SOC,filename = here('indata/CSV/SOC0.csv'),
          'p','cost','deaths','lives','refers','dxc','dxb','att','check')

## create version with probs/costs
fn <- here('indata/CSV/SOC1.csv')
if(file.exists(fn)){
  ## read
  labz <- fread(fn)
  SOC$Set(p=labz$p)
  SOC$Set(cost=labz$cost)
  tree2file(SOC,filename = here('indata/CSV/SOC2.csv'),
            'p','cost','deaths','lives','refers','dxc','dxb','att','check')
}

## === IPD
IPD <- MSorg2tree(here('indata/dIPD.txt'))
IPD <- top(IPD)
print(IPD)
IPD <- AddOutcomes(IPD)
tree2file(IPD,filename = here('indata/CSV/IPD0.csv'),
          'p','cost','deaths','lives','refers','dxc','dxb','att','check')

## create version with probs/costs
fn <- here('indata/CSV/IPD1.csv')
if(file.exists(fn)){
    ## read
    labz <- fread(fn)
    IPD$Set(p=labz$p)
    IPD$Set(cost=labz$cost)
  tree2file(IPD,filename = here('indata/CSV/IPD2.csv'),
            'p','cost','deaths','lives','refers','dxc','dxb','att','check')
}

## === IDH
IDH <- MSorg2tree(here('indata/dIDH.txt'))
IDH <- top(IDH)
print(IDH)
IDH <- AddOutcomes(IDH)
tree2file(IDH,filename = here('indata/CSV/IDH0.csv'),
          'p','cost','deaths','lives','refers','dxc','dxb','att','check')

## create version with probs/costs
fn <- here('indata/CSV/IDH1.csv')
if(file.exists(fn)){
    ## read
    labz <- fread(fn)
    IDH$Set(p=labz$p)
    IDH$Set(cost=labz$cost)
  tree2file(IDH,filename = here('indata/CSV/IDH2.csv'),
            'p','cost','deaths','lives','refers','dxc','dxb','att','check')
}

## === IPH
IPH <- MSorg2tree(here('indata/dIPH.txt'))
IPH <- top(IPH)
print(IPH)
IPH <- AddOutcomes(IPH)
tree2file(IPH,filename = here('indata/CSV/IPH0.csv'),
          'p','cost','deaths','lives','refers','dxc','dxb','att','check')

## create version with probs/costs
fn <- here('indata/CSV/IPH1.csv')
if(file.exists(fn)){
    ## read
    labz <- fread(fn)
    IPH$Set(p=labz$p)
    IPH$Set(cost=labz$cost)
  tree2file(IPH,filename = here('indata/CSV/IPH2.csv'),
            'p','cost','deaths','lives','refers','dxc','dxb','att','check')
}


## make functions
SOC.F <- makeTfuns(SOC,c('check','cost','deaths','att',
                         'lives','refers','dxc','dxb'))
IPD.F <- makeTfuns(IPD,c('check','cost','deaths','att',
                         'lives','refers','dxc','dxb'))
IDH.F <- makeTfuns(IDH,c('check','cost','deaths','att',
                         'lives','refers','dxc','dxb'))
IPH.F <- makeTfuns(IPH,c('check','cost','deaths','att',
                         'lives','refers','dxc','dxb'))

runallfuns <- function(Data,arm='all'){
    done <- FALSE
    if(arm=='SOC'){
        print('Running functions for SOC arm!')
        Data$check.soc <- SOC.F$checkfun(Data); print('check run!')
        Data$lives.soc <- SOC.F$livesfun(Data); print('lives run!')
        Data$deaths.soc <- SOC.F$deathsfun(Data); print('deaths run!')
        Data$cost.soc <- SOC.F$costfun(Data); print('costs run!')
        Data$refers.soc <- SOC.F$refersfun(Data); print('refers run!')
        Data$dxb.soc <- SOC.F$dxbfun(Data); print('dxb run!')
        Data$dxc.soc <- SOC.F$dxcfun(Data); print('dxc run!')
        done <- TRUE
    }
    if(arm=='IPD'){
        print('Running functions for IPD arm!')
        Data$check.spd <- IPD.F$checkfun(Data); print('check run!')
        Data$lives.spd <- IPD.F$livesfun(Data); print('lives run!')
        Data$deaths.spd <- IPD.F$deathsfun(Data); print('deaths run!')
        Data$cost.spd <- IPD.F$costfun(Data); print('costs run!')
        Data$refers.spd <- IPD.F$refersfun(Data); print('refers run!')
        Data$dxb.spd <- IPD.F$dxbfun(Data); print('dxb run!')
        Data$dxc.spd <- IPD.F$dxcfun(Data); print('dxc run!')
        done <- TRUE
    }
    if(arm=='IDH'){
        print('Running functions for IDH arm!')
        Data$check.idh <- IDH.F$checkfun(Data); print('check run!')
        Data$lives.idh <- IDH.F$livesfun(Data); print('lives run!')
        Data$deaths.idh <- IDH.F$deathsfun(Data); print('deaths run!')
        Data$cost.idh <- IDH.F$costfun(Data); print('costs run!')
        Data$refers.idh <- IDH.F$refersfun(Data); print('refers run!')
        Data$dxb.idh <- IDH.F$dxbfun(Data); print('dxb run!')
        Data$dxc.idh <- IDH.F$dxcfun(Data); print('dxc run!')
        done <- TRUE
    }
    if(arm=='IPH'){
        print('Running functions for IPH arm!')
        Data$check.iph <- IPH.F$checkfun(Data); print('check run!')
        Data$lives.iph <- IPH.F$livesfun(Data); print('lives run!')
        Data$deaths.iph <- IPH.F$deathsfun(Data); print('deaths run!')
        Data$cost.iph <- IPH.F$costfun(Data); print('costs run!')
        Data$refers.iph <- IPH.F$refersfun(Data); print('refers run!')
        Data$dxb.iph <- IPH.F$dxbfun(Data); print('dxb run!')
        Data$dxc.iph <- IPH.F$dxcfun(Data); print('dxc run!')
        done <- TRUE
    }
    if(arm=='all'){
        Data <- runallfuns(Data,arm='SOC') #standard of care
        Data <- runallfuns(Data,arm='IPD') #soc, partially decent'd
        Data <- runallfuns(Data,arm='IDH') #DH-focused intervention
        Data <- runallfuns(Data,arm='IPH') #PHC-focused intervention
        done <- TRUE
    }
    if(!done)stop('Functions not all run! Likely unrecognised arm supplied.')
    return(Data)
}





## --- make functions



## checking
vrz.soc <- showParmz(SOC)$vars
vrz.spd <- showParmz(IPD)$vars
vrz.idh <- showParmz(IDH)$vars
vrz.iph <- showParmz(IPH)$vars
vrz <- c(vrz.soc,
         vrz.spd,
         vrz.idh,
         vrz.iph
         )
vrz <- unique(vrz)
vrz <- c(vrz,
         'dh.prsmptv',
         'phc.prsmptv',
         'soc.dh.presumed',
         'soc.dh.test',
         'soc.dh.fracsp',
         'soc.phc.presumed',
         'soc.phc.referDH2',
         'soc.phc.referDH',
         'ipd.dh.presumed',
         'ipd.dh.ptbc',
         'ipd.dhref.ptbc',
         'ipd.phc.presumed',
         'ipd.phc.referDH',
         'idh.dh.ptbc',
         'idh.phc.presumed',
         'idh.dh.presumed',
         'iph.phc.presumed',
         'iph.dh.presumed',
         'iph.dh.ptbc',
         'iph.phc.ptbc'
         )



## 
ncheck <- 100
A <- data.table(vrz,value=runif(length(vrz)))
A <- A[rep(1:length(vrz),each=ncheck)]
idz <- rep(1:ncheck,length(vrz))
A[,id:=idz]
A[,value:=runif(nrow(A))]
## 
A[,.N,by=vrz]
A[,min(id),by=vrz]
A[,max(id),by=vrz]
## 
A <- dcast(A,id~vrz,value.var = 'value')
dim(A)
## 
## random test variables
A <- data.table(vrz,value=runif(length(vrz)))
## A <- dcast(A,.~vrz,value.var = 'value') #TODO
A <- transpose(A,make.names = 'vrz')

## A$ipd.p.phc <- 1 #bug when @ PHC
## A$ipd.phc.presumed <- 1

## A$dh.prsmptv <- 0



## checks
IPD.F$checkfun(A) #NOTE OK

IPH.F$checkfun(A) #NOTE OK

IDH.F$checkfun(A) #NOTE OK

SOC.F$checkfun(A) #NOTE OK


## TODO port over model-runner, include model functions, simplify multirun; done SOC rest to do before removing SOC etc

## ## running all function
## runallfuns <- function(D,arm='both'){
##   done <- FALSE
##   if(arm=='INT' | arm=='both'){
##     cat('Running functions for intervention arm:\n')
##     for(nm in names(treefunsINT)){
##       snm <- gsub('fun','',nm)
##       snma <- paste0(snm,'INT')
##       D[[snma]] <- treefunsINT[[nm]](D)
##       cat(snm,' run...\n')
##       done <- TRUE
##     }
##   }
##   if(arm=='SOC' | arm=='both'){
##     cat('Running functions for SOC arm:\n')
##     for(nm in names(treefunsSOC)){
##       snm <- gsub('fun','',nm)
##       snma <- paste0(snm,'SOC')
##       D[[snma]] <- treefunsSOC[[nm]](D)
##       cat(snm,' run...\n')
##       done <- TRUE
##     }
##   }
##   if(!done)stop('Functions not run! Likely unrecognised arm supplied.')
##   return(D)
## }


