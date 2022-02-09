## looking at tree 7
rm(list=ls())
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
  labz$p <- gsub("p\\.rr","prr",labz$p) #NOTE fixing typo
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
  labz$p <- gsub("p\\.rr","prr",labz$p) #NOTE fixing typo
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
  labz$p <- gsub("p\\.rr","prr",labz$p) #NOTE fixing typo
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
  labz$p <- gsub("p\\.rr","prr",labz$p) #NOTE fixing typo
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


## running all function
runallfuns <- function(D,arm='both'){
  done <- FALSE
  if(arm=='SOC' | arm=='all'){
    cat('Running functions for SOC:\n')
    for(nm in names(SOC.F)){
      snm <- gsub('fun','',nm)
      snma <- paste0(snm,'.soc')
      D[[snma]] <- SOC.F[[nm]](D)
      cat(snm,' run...\n')
      done <- TRUE
    }
  }
  if(arm=='IPD' | arm=='all'){
    cat('Running functions for IPD:\n')
    for(nm in names(IPD.F)){
      snm <- gsub('fun','',nm)
      snma <- paste0(snm,'.ipd')
      D[[snma]] <- IPD.F[[nm]](D)
      cat(snm,' run...\n')
      done <- TRUE
    }
  }
  if(arm=='IDH' | arm=='all'){
    cat('Running functions for IDH:\n')
    for(nm in names(IDH.F)){
      snm <- gsub('fun','',nm)
      snma <- paste0(snm,'.idh')
      D[[snma]] <- IDH.F[[nm]](D)
      cat(snm,' run...\n')
      done <- TRUE
    }
  }
  if(arm=='IPH' | arm=='all'){
    cat('Running functions for IPH:\n')
    for(nm in names(IPH.F)){
      snm <- gsub('fun','',nm)
      snma <- paste0(snm,'.iph')
      D[[snma]] <- IPH.F[[nm]](D)
      cat(snm,' run...\n')
      done <- TRUE
    }
  }
  if(!done)stop('Functions not run! Likely unrecognised arm supplied.')
  return(D)
}



## --- CHECKS
showAllParmz <- function(TREE){
  B <- showParmz(TREE)
  ## get calx
  cx <- B$calcs
  cx <- gsub("\\*|\\+|-|\\/|\\(|\\)"," ",cx)
  cx <- paste(cx,collapse=" ")
  cx <- strsplit(cx,split=" ")[[1]]
  cx <- cx[cx!=""]
  cx <- cx[cx!="1"]
  ## get non calcs
  cx <- c(cx,B$vars)
  unique(cx)
}

makeTestData <- function(ncheck,vnames){
  A <- data.table(vnames,value=runif(length(vnames)))
  A <- A[rep(1:length(vnames),each=ncheck)]
  idz <- rep(1:ncheck,length(vnames))
  A[,id:=idz]
  A[,value:=runif(nrow(A))]
  dcast(A,id~vnames,value.var = 'value')
}


## checking
vrz.soc <- showAllParmz(SOC)
vrz.spd <- showAllParmz(IPD)
vrz.idh <- showAllParmz(IDH)
vrz.iph <- showAllParmz(IPH)
vrz <- c(vrz.soc,
         vrz.spd,
         vrz.idh,
         vrz.iph
         )
vrz <- unique(vrz)
A <- makeTestData(50,vrz)


## checks
IPD.F$checkfun(A) #NOTE OK
IPH.F$checkfun(A) #NOTE OK
IDH.F$checkfun(A) #NOTE OK
SOC.F$checkfun(A) #NOTE OK

## TODO port over model-runner, include model functions
## TODO add showAllParmz to HE dtree, also tester and runall

