## TODO questions:
## - NPA and Stool poss
## - check this is or
## - costing for this
## - stool/sputum
## - flag assumption = groups for SA or pending data
## 

## ========= UTILITIES ===============
logit <- function(x) log(odds(x))
ilogit <- function(x) iodds(exp(x))
odds <- function(x) x/(1-x)
iodds <- function(x) x/(1+x)


## ========= DIAGNOSIS ===============

## function for combining sample modality with bacteriological test
## to calculate the probablility bac+ TB is diagnosed
TBbacsampletest <- function(samplepossible,testpos){
  samplepossible * testpos
}


## OR function - with NPA/stool in mind
## positive if either sample/test combo is
TBbac4ST1orST2 <- function(S1,T1,S2,T2){
  1 - TBbacsampletest(S1,T1) * TBbacsampletest(S2,T2)
}


## function to add in the Sample/Test combined probabilities of TB dx
## works by side effect
AddSampleTests <- function(D){

  ## ------- X on NPA or stool/sputum ------- NOTE working with stool for now; also assuming all same
  ## IDH, at DH: TB dx bac+ on Xpert Ultra on NPA & stool/sputum, in people identified as having presumptive TB
  D[,d.idh.dh.ptbxns:=TBbac4ST1orST2(ST.poss,ifelse(tb=='TB+',sens.xnpa,1-spec.xnpa),
                                  NPA.poss,ifelse(tb=='TB+',sens.xstool,1-spec.xstool))]
  ## IPH, at DH: TB dx bac+ on Xpert Ultra on NPA & stool/sputum, in people identified as having presumptive TB
  D[,d.iph.dh.ptbxns:=TBbac4ST1orST2(ST.poss,ifelse(tb=='TB+',sens.xnpa,1-spec.xnpa),
                                  NPA.poss,ifelse(tb=='TB+',sens.xstool,1-spec.xstool))]
  ## IPH, at PHC: TB dx bac+ on Xpert Ultra on NPA & stool/sputum, in people identified as having presumptive TB
  D[,d.iph.phc.ptbxns:=TBbac4ST1orST2(ST.poss,ifelse(tb=='TB+',sens.xnpa,1-spec.xnpa),
                                   NPA.poss,ifelse(tb=='TB+',sens.xstool,1-spec.xstool))]
  ## IPD, at DH: TB dx bac+ on Xpert Ultra on NPA & stool/sputum, in people identified as having presumptive TB
  D[,d.ipd.dh.ptbxns:=TBbac4ST1orST2(ST.poss,ifelse(tb=='TB+',sens.xnpa,1-spec.xnpa),
                                  NPA.poss,ifelse(tb=='TB+',sens.xstool,1-spec.xstool))]
  ## IPD, at DH: TB dx bac+ on Xpert Ultra on NPA & stool/sputum, in people chosen to be referred from PHC who were bac- on sputum test at PHC
  D[,d.ipd.dhreftest.ptbxns:=TBbac4ST1orST2(ST.poss,ifelse(tb=='TB+',sens.xnpa,1-spec.xnpa),
                                         NPA.poss,ifelse(tb=='TB+',sens.xstool,1-spec.xstool))]
  ## IPD, at DH: TB dx bac+ on Xpert Ultra on NPA & stool/sputum, in people chosen to be referred from PHC who were not tested at PHC
  D[,d.ipd.dhrefnotest.ptbxns:=TBbac4ST1orST2(ST.poss,ifelse(tb=='TB+',sens.xnpa,1-spec.xnpa),
                                           NPA.poss,ifelse(tb=='TB+',sens.xstool,1-spec.xstool))]

  ## ------- X on sputum ------- #NOTE hardcodes 0 sample availability for u5
  ## IPH, at PHC: TB dx bac+ on Xpert Ultra on sputum, in people identified as having presumptive TB
  D[,d.ipd.phc.ptbxsp:=TBbacsampletest(ifelse(age=='5-14',ES.poss.phc,0),
                                    ifelse(tb=='TB+',sens.xsputum,1-spec.xsputum))]
  ## SOC, at DH: TB dx bac+ on Xpert Ultra on sputum, in people identified as having presumptive TB
  D[,d.soc.dh.ptbxsp:=TBbacsampletest(ifelse(age=='5-14',ES.poss.dh,0),
                                   ifelse(tb=='TB+',sens.xsputum,1-spec.xsputum))]
  ## SOC, at PHC: TB dx bac+ on Xpert Ultra on sputum, in people identified as having presumptive TB
  D[,d.soc.phc.ptbxsp:=TBbacsampletest(ifelse(age=='5-14',ES.poss.phc,0),
                                    ifelse(tb=='TB+',sens.xsputum,1-spec.xsputum))]
  ## SOC, at DH: TB dx bac+ on Xpert Ultra on sputum, in people chosen to be referred from PHC who were bac- on sputum test at PHC
  D[,d.soc.dhreftest.ptbxsp:=TBbacsampletest(ifelse(age=='5-14',ES.poss.dh,0),
                                          ifelse(tb=='TB+',sens.xsputum,1-spec.xsputum))]
  ## SOC, at DH: TB dx bac+ on Xpert Ultra on sputum, in people chosen to be referred from PHC who were not tested at PHC
  D[,d.soc.dhrefnotest.ptbxsp:=TBbacsampletest(ifelse(age=='5-14',ES.poss.dh,0),
                                            ifelse(tb=='TB+',sens.xsputum,1-spec.xsputum))]

  ## ------- X on GA -------
  ## SOC, at DH: TB dx bac+ on Xpert Ultra on GA, in people identified as having presumptive TB
  D[,d.soc.dh.ptbxga:=TBbacsampletest(ifelse(age=='5-14',GA.poss.dh.o5,GA.poss.dh.u5),
                                   ifelse(tb=='TB+',sens.xga,1-spec.xga))]
  ## SOC, at DH: TB dx bac+ on Xpert Ultra on GA, in people chosen to be referred from PHC who were bac- on sputum test at PHC
  D[,d.soc.dhreftest.ptbxga:=TBbacsampletest(ifelse(age=='5-14',GA.poss.dh.o5,GA.poss.dh.u5),
                                          ifelse(tb=='TB+',sens.xga,1-spec.xga))]
  ## SOC, at DH: TB dx bac+ on Xpert Ultra on GA, in people chosen to be referred from PHC who were not tested at PHC
  D[,d.soc.dhrefnotest.ptbxga:=TBbacsampletest(ifelse(age=='5-14',GA.poss.dh.o5,GA.poss.dh.u5),
                                            ifelse(tb=='TB+',sens.xga,1-spec.xga))]

  ## ------- clinical -------
  ## IDH, at DH: TB dx clinical, in bac- people identified as having presumptive TB
  D[,d.idh.dh.test.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]
  ## IPH, at DH: TB dx clinical, in bac- people identified as having presumptive TB
  D[,d.iph.dh.test.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]
  ## IPH, at PHC: TB dx clinical, in bac- people identified as having presumptive TB
  D[,d.iph.phc.test.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]
  ## IPD, at DH: TB dx clinical, in bac- people identified as having presumptive TB
  D[,d.ipd.dh.test.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]
  ## IPD, at PHC: TB dx clinical, in bac- people identified as having presumptive TB
  D[,d.ipd.phc.test.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]
  ## SOC, at DH: TB dx clinical, in bac- people identified as having presumptive TB
  D[,d.soc.dh.test.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]
  ## SOC, at PHC: TB dx clinical, in bac- people identified as having presumptive TB
  D[,d.soc.phc.test.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]
  ## IPD, at PHC: TB dx clinical, in untested people identified as having presumptive TB
  D[,d.ipd.phc.notest.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]
  ## SOC, at DH: TB dx clinical, in untested people identified as having presumptive TB
  D[,d.soc.dh.notest.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]
  ## SOC, at PHC: TB dx clinical, in untested people identified as having presumptive TB
  D[,d.soc.phc.notest.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]
  ## IPD, at DH: TB dx clinical, in bac- people chosen to be referred from PHC who were bac- on sputum test at PHC
  D[,d.ipd.dhreftest.test.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]
  ## SOC, at DH: TB dx clinical, in bac- people chosen to be referred from PHC who were bac- on sputum test at PHC
  D[,d.soc.dhreftest.test.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]
  ## SOC, at DH: TB dx clinical, in untested people chosen to be referred from PHC who were bac- on sputum test at PHC
  D[,d.soc.dhreftest.notest.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]
  ## IDH, at DH: TB dx clinical, in bac- people chosen to be referred from PHC who were not tested at PHC
  D[,d.ipd.dhrefnotest.test.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]
  ## SOC, at DH: TB dx clinical, in bac- people chosen to be referred from PHC who were not tested at PHC
  D[,d.soc.dhrefnotest.test.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]
  ## SOC, at DH: TB dx clinical, in untested people chosen to be referred from PHC who were not tested at PHC either
  D[,d.soc.dhrefnotest.notest.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]
  ## IDH, at DH: TB dx clinical at 7-day reassessment, in people bac- clin- in previous DH testing
  D[,d.idh.dh.test7.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]
  ## IPH, at DH: TB dx clinical at 7-day reassessment, in people bac- clin- in previous DH testing
  D[,d.iph.dh.test7.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]
  ## IPD, at DH: TB dx clinical at 7-day reassessment, in people bac- clin- in previous DH testing
  D[,d.ipd.dh.test7.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]
  ## IPH, at DH: TB dx clinical following CXR, referred to DH due to ongoing symptoms at 7-day reassessment following bac- clin- test at PHC
  D[,d.iph.dhreftest7.ptbcxr:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]

}


## ========= OUTCOMES ===============
CFRtxY <- function(age,hiv=0,art=0){#NB optimized for clarity not speed
  if(length(age)>1 & length(hiv)==1) hiv <- rep(hiv,length(age))
  if(length(age)>1 & length(art)==1) art <- rep(art,length(age))
  tmp <- PZ$ontxY$r(length(age))
  tmp[age=='5-14'] <- PZ$ontxO$r(sum(age=='5-14'))  #NB this could be achieved in  the tree model
  ## hivartOR
  Z <- PZ$hivartOR$r(length(age))
  hor <- rep(1,length(age))
  tmp <- logit(tmp)                     #transform
  tmp[hiv>0] <- tmp[hiv>0]+Z[hiv>0,1]
  tmp[art>0] <- tmp[art>0]+Z[art>0,2]
  tmp <- ilogit(tmp)                    #inverse transform
  tmp
}
## CFRtxY(1:10)                            #test
## summary(CFRtxY(1:1e3))
## summary(CFRtxY(1:1e3,hiv=1))
## summary(CFRtxY(1:1e3,hiv=1,art=1))


## == CFR off tx
CFRtxN <- function(age,hiv=0,art=0){
  if(length(age)>1 & length(hiv)==1) hiv <- rep(hiv,length(age))
  if(length(age)>1 & length(art)==1) art <- rep(art,length(age))
  tmp <- PZ$notxY$r(length(age))          #default a<5 and hiv=art=0
  tmp[age!='5-14' & hiv>0 & art==0] <- PZ$notxHY$r(sum(age!='5-14' & hiv>0 & art==0)) #u5,HIV+,ART-
  tmp[age!='5-14' & hiv>0 & art>0] <- PZ$notxHAY$r(sum(age!='5-14' & hiv>0 & art>0)) #u5,HIV+,ART+
  tmp[age=='5-14'] <- PZ$notxO$r(sum(age=='5-14'))    #o5, HIV-ve
  tmp[age=='5-14' & hiv>0 & art==0] <- PZ$notxHO$r(sum(age=='5-14' & hiv>0 & art==0)) #o5,HIV+,ART-
  tmp[age=='5-14' & hiv>0 & art>0] <- PZ$notxHAO$r(sum(age=='5-14' & hiv>0 & art>0)) #o5,HIV+,ART+
  tmp
}
## CFRtxN(1:10)                            #test
## summary(CFRtxN(1:1e3))
## summary(CFRtxN(1:1e3,hiv=1))
## summary(CFRtxN(1:1e3,hiv=1,art=1))




## add CFRs to data by side-effect
AddCFRs <- function(D){
  ## d.cfr.notx & d.cfr.tx
  D[tb=='noTB',c('d.cfr.notx','d.cfr.tx'):=0] #NOTE neglect non-TB mortality
  ## CFR on  ATT
  D[tb!="noTB",d.cfr.tx:=CFRtxY(age,hiv,art)]
  ## CFR w/o ATT
  D[tb!="noTB",p.cfr.notx:=CFRtxN(age,hiv,art)]

}


## ======= COMBINED LABELLER ===========

## combined function to add the labeks to the tree prior to calculations
MakeTreeParms <- function(D){
  ## -- use of other functions
  AddSampleTests(D) #samples/tests
  AddCFRs(D) #outcomes
  ## -- other not covered above
  ## some parms that are only !=0 for older children:
  D[,d.ipd.phc.test:=ifelse(age=='5-14',d.ipd.phc.test.514,0)]
  D[,d.soc.dh.fracsp:=ifelse(age=='5-14',d.soc.dh.fracsp.514,0)]
  D[,d.soc.phc.test:=ifelse(age=='5-14',d.soc.phc.test.514,0)]
}

## ======= EPIDEMIOLOGY ===========
## TODO - HIV/ART IRRs
makeAttributes <- function(D){
    nrep <- nrow(D)
    D[,id:=1:nrep]
    fx <- list(age=agelevels,tb=tblevels,hiv=hivlevels,art=artlevels)
    cofx <- expand.grid(fx)
    cat('Attribute combinations used:\n')
    print(cofx)
    D <- D[rep(1:nrow(D),each=nrow(cofx))] #expand out data
    D[,names(cofx):=cofx[rep(1:nrow(cofx),nrep),]]
    ## age
    D[,value:=ifelse(age=='5-14',1-d.F.u5,d.F.u5)] #NOTE value first set here
    ## HIV/ART
    D[age!='5-14',value:=value*ifelse(hiv==1,d.hivprev.u5,d.hivprev.u5)]
    D[age=='5-14',value:=value*ifelse(hiv==1,d.hivprev.o5,d.hivprev.u5)]
    D[hiv==1,value:=value*ifelse(art==1,d.artcov,1-d.artcov)]
    D[hiv==0 & art==1,value:=0]
    ## TB
    D[tb=='noTB',value:=value*d.dh.tbinprsmptv]
    D[tb=='TB-',value:=value*d.dh.tbinprsmptv*(1-Fbc.u5)] #NOTE assuming no TB outside of presumptive?
    D[tb=='TB+',value:=value*d.dh.tbinprsmptv*(Fbc.u5)]   #TODO need value swapouts to handle different cohorts
    return(D)
}

## TODO bring runner from tree into this file to build in swap outs
