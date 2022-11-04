## TODO questions:
## 
## - stool/sputum
## - flag assumption = groups for SA or pending data
## 

## ========= UTILITIES ===============
logit <- function(x) log(odds(x))
ilogit <- function(x) iodds(exp(x))
odds <- function(x) x/(1-x)
iodds <- function(x) x/(1+x)
lo <- function(x) quantile(x,probs = 0.025)
hi <- function(x) quantile(x,probs = 1-0.025)
rot45 <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

brkt <- function(M,L,H,ndp=0) paste0(round(M,ndp),' (',
                                     round(L,ndp),' - ',
                                     round(H,ndp),')')
gm <- function(x) exp(mean(log(x))) #geometric mean
gh <- function(x) glue(here(x))

## ========= DIAGNOSIS ===============

## function for combining sample modality with bacteriological test
## to calculate the probablility bac+ TB is diagnosed
TBbacsampletest <- function(samplepossible,testpos){
  samplepossible * testpos
}


## OR function - with NPA/stool in mind
## positive if either sample/test combo is
TBbac4ST1orST2 <- function(S1,T1,S2,T2){
  1 - (1-TBbacsampletest(S1,T1)) * (1-TBbacsampletest(S2,T2))
}


## function to add in the Sample/Test combined probabilities of TB dx
## works by side effect
AddSampleTests <- function(D){

  ## ------- X on NPA or stool/sputum ------- NOTE working with stool for now; also assuming all same
  ## IDH, at DH: TB dx bac+ on Xpert Ultra on NPA & stool/sputum, in people identified as having presumptive TB
  D[,d.idh.dh.ptbxns:=TBbac4ST1orST2(ST.poss.dh,ifelse(tb=='TB+',sens.xnpa,1-spec.xnpa),
                                  NPA.poss.dh,ifelse(tb=='TB+',sens.xstool,1-spec.xstool))]
  ## IPH, at DH: TB dx bac+ on Xpert Ultra on NPA & stool/sputum, in people identified as having presumptive TB
  D[,d.iph.dh.ptbxns:=TBbac4ST1orST2(ST.poss.dh,ifelse(tb=='TB+',sens.xnpa,1-spec.xnpa),
                                  NPA.poss.dh,ifelse(tb=='TB+',sens.xstool,1-spec.xstool))]
  ## IPH, at PHC: TB dx bac+ on Xpert Ultra on NPA & stool/sputum, in people identified as having presumptive TB
  D[,d.iph.phc.ptbxns:=TBbac4ST1orST2(ST.poss.phc,ifelse(tb=='TB+',sens.xnpa,1-spec.xnpa),
                                   NPA.poss.phc,ifelse(tb=='TB+',sens.xstool,1-spec.xstool))]
  ## IPD, at DH: TB dx bac+ on Xpert Ultra on NPA & stool/sputum, in people identified as having presumptive TB
  D[,d.ipd.dh.ptbxns:=TBbac4ST1orST2(ST.poss.dh,ifelse(tb=='TB+',sens.xnpa,1-spec.xnpa),
                                  NPA.poss.dh,ifelse(tb=='TB+',sens.xstool,1-spec.xstool))]
  ## IPD, at DH: TB dx bac+ on Xpert Ultra on NPA & stool/sputum, in people chosen to be referred from PHC who were bac- on sputum test at PHC
  D[,d.ipd.dhreftest.ptbxns:=TBbac4ST1orST2(ST.poss.dh,ifelse(tb=='TB+',sens.xnpa,1-spec.xnpa),
                                         NPA.poss.dh,ifelse(tb=='TB+',sens.xstool,1-spec.xstool))]
  ## IPD, at DH: TB dx bac+ on Xpert Ultra on NPA & stool/sputum, in people chosen to be referred from PHC who were not tested at PHC
  D[,d.ipd.dhrefnotest.ptbxns:=TBbac4ST1orST2(ST.poss.dh,ifelse(tb=='TB+',sens.xnpa,1-spec.xnpa),
                                           NPA.poss.dh,ifelse(tb=='TB+',sens.xstool,1-spec.xstool))]

  ## ------- X on sputum/GA -------
  ## SOC, at DH: people receiving Xpert Ultra testing [either sputum or GA], in those identified as having presumptive TB
  D[,d.soc.dh.test:=ifelse(age=='5-14',d.soc.dh.test.o5,d.soc.dh.test.u5)]

  ## ------- X on sputum ------- #NOTE hardcodes 0 sample availability for u5
  ## IPH, at PHC: TB dx bac+ on Xpert Ultra on sputum, in people identified as having presumptive TB
  D[,d.ipd.phc.ptbxsp:=TBbacsampletest(ifelse(age=='5-14',ES.poss.phc.o5,0),
                                    ifelse(tb=='TB+',sens.xsputum,1-spec.xsputum))]
  ## SOC, at DH: TB dx bac+ on Xpert Ultra on sputum, in people identified as having presumptive TB
  D[,d.soc.dh.ptbxsp:=TBbacsampletest(ifelse(age=='5-14',ES.poss.dh.o5,0),
                                   ifelse(tb=='TB+',sens.xsputum,1-spec.xsputum))]
  ## SOC, at PHC: TB dx bac+ on Xpert Ultra on sputum, in people identified as having presumptive TB
  D[,d.soc.phc.ptbxsp:=TBbacsampletest(ifelse(age=='5-14',ES.poss.phc.o5,0),
                                    ifelse(tb=='TB+',sens.xsputum,1-spec.xsputum))]
  ## SOC, at DH: TB dx bac+ on Xpert Ultra on sputum, in people chosen to be referred from PHC who were bac- on sputum test at PHC
  D[,d.soc.dhreftest.ptbxsp:=TBbacsampletest(ifelse(age=='5-14',ES.poss.dh.o5,0),
                                          ifelse(tb=='TB+',sens.xsputum,1-spec.xsputum))]
  ## SOC, at DH: TB dx bac+ on Xpert Ultra on sputum, in people chosen to be referred from PHC who were not tested at PHC
  D[,d.soc.dhrefnotest.ptbxsp:=TBbacsampletest(ifelse(age=='5-14',ES.poss.dh.o5,0),
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
## TODO - remove excess RNG here
CFRtxY <- function(age,hiv=0,art=0,P){#NB optimized for clarity not speed
  if(length(age)>1 & length(hiv)==1) hiv <- rep(hiv,length(age))
  if(length(age)>1 & length(art)==1) art <- rep(art,length(age))
  tmp <- P$ontx.u5$r(length(age))
  tmp[age=='5-14'] <- P$ontx.o5$r(sum(age=='5-14'))  #NB this could be achieved in  the tree model
  ## hivartOR
  Z <- P$hivartOR$r(length(age))
  hor <- rep(1,length(age))
  tmp <- logit(tmp)                     #transform
  tmp[hiv>0] <- tmp[hiv>0]+Z[hiv>0,1]
  tmp[art>0] <- tmp[art>0]+Z[art>0,2]
  tmp <- ilogit(tmp)                    #inverse transform
  tmp
}
## CFRtxY(1:10,P)                            #test
## summary(CFRtxY(1:1e3,P))
## summary(CFRtxY(1:1e3,hiv=1))
## summary(CFRtxY(1:1e3,hiv=1,art=1,P))


## == CFR off tx
CFRtxN <- function(age,hiv=0,art=0,P){
  if(length(age)>1 & length(hiv)==1) hiv <- rep(hiv,length(age))
  if(length(age)>1 & length(art)==1) art <- rep(art,length(age))
  tmp <- P$notx.u5$r(length(age))          #default a<5 and hiv=art=0
  tmp[age!='5-14' & hiv>0 & art==0] <- P$notxH.u5$r(sum(age!='5-14' & hiv>0 & art==0)) #u5,HIV+,ART-
  tmp[age!='5-14' & hiv>0 & art>0] <- P$notxHA.u5$r(sum(age!='5-14' & hiv>0 & art>0)) #u5,HIV+,ART+
  tmp[age=='5-14'] <- P$notx.o5$r(sum(age=='5-14'))    #o5, HIV-ve
  tmp[age=='5-14' & hiv>0 & art==0] <- P$notxH.o5$r(sum(age=='5-14' & hiv>0 & art==0)) #o5,HIV+,ART-
  tmp[age=='5-14' & hiv>0 & art>0] <- P$notxHA.o5$r(sum(age=='5-14' & hiv>0 & art>0)) #o5,HIV+,ART+
  tmp
}
## CFRtxN(1:10,P)                            #test
## summary(CFRtxN(1:1e3,P))
## summary(CFRtxN(1:1e3,hiv=1,P))
## summary(CFRtxN(1:1e3,hiv=1,art=1,P))




## add CFRs to data by side-effect
AddCFRs <- function(D,P){
  ## d.cfr.notx & d.cfr.tx
  D[tb=='noTB',c('d.cfr.notx','d.cfr.tx'):=0] #NOTE neglect non-TB mortality
  ## CFR on  ATT
  D[tb!="noTB",d.cfr.tx:=CFRtxY(age,hiv,art,P)]
  ## CFR w/o ATT
  D[tb!="noTB",d.cfr.notx:=CFRtxN(age,hiv,art,P)]
}


## ======= COMBINED LABELLER ===========

## additional labels from data (may overwrite some initial version currently)
AddDataDrivenLabels <- function(D){

  ## initial care seeking
  D[,d.soc.pphc:=ifelse(age=='0-4',d.soc.pphc.u5,d.soc.pphc.o5)]
  D[,d.idh.pphc:=ifelse(age=='0-4',d.idh.pphc.u5,d.idh.pphc.o5)]
  D[,d.iph.pphc:=ifelse(age=='0-4',d.iph.pphc.u5,d.iph.pphc.o5)]

  ## TB more likely to go to DH
  ## D[tb!='noTB',d.soc.pphc:=1-iodds(d.OR.dh.if.TB*odds(1-d.soc.pphc))]
  ## D[tb!='noTB',d.idh.pphc:=1-iodds(d.OR.dh.if.TB*odds(1-d.idh.pphc))]
  ## D[tb!='noTB',d.iph.pphc:=1-iodds(d.OR.dh.if.TB*odds(1-d.iph.pphc))]
  ## version with arm dependence:
  D[tb!='noTB' & age=='0-4',d.soc.pphc:=1-iodds(d.OR.dh.if.TB.soc.u5*odds(1-d.soc.pphc))]
  D[tb!='noTB' & age=='0-4',d.idh.pphc:=1-iodds(d.OR.dh.if.TB.idh.u5*odds(1-d.idh.pphc))]
  D[tb!='noTB' & age=='0-4',d.iph.pphc:=1-iodds(d.OR.dh.if.TB.iph.u5*odds(1-d.iph.pphc))]
  D[tb!='noTB' & age!='0-4',d.soc.pphc:=1-iodds(d.OR.dh.if.TB.soc.o5*odds(1-d.soc.pphc))]
  D[tb!='noTB' & age!='0-4',d.idh.pphc:=1-iodds(d.OR.dh.if.TB.idh.o5*odds(1-d.idh.pphc))]
  D[tb!='noTB' & age!='0-4',d.iph.pphc:=1-iodds(d.OR.dh.if.TB.iph.o5*odds(1-d.iph.pphc))]

  ## screening/assessment coverage
  D[,d.soc.dh.assess:=ifelse(age=='0-4',d.soc.dh.assess.u5,d.soc.dh.assess.o5)]
  D[,d.soc.phc.assess:=ifelse(age=='0-4',d.soc.phc.assess.u5,d.soc.phc.assess.o5)]
  D[,d.idh.dh.assess:=ifelse(age=='0-4',d.idh.dh.assess.u5,d.idh.dh.assess.o5)]
  D[,d.idh.phc.assess:=ifelse(age=='0-4',d.idh.phc.assess.u5,d.idh.phc.assess.o5)]
  D[,d.iph.dh.assess:=ifelse(age=='0-4',d.iph.dh.assess.u5,d.iph.dh.assess.o5)]
  D[,d.iph.phc.assess:=ifelse(age=='0-4',d.iph.phc.assess.u5,d.iph.phc.assess.o5)]

  ## ## prevalence set in make attributes
  ## d.TBprev.ICS.u5
  ## d.TBprev.ICS.o5

  ## specificity of presuming
  D[,c('d.dh.prsmptv','d.phc.prsmptv'):=1.0] #TODO remove this parameter if continuing like this
  D[age=='0-4',d.soc.dh.presumed:=ifelse(tb!='noTB',1.0,1-d.soc.dh.presumesp.u5)] #SOC
  D[age!='0-4',d.soc.dh.presumed:=ifelse(tb!='noTB',1.0,1-d.soc.dh.presumesp.o5)]
  D[age=='0-4',d.soc.phc.presumed:=ifelse(tb!='noTB',1.0,1-d.soc.phc.presumesp.u5)]
  D[age!='0-4',d.soc.phc.presumed:=ifelse(tb!='noTB',1.0,1-d.soc.phc.presumesp.o5)]
  D[age=='0-4',d.idh.dh.presumed:=ifelse(tb!='noTB',1.0,1-d.idh.dh.presumesp.u5)] #IDH
  D[age!='0-4',d.idh.dh.presumed:=ifelse(tb!='noTB',1.0,1-d.idh.dh.presumesp.o5)]
  D[age=='0-4',d.idh.phc.presumed:=ifelse(tb!='noTB',1.0,1-d.idh.phc.presumesp.u5)]
  D[age!='0-4',d.idh.phc.presumed:=ifelse(tb!='noTB',1.0,1-d.idh.phc.presumesp.o5)]
  D[age=='0-4',d.iph.dh.presumed:=ifelse(tb!='noTB',1.0,1-d.iph.dh.presumesp.u5)] #IPH
  D[age!='0-4',d.iph.dh.presumed:=ifelse(tb!='noTB',1.0,1-d.iph.dh.presumesp.o5)]
  D[age=='0-4',d.iph.phc.presumed:=ifelse(tb!='noTB',1.0,1-d.iph.phc.presumesp.u5)]
  D[age!='0-4',d.iph.phc.presumed:=ifelse(tb!='noTB',1.0,1-d.iph.phc.presumesp.o5)]

}



## combined function to add the labels to the tree prior to calculations
MakeTreeParms <- function(D,P){
  ## -- use of other functions
  AddSampleTests(D) #samples/tests
  AddCFRs(D,P) #outcomes
  ## -- other not covered above
  ## some parms that are only !=0 for older children:
  D[,d.ipd.phc.test:=ifelse(age=='5-14',d.ipd.phc.test.o5,0)]
  D[,d.soc.dh.fracsp:=ifelse(age=='5-14',d.soc.dh.fracsp.o5,0)]
  D[,d.soc.phc.test:=ifelse(age=='5-14',d.soc.phc.test.o5,0)]
  ## new labels from data
  AddDataDrivenLabels(D)
}

## ======= EPIDEMIOLOGY ===========

makeAttributes <- function(D){
    nrep <- nrow(D)
    D[,id:=1:nrep]
    fx <- list(age=agelevels,hiv=hivlevels,art=artlevels,tb=tblevels)
    cofx <- expand.grid(fx)
    cat('Attribute combinations used:\n')
    print(cofx)
    D <- D[rep(1:nrow(D),each=nrow(cofx))] #expand out data
    D[,names(cofx):=cofx[rep(1:nrow(cofx),nrep),]]
    ## --- age
    D[,value:=ifelse(age=='5-14',1-d.F.u5,d.F.u5)] #NOTE value first set here
    ## --- HIV/ART
    D[,h01:=0]
    D[age!='5-14',h10:=d.hivprev.u5*(1-d.artcov)]
    D[age=='5-14',h10:=d.hivprev.o5*(1-d.artcov)]
    D[age!='5-14',h00:=1-d.hivprev.u5]
    D[age=='5-14',h00:=1-d.hivprev.o5]
    D[age!='5-14',h11:=d.hivprev.u5*d.artcov]
    D[age=='5-14',h11:=d.hivprev.o5*d.artcov]
    D[hiv==0 & art==0,value:=value*h00]
    D[hiv==0 & art==1,value:=value*h01]
    D[hiv==1 & art==0,value:=value*h10]
    D[hiv==1 & art==1,value:=value*h11]
    D[,c('h00','h01','h10','h11'):=NULL]
    ## --- TB
    ## ## (old version) calculate true TB prev among initial care seeking as:
    ## ## tbi = f x tbp + (1-f) x tbd
    ## ## where: f=fraction initially seeking care at PHC; tbp=prev @ phc; tbd=prev @ dh
    ## ## NOTE the 'underlying' TB prev in care-seekers in controlled by soc parms
    ## D[,tbi:= d.soc.pphc * d.phc.tbinprsmptv + (1-d.soc.pphc) * d.dh.tbinprsmptv]
    ## D[tb=='noTB',value:=value*(1-tbi)]
    ## D[tb=='TB-',value:=value*tbi*ifelse(age=='5-14',1-Fbc.o5,1-Fbc.u5)] #NOTE assuming no TB outside of presumptive?
    ## D[tb=='TB+',value:=value*tbi*ifelse(age=='5-14',Fbc.o5,Fbc.u5)]
    ## D[,tbi:=NULL]                            #remove temporary variable
    ## (new version) based on data:
    D[,tbi:=ifelse(age=='5-14',d.TBprev.ICS.o5,d.TBprev.ICS.u5)]
    D[tb=='noTB',value:=value*(1-tbi)]
    D[tb=='TB-',value:=value*tbi*ifelse(age=='5-14',1-Fbc.o5,1-Fbc.u5)] #NOTE assuming no TB outside of presumptive?
    D[tb=='TB+',value:=value*tbi*ifelse(age=='5-14',Fbc.o5,Fbc.u5)]
    D[,tbi:=NULL]                            #remove temporary variable
    return(D)
}



## function for generating random sample of costs
MakeCostData <- function(csts,          #base data table of cost data
                         nrep,          #number of replicates being used in PSA
                         anmz=NULL      #attribute names (if any)
                         ){
  if(nrow(csts[cost.sd>0 & cost.m==0])>0) warning(paste0('Some cost input variables have zero mean & SD>0. These will be treated as fixed variables:\n',paste0(csts[cost.sd>0 & cost.m==0,cost],collapse='\n')))
  if(is.null(anmz)& any(csts[,table(cost)]>1)) warning('Some cost names occur >1 times, but no attributes have been specified! This is unlikely to do what you want.')
  csts[cost.m>0,gmsc:=cost.sd^2/cost.m]
  csts[!is.na(gmsc) & gmsc > 0, gmk:=cost.m/gmsc]
  NR <- nrow(csts)
  csts <- csts[rep(1:NR,nrep)]
  csts[,id:=rep(1:nrep,each=NR)]
  csts[,rnd:=!is.na(gmsc) & !is.na(gmk) & gmk>0 & gmsc > 0]
  csts[rnd==TRUE,value:=rgamma(sum(rnd),shape=gmk,scale = gmsc)] #random sample from gamma distribution
  csts[rnd!=TRUE,value:=cost.m]                                  #fixed values
  ## csts[,cnms:=paste0('c_',cost)]
  csts[,cnms:=paste0(cost)]
  F <- 'id '
  if(!is.null(anmz)) F <- paste0(F,'+ ',paste(anmz,collapse='+')) #split out by attributes if included
  F <- paste0(F, ' ~ cnms')
  dcast(csts,as.formula(F),value.var = 'value')      #id ~ cnms
}
## NOTE
## if attributes are included, all costs need to be specified by them even if this means duplicating those without dependence





## making life years
GetLifeYears <- function(isolist,discount.rate,yearfrom){
    ## template:
    LYT <- data.table(age=0:14,
                      age_group=c(rep('0-4',5),rep('5-14',10)),
                      LYS=0.0)
    ## make country/age key
    LYK <- list()
    for(iso in isolist){
        ## iso <- cn
        tmp <- copy(LYT)
        tmp[,iso3:=iso]
        for(ag in tmp$age)
            tmp[age==ag,LYS:=discly::discly(iso3=iso,
                                            age=ag,
                                            yearnow=yearfrom,
                                            sex='Total',
                                            endyear = 2098,
                                            HR=1,
                                            dr=discount.rate,
                                            hiv='both'
                                            )]
        LYK[[iso]] <- tmp
    }
    LYK <- rbindlist(LYK)
    ## assume unweighted & collapse
    LYK <- LYK[,.(LYS=mean(LYS)),by=.(iso3,age=age_group)]
    setkey(LYK,age)
    LYK
}

## ## scraps for development of below fn
## data <- merge(D,LYK,by='age') #add age
## Kmax <- 1e3
## file.id <- 'test'
## wtp <- 500

## NOTE this is more illustrative for now
## NOTE needs a folder called graphs/ creating (which is currently excluded from the repo)
## some automatic CEA outputs
MakeCEAoutputs <- function(data,LY,
                           file.id='',Kmax=5e3,wtp=5e3,
                           arms=c('SOC','IPD','IDH','IPH')){
  data <- merge(data,LY,by='age') #add age
  DS <- data[,.(cost.SOC=sum(cost.soc*value),
                cost.IPD=sum(cost.ipd*value),
                cost.IDH=sum(cost.idh*value),
                cost.IPH=sum(cost.iph*value),
                lyl.SOC=sum(deaths.soc*value*LYS),
                lyl.IPD=sum(deaths.ipd*value*LYS),
                lyl.IDH=sum(deaths.idh*value*LYS),
                lyl.IPH=sum(deaths.iph*value*LYS)),
             by=id] #PSA summary

  ## prep for BCEA
  LYS <- CST <- matrix(nrow=nreps,ncol=4)
  LYS[,1] <- 1-DS$lyl.SOC #NOTE this is life years lost
  LYS[,2] <- 1-DS$lyl.IPD
  LYS[,3] <- 1-DS$lyl.IDH
  LYS[,4] <- 1-DS$lyl.IPH
  CST[,1] <- DS$cost.SOC
  CST[,2] <- DS$cost.IPD
  CST[,3] <- DS$cost.IDH
  CST[,4] <- DS$cost.IPH
  ## BCEA outputs
  M <- bcea(e=LYS,c=CST,ref=1,interventions = arms,Kmax=Kmax)
  print(summary(M))

  fn <- paste0(here('outdata/kstar_'),file.id,'.txt')
  cat(M$kstar,file = fn)
  fn <- paste0(here('outdata/ICER_'),file.id,'.txt')
  cat(M$ICER,file = fn)

  ## NOTE may need more configuration
  ceac.plot(M,graph='ggplot2') +
    scale_x_continuous(label=comma) +
    theme_classic() + ggpubr::grids()
  fn <- paste0(here('graphs/CEAC_'),file.id,'.png')
  ggsave(file=fn,w=7,h=7)

  ceplane.plot(M,graph='ggplot2',wtp=wtp)+
    scale_x_continuous(label=comma) +
    theme_classic() +
    theme(legend.position = 'top') + ggpubr::grids()
  fn <- paste0(here('graphs/CE_'),file.id,'.png')
  ggsave(file=fn,w=7,h=7)

  eib.plot(M,graph='ggplot2',wtp=wtp) +
    scale_x_continuous(label=comma) +
    theme_classic() + ggpubr::grids()
  fn <- paste0(here('graphs/EIB_'),file.id,'.png')
  ggsave(file=fn,w=7,h=7)

  evi.plot(M,graph='ggplot2',wtp=wtp) +
    scale_x_continuous(label=comma) +
    theme_classic() + ggpubr::grids()
  fn <- paste0(here('graphs/EVI_'),file.id,'.png')
  ggsave(file=fn,w=7,h=7)

}

## --- for reformatting costs
reformatCosts <- function(rcsts){
  iextra <- outer(isoz,c('.lo','.hi','drop'),paste0)
  iextra <- c(t(iextra)); iextra <- rev(rev(iextra)[-1])
  nnmz <- c('drop','DESCRIPTION','NAME',iextra)
  names(rcsts)[1:length(nnmz)] <- nnmz
  drop <- grep('drop',names(rcsts),value=TRUE)
  rcsts[,c(drop):=NULL]
  rcsts[is.na(rcsts)] <- 0.1 #dummy
  rcsts <- melt(rcsts,id=c('NAME','DESCRIPTION'))
  rcsts[,DESCRIPTION:=NULL]
  rcsts[,c('iso3','hilo'):=tstrsplit(variable,split="\\.")]
  rcsts <- dcast(rcsts,iso3 + NAME ~ hilo,value.var = 'value')
  rcsts[,c('cost.m','cost.sd'):=.((lo+hi)/2,(hi-lo)/3.92)]
  rcsts <- rcsts[,.(iso3,cost=NAME,cost.m,cost.sd)]
  rcsts
}


## --- for computing cascades
computeCascadeData <- function(D){
  nmz <- names(D)
  rnmz0 <- grep("DH\\.|PHC\\.",nmz,value=TRUE)
  rnmz <- c('id','age','tb','value',rnmz0)
  A <- D[,..rnmz]
  ## A[,sum(value),by=id] #CHECK
  A[,pop:=value]       #rename for melting
  A[,value:=NULL]
  A[,TB:=ifelse(tb!='noTB','TB','not TB')] #simple version of TB indicator
  A[,tb:=NULL]                             #drop
  ## nrow(A) #240K (attributes x nreps)
  ## A[,sum(pop),by=id] #CHECK
  ## A[TB=='TB',1e2*sum(pop),by=id] #CHECK
  ## A[,sum(pop),by=.(id,age)] #CHECK

  ## population scaling
  A[,c(rnmz0):=lapply(.SD,function(x) x*pop),
    .SDcols=rnmz0] #multiply variables by population
  A[,pop:=NULL]                                             #can drop now

  ## melt
  AM <- melt(A,id=c('id','age','TB'))
  AM[,c('location','stage','arm'):=tstrsplit(variable,split='\\.')]
  ## sum over other attributes
  AM <- AM[,.(value=sum(value)),by=.(id,arm,age,location,stage,TB)]
  ## nrow(AM) # CHECK
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

  ## return
  list(woTB=AS2,wTB=AS)

}

## --- function for preliminary run throughts to compute mean sense/spec
computeDxAccuracy <- function(PD0,PD1,C,nreps){
  PD1[,2] <- "0.5" #not relevant but needed to generate answers
  P1 <- parse.parmtable(PD0)             #convert into parameter object
  P2 <- parse.parmtable(PD1)             #convert into parameter object
  notIPD <- c('SOC','IDH','IPH')
  ## ## NOTE for computing sense
  P2$d.TBprev.ICS.o5 <- 0.9999
  P2$d.TBprev.ICS.u5 <- 0.9999
  P <- c(P1,P2)
  cat('calculating...\n')
  D <- makePSA(nreps,P,dbls = list(c('cfrhivor','cfrartor')))
  invisible(capture.output(D <- makeAttributes(D) ))
  invisible(MakeTreeParms(D,P))
  D <- merge(D,C,by='id',all.x = TRUE)        #merge into PSA
  invisible(capture.output( D <- runallfuns(D,arm=notIPD)))
  AS2 <- computeCascadeData(D)$woTB
  tpl2 <- AS2[stage=='presumed']
  AS2 <- merge(AS2,tpl2[,.(arm,age,location,prmid=mid)],
               by=c('arm','age','location'),all.x=TRUE)
  sense <- AS2[stage=='treated',.(sense=round(1e2*mid/prmid,1),
                                  arm,age,location)]
  fwrite(sense,file=here('graphs/sense.csv'))
  cat('mean dx sensitivity computed outputed to graphs/sense.csv:\n')
  print(sense)

  ## NOTE for computing spec
  P2$d.TBprev.ICS.o5 <- 1e-5
  P2$d.TBprev.ICS.u5 <- 1e-5
  P <- c(P1,P2)
  cat('calculating...\n')
  D <- makePSA(nreps,P,dbls = list(c('cfrhivor','cfrartor')))
  invisible(capture.output( D <- makeAttributes(D) ))
  invisible(MakeTreeParms(D,P))
  D <- merge(D,C,by='id',all.x = TRUE)        #merge into PSA
  invisible(capture.output( D <- runallfuns(D,arm=notIPD)))
  AS2 <- computeCascadeData(D)$woTB
  tpl2 <- AS2[stage=='presumed']
  AS2 <- merge(AS2,tpl2[,.(arm,age,location,prmid=mid)],
               by=c('arm','age','location'),all.x=TRUE)
  spec <- AS2[stage=='treated',.(spec=round(1e2*(1-mid/prmid),1),
                                 arm,age,location)]
  fwrite(spec,file=here('graphs/spec.csv'))
  cat('mean dx specificity computed outputed to graphs/spec.csv:\n')
  print(spec)

  ## output
  merge(sense,spec,by=c('arm','age','location'))

}


## ------------- function to compute parameters derived from cascades
computeCascadeParameters <- function(DD,ICS,DxA,using='mean'){
  aal <- c('arm','age','location')

  ## prevalence of TB in presumed
  attinpr <- dcast(DD[stage %in% c('presumed','treated')],
                   arm+age+location~stage,value.var = 'vpl')
  attinpr[,attinpr:=treated/presumed]
  attinpr <- merge(attinpr,DxA,by=aal)
  attinpr[,tbinpr:=(attinpr+spec/1e2-1)/(spec/1e2+sense/1e2-1)]
  while(any(attinpr$tbinpr<0)){
    cat('** specificity problems **:\n')
    print(attinpr[tbinpr<0,.(arm,age,location,spec,tbinpr)])
    cat('halving FP rate\n')
    attinpr[tbinpr<0,spec:=100-(100-spec)/2]
    attinpr[,tbinpr:=(attinpr+spec/1e2-1)/(spec/1e2+sense/1e2-1)]
  }

  ## fraction of screened presumed
  prins <- dcast(DD[stage %in% c('presumed','screened')],
                   arm+age+location~stage,value.var = 'vpl')
  prins[,prins:=presumed/screened]
  prins <- merge(prins,attinpr[,.(arm,age,location,attinpr)],by=aal)
  prins[,tbinit:=prins * attinpr] #TB at screen/ICS
  prins[,ssp:=1-(prins-tbinit)/(1-tbinit)] #specificity of screen
  ssps <- prins[,.(arm,age,location,ssp)]
  ssps[,NAME:=paste0('d.',
                     arm,'.',
                     tolower(location),'.',
                     'presumesp.',
                     ifelse(age=='0-4','u5','o5')
                     )]
  ## age split ICS
  tmp <- ICS[age=='0-4',.(p=sum(icsp)),by=arm] #similar by arm
  d.F.u5 <- tmp[,mean(p)]
  agesp <- data.table(NAME='d.F.u5',DISTRIBUTION=d.F.u5)
  ICS[,agetotal:=sum(icsp),by=.(arm,age)]
  agedest <- ICS[location=='PHC',.(pphc=icsp/agetotal),by=.(arm,age)]
  agedest[,NAME:=paste0('d.',
                        arm,'.',
                        'pphc.',
                        ifelse(age=='0-4','u5','o5')
                        )]

  alldest <- rbind(agedest[,.(arm,age,location='PHC',pinit=pphc)],
                   agedest[,.(arm,age,location='DH',pinit=1-pphc)])
  alldest <- merge(alldest,prins[,.(arm,age,location,tbinit)],by=aal)
  tbagearm <- alldest[,.(tbICS=sum(pinit*tbinit)),by=.(arm,age)]
  cat('Implied prevalences at ICS (using ',using,'):\n')
  print(tbagearm[order(age),.(`TB prev @ ICE`=1e5*tbICS),by=.(age,arm)])
  if(using=='min')
    tbbyage <- tbagearm[,.(tbICS=min(tbICS)),by=age]
  if(using=='max')
    tbbyage <- tbagearm[,.(tbICS=max(tbICS)),by=age]
  if(using=='mean')
    tbbyage <- tbagearm[,.(tbICS=mean(tbICS)),by=age]
  if(using=='median')
    tbbyage <- tbagearm[,.(tbICS=median(tbICS)),by=age]
  if(using=='gm')
    tbbyage <- tbagearm[,.(tbICS=gm(tbICS)),by=age]
  tbbyage[,NAME:=paste0("d.TBprev.ICS.",ifelse(age=='0-4','u5','o5'))]
  print(tbbyage[,.(NAME,age,tbICS)])

  ## OR for DH calculations
  ORcalc <- dcast(data=alldest,arm+age~location,value.var = c('pinit','tbinit'))
  ORcalc <- merge(ORcalc,tbagearm,by=c('arm','age'))
  ORcalc[,probPHCifTB:=tbinit_PHC * pinit_PHC / tbICS]
  ORcalc[,probPHCifNotTB:=(1-tbinit_PHC) * pinit_PHC / (1-tbICS)]
  ORcalc[,ORDH:=(1/probPHCifTB-1)/(1/probPHCifNotTB-1)]
  cat('Implied OR for CS @ DH if TB+ :\n')
  print(ORcalc[,.(arm,age,ORDH)]) #NOTE could make arm dependent
  ORcalc[,NAME:=paste0('d.OR.dh.if.TB.',
                       arm,'.',
                       ifelse(age=='0-4','u5','o5')
                       )]
  ## d.OR.dh.if.TB <- ORcalc[,mean(ORDH)] #single OR

  ## and assessent coverage
  assess <- DD[stage %in% 'screened',.(assess=vpl/1e5,arm,age,location)]
  assess[,NAME:=paste0('d.',
                     arm,'.',
                     tolower(location),'.',
                     'assess.',
                     ifelse(age=='0-4','u5','o5')
                     )]


  ## combine answer
  ## ssps,assess,agesp,agedest,tbbyage,d.OR.dh.if.TB
  NP <- rbindlist(list(
    ssps[,.(NAME,DISTRIBUTION=ssp)],
    assess[,.(NAME,DISTRIBUTION=assess)],
    agesp,
    agedest[,.(NAME,DISTRIBUTION=pphc)],
    tbbyage[,.(NAME,DISTRIBUTION=tbICS)],
    ORcalc[,.(NAME,DISTRIBUTION=ORDH)]## new version with arm dependence
    ## data.table(NAME='d.OR.dh.if.TB',DISTRIBUTION=d.OR.dh.if.TB)
  ))

  cat('saving calculated parameters to indata/calcparms_',using,'.csv\n')
  fn <- paste0(here('indata/calcparms_'),using,'.csv')
  fwrite(NP,file=fn)

  ## return value
  as.data.frame(NP)
}


## calculate mid/lo/hi
MLH <- function(dat){
  nnmz <- names(dat)
  lnmz <- paste0(nnmz,'.lo')
  hnmz <- paste0(nnmz,'.hi')
  mnmz <- paste0(nnmz,'.mid')
  L <- dat[,lapply(.SD,lo),.SDcols=nnmz]
  M <- dat[,lapply(.SD,mean),.SDcols=nnmz]
  H <- dat[,lapply(.SD,hi),.SDcols=nnmz]
  setnames(L,nnmz,lnmz); setnames(M,nnmz,mnmz); setnames(H,nnmz,hnmz);
  list(L=L,M=M,H=H)
}
## MLH(out[,.(DcostperATT.iph,DcostperATT.soc)]) #test


## =========== output formatters
outsummary <- function(out){

  ## mid/lo/hi
  outa <- MLH(out[,.(DcostperATT.iph,DcostperATT.idh,DcostperATT.soc,
                     Ddeaths.iph,Ddeaths.idh,
                     DLYL.iph,DLYL.idh,
                     DLYL0.iph,DLYL0.idh,
                     DcostperLYS0.iph,DcostperLYS0.idh,
                     DcostperLYS.iph,DcostperLYS.idh,
                     Dcostperdeaths.iph,Dcostperdeaths.idh,
                     Dcost.iph,Dcost.idh)])

  ## more bespoke statistics
  outi <- out[,.(ICER.iph= -mean(Dcost.iph) / mean(DLYL.iph),
                 ICER.idh= -mean(Dcost.idh) / mean(DLYL.idh))]

  ## join
  outs <- do.call(cbind,list(outa$M,outa$L,outa$H,outi)) #combine

  ## pretty version
  pouts <- outs[,.(DcostperATT.iph = brkt(DcostperATT.iph.mid,DcostperATT.iph.lo,DcostperATT.iph.hi),
                  DcostperATT.idh = brkt(DcostperATT.idh.mid,DcostperATT.idh.lo,DcostperATT.idh.hi),
                  DcostperLYS0.iph = brkt(DcostperLYS0.iph.mid,DcostperLYS0.iph.lo,DcostperLYS0.iph.hi),
                  DcostperLYS0.idh = brkt(DcostperLYS0.idh.mid,DcostperLYS0.idh.lo,DcostperLYS0.idh.hi),
                  DcostperLYS.iph = brkt(DcostperLYS.iph.mid,DcostperLYS.iph.lo,DcostperLYS.iph.hi),
                  DcostperLYS.idh = brkt(DcostperLYS.idh.mid,DcostperLYS.idh.lo,DcostperLYS.idh.hi),
                  Dcostperdeaths.iph = brkt(Dcostperdeaths.iph.mid,
                                            Dcostperdeaths.iph.lo,Dcostperdeaths.iph.hi),
                  Dcostperdeaths.idh = brkt(Dcostperdeaths.idh.mid,
                                            Dcostperdeaths.idh.lo,Dcostperdeaths.idh.hi),
                  DcostPerOPD.iph = brkt(Dcost.iph.mid,Dcost.iph.lo,Dcost.iph.hi),
                  DcostPerOPD.idh = brkt(Dcost.idh.mid,Dcost.idh.lo,Dcost.idh.hi),
                  DdeathsPer100kOPD.iph = brkt(-1e5*Ddeaths.iph.mid,
                                               -1e5*Ddeaths.iph.hi,-1e5*Ddeaths.iph.lo),
                  DdeathsPer100kOPD.idh = brkt(-1e5*Ddeaths.idh.mid,
                                               -1e5*Ddeaths.idh.hi,-1e5*Ddeaths.idh.lo),
                  DLYS0Per100kOPD.iph = brkt(-1e5*DLYL0.iph.mid,
                                               -1e5*DLYL0.iph.hi,-1e5*DLYL0.iph.lo),
                  DLYS0Per100kOPD.idh = brkt(-1e5*DLYL0.idh.mid,
                                            -1e5*DLYL0.idh.hi,-1e5*DLYL0.idh.lo),
                  DLYSPer100kOPD.iph = brkt(-1e5*DLYL.iph.mid,
                                            -1e5*DLYL.iph.hi,-1e5*DLYL.iph.lo),
                  DLYSPer100kOPD.idh = brkt(-1e5*DLYL.idh.mid,
                                            -1e5*DLYL.idh.hi,-1e5*DLYL.idh.lo),
                  ICER.iph=round(ICER.iph,0),ICER.idh=round(ICER.idh,0))]

  ## return value
  list(outs=outs,pouts=pouts)
}


## ---- utilities for making CEACs
make.ceac <- function(CEA,lamz){
    crv <- lamz
    for(i in 1:length(crv)) crv[i] <- CEA[,mean(lamz[i]*Q-P>0)]
    crv
}



## ## ============ experiments with reweighting


## getWeights <- function(X,V=1,W=rep(1,4)){
##   lnmz <- c('id','age','value',rnmz0)
##   XL <- X[,..lnmz]
##   XL[,c(rnmz0):=lapply(.SD,function(x) x*value),
##     .SDcols=rnmz0] #multiply variables by population
##   XL[,value:=NULL]                                             #can drop now
##   XL <- XL[,lapply(.SD,sum),.SDcols=rnmz0,by=.(id,age)]
##   XL <- merge(XL,DDW,by='age',all.x=TRUE)
##   if(TRUE){
##     ## SSE error
##     XL[,SSE:=
##           W[1]*(DH.presented.idh-DH_presented_idh/1e5)^2+
##           W[1]*(DH.presented.iph-DH_presented_iph/1e5)^2+
##           W[1]*(DH.presented.soc-DH_presented_soc/1e5)^2+
##           W[2]*(DH.presumed.idh-DH_presumed_idh/1e5)^2+
##           W[2]*(DH.presumed.iph-DH_presumed_iph/1e5)^2+
##           W[2]*(DH.presumed.soc-DH_presumed_soc/1e5)^2+
##           W[3]*(DH.screened.idh-DH_screened_idh/1e5)^2+
##           W[3]*(DH.screened.iph-DH_screened_iph/1e5)^2+
##           W[3]*(DH.screened.soc-DH_screened_soc/1e5)^2+
##           W[4]*(DH.treated.idh-DH_treated_idh/1e5)^2+
##           W[4]*(DH.treated.iph-DH_treated_iph/1e5)^2+
##           W[4]*(DH.treated.soc-DH_treated_soc/1e5)^2+
##           W[1]*(PHC.presented.idh-PHC_presented_idh/1e5)^2+
##           W[1]*(PHC.presented.iph-PHC_presented_iph/1e5)^2+
##           W[1]*(PHC.presented.soc-PHC_presented_soc/1e5)^2+
##           W[2]*(PHC.presumed.idh-PHC_presumed_idh/1e5)^2+
##           W[2]*(PHC.presumed.iph-PHC_presumed_iph/1e5)^2+
##           W[2]*(PHC.presumed.soc-PHC_presumed_soc/1e5)^2+
##           W[3]*(PHC.screened.idh-PHC_screened_idh/1e5)^2+
##           W[3]*(PHC.screened.iph-PHC_screened_iph/1e5)^2+
##           W[3]*(PHC.screened.soc-PHC_screened_soc/1e5)^2+
##           W[4]*(PHC.treated.idh-PHC_treated_idh/1e5)^2+
##           W[4]*(PHC.treated.iph-PHC_treated_iph/1e5)^2+
##           W[4]*(PHC.treated.soc-PHC_treated_soc/1e5)^2]
##     XL <- XL[,.(LL=-sum(SSE)/V),by=id]
##     XL[,wts:=LL-max(LL)]
##     XL[,wts:=exp(wts)]
##     XL[,wts:=wts/sum(wts)]
##     cat('ESS=',XL[,1/sum(wts^2)],'\n')
##     }
##   XL
## }





## ## tmp <- getWeights(D,V=1e-6,W=c(0,0,0,1))


## ## tmp[,.(id,age,
## ##   (DH.presented.idh-DH_presented_idh/1e5),
## ##   (DH.presented.iph-DH_presented_iph/1e5),
## ##   (DH.presented.soc-DH_presented_soc/1e5),
## ##   (DH.presumed.idh-DH_presumed_idh/1e5),
## ##   (DH.presumed.iph-DH_presumed_iph/1e5),
## ##   (DH.presumed.soc-DH_presumed_soc/1e5),
## ##   (DH.screened.idh-DH_screened_idh/1e5),
## ##   (DH.screened.iph-DH_screened_iph/1e5),
## ##   (DH.screened.soc-DH_screened_soc/1e5),
## ##   (DH.treated.idh-DH_treated_idh/1e5),
## ##   (DH.treated.iph-DH_treated_iph/1e5),
## ##   (DH.treated.soc-DH_treated_soc/1e5),
## ##   (PHC.presented.idh-PHC_presented_idh/1e5),
## ##   (PHC.presented.iph-PHC_presented_iph/1e5),
## ##   (PHC.presented.soc-PHC_presented_soc/1e5),
## ##   (PHC.presumed.idh-PHC_presumed_idh/1e5),
## ##   (PHC.presumed.iph-PHC_presumed_iph/1e5),
## ##   (PHC.presumed.soc-PHC_presumed_soc/1e5),
## ##   (PHC.screened.idh-PHC_screened_idh/1e5),
## ##   (PHC.screened.iph-PHC_screened_iph/1e5),
## ##   (PHC.screened.soc-PHC_screened_soc/1e5),
## ##   (PHC.treated.idh-PHC_treated_idh/1e5),
## ##   (PHC.treated.iph-PHC_treated_iph/1e5),
## ##   (PHC.treated.soc-PHC_treated_soc/1e5))]

