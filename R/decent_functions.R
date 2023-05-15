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

brkt2 <- function(M,L,H,ndp=0) paste0(format(round(M,ndp),big.mark = ','),' (',
                                      format(round(L,ndp),big.mark = ','),' to ',
                                      format(round(H,ndp),big.mark = ','),')')


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
  rcsts[is.na(rcsts)] <- 0.0
  ## drop the extra columns if existing
  if(ncol(rcsts)>(which(names(rcsts)=='ZMB.hi')+1)){
    drop <- (which(names(rcsts)=='ZMB.hi')+1):ncol(rcsts)
    drop <- names(rcsts)[drop]
    rcsts[,c(drop):=NULL]
  }
  rcsts <- melt(rcsts,id=c('NAME','DESCRIPTION'))
  rcsts[,DESCRIPTION:=NULL]
  rcsts[,c('iso3','hilo'):=tstrsplit(variable,split="\\.")]
  rcsts <- rcsts[NAME!=''] #drop extra rows if existing
  rcsts <- dcast(rcsts,iso3 + NAME ~ hilo,value.var = 'value')
  rcsts[,c('cost.m','cost.sd'):=.((lo+hi)/2,(hi-lo)/3.92)]
  rcsts <- rcsts[,.(iso3,cost=NAME,cost.m,cost.sd)]
  rcsts
}


## --- sampling spec and tbprev condiotional on cascade (by sideeffect)
## x = prev; y = phi
## NOTE this requires priors for phi/tbprev to be globally defined
## NOTE this needs the DoC cascade parameter, as well as calculated sense, a, b
BayesSpecPrev <- function(DI){
  ## to generate conditional sample, need treated/presumed
  ## Y = l - (se-l) * (a+b*x)
  ## Y = [l - (se-l) * a]  - (se-l) * b* x
  ## exp[ -(x-mx)^2/(2*Sx^2)) -(y-my)^2/(2*Sy^2)) ] = 
  ## exp[ -(x-mx)^2/(2*Sx^2)) -([l - (se-l)*a]-(se-l)*b*x-my)^2/(2*Sy^2)) ] =
  ## exp[ -(x-mx)^2/(2*Sx^2)) -(x-(l-(se-l)*a-my)/((se-l)*b))^2/(2*Sy^2/((se-l)*b)^2) ] =
  ## 1/V = 1/Sx^2 + ((se-l)*b)^2/Sy^2
  ## M =  { mx/Sx^2 + [(l-(se-l)*a-my)*((se-l)*b)]/Sy^2 } * V
  DI[,V:=1/(1/prior.tbprev$sd^2 + ((sense-DoC)/b)^2/prior.phi$sd^2)]
  DI[,M:=V * (prior.tbprev$mean/prior.tbprev$sd^2 +
                (DoC-a-b*prior.phi$mean)*(sense-DoC)/(b^2*prior.phi$sd^2))]

  DI[,tbprev:=rnorm(nrow(DI),mean=M,sd=sqrt(V))]
  DI[,phi:=abs((DoC-a)/b-((sense-DoC)/b)*tbprev)]
  DI[,tbprev:=iodds(abs(tbprev))]
}


## --- for computing cascades
computeCascadeData <- function(D,PSA=FALSE){
  nmz <- names(D)
  rnmz0 <- grep("DH\\.|PHC\\.",nmz,value=TRUE)
  rnmz0 <- rnmz0[!grepl('cost',rnmz0)] #don't include the cost-by-stage data
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
  ## version with PSA id
  if(PSA==TRUE) return(AM[,.(value=sum(value)),by=.(id,arm,age,location,stage)]) 

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


computeDxAccuracy4PSA <- function(DI){
  cat('calculating sense...\n')
  ## ## NOTE for computing sense
  DI$d.TBprev.ICS.o5 <- DI$d.TBprev.ICS.u5 <- 0.9999
  invisible(capture.output(D <- makeAttributes(DI) ))
  invisible(MakeTreeParms(D,P))
  D <- merge(D,C,by='id',all.x = TRUE)        #merge into PSA
  invisible(capture.output( D <- runallfuns(D,arm=notIPD)))
  AM <- computeCascadeData(D,PSA=TRUE)
  AM <- dcast(AM,id+arm+age+location~stage,value.var = 'value')
  ANS <- AM[,.(id,arm,age,location,sense=treated/presumed)]

  cat('calculating spec...\n') #computes spec with given parameter vals
  DI$d.TBprev.ICS.o5 <- DI$d.TBprev.ICS.u5 <- 1e-5
  invisible(capture.output( D <- makeAttributes(DI) ))
  invisible(MakeTreeParms(D,P))
  D <- merge(D,C,by='id',all.x = TRUE)        #merge into PSA
  invisible(capture.output( D <- runallfuns(D,arm=notIPD)))
  AM <- computeCascadeData(D,PSA=TRUE)
  AM <- dcast(AM,id+arm+age+location~stage,value.var = 'value')
  ANS2 <- AM[,.(id,arm,age,location,spec=1-treated/presumed)]

  cat('calculating spec1...\n') #computes spec with raw parm=1
  DI$d.TBprev.ICS.o5 <- DI$d.TBprev.ICS.u5 <- 1e-5
  DI$spec.clin <- DI$spec.clinCXR.soc <- 1-1e-5
  invisible(capture.output( D <- makeAttributes(DI) ))
  invisible(MakeTreeParms(D,P))
  D <- merge(D,C,by='id',all.x = TRUE)        #merge into PSA
  invisible(capture.output( D <- runallfuns(D,arm=notIPD)))
  AM <- computeCascadeData(D,PSA=TRUE)
  AM <- dcast(AM,id+arm+age+location~stage,value.var = 'value')
  ANS3 <- AM[,.(id,arm,age,location,spec1=1-treated/presumed)]

  ## output
  ANS <- merge(ANS,ANS2,by=c('id','arm','age','location'))
  ANS <- merge(ANS,ANS3,by=c('id','arm','age','location'))
  ANS
}

## for reading the cascade data meta-analyses
cd.read <- function(f){
  fn <- here(glue('graphs/cascades/data/{f}.csv'))
  tmp <- fread(fn)[location=='All']
  if('alpha' %in% names(tmp)){ #for alpha
    setnames(tmp,old=c('alpha','alpha.lo','alpha.hi'),new=c('mid','lo','hi'))
  }
  if('Age' %in% names(tmp)){ #for T/P = D/C
    tmp <- tmp[Arm=='IPH']   #NOTE we use IPH for TB prev calculating
    ans <- tmp[,.(nm=paste0('DoC_',Age),mid,lo,hi)]
  } else{
    ans <- tmp[,.(nm=f,mid,lo,hi)]
  }
  ans
}

## F_refs             : referrals
## F_refu             : refusals=ltfu
## F_presume          : C/B
## F_ICS              : corrected initial care seeking
## F_omega_flat       : Bd/Bp ratio see notes
## F_alpha            : INT coverage alpha
## F_SOC_BoA_modelled : SOC coverage alpha
## F_treat            : not needed

## read and parse M/A cascade data
cdnmz <- c('F_refs','F_refu','F_presume','F_ICS','F_omega_flat',
           'F_alpha','F_SOC_BoA_modelled',## 'F_treat',
           'F_presume',
           'TPS' #NOTE odd one out as Arm x Age
           )
CDL <- list()
for(nm in cdnmz) CDL[[nm]] <- cd.read(nm)
CDD <- rbindlist(CDL) #CascaDe Data


## TB-indept cascade things: alpha, omega, Fu5, pphc, lambda, rho
TBinptCascadeParms <- function(CDD,ICS){
  ## ICS
  tmp <- ICS
  tmp[,tot:=sum(icsp),by=.(arm,age)]
  Fu5 <-  tmp[age=='0-4',.(nm='Fu5',mid=mean(tot),sd=sd(tot))]
  Fu5[,c('a','b'):=getAB(mid,sd^2)]

  ## rho,lambda
  prps <- CDD[nm %in% c('F_ICS', #pphc
                        'F_presume',                    #presume/screened = CoB
                        'F_SOC_BoA_modelled','F_alpha', #alpha
                        'F_refs',# rho = refs NOTE graph mislabelled
                        'F_refu', #LTFU
                        'DoC_0-4',
                        'DoC_5-14'
                       )]
  prps[,c('a','b'):=getAB(mid,(hi-lo)^2/3.92^2)]

  ## omega as gamma
  pstv <- CDD[nm=='F_omega_flat']
  pstv[,gmsc:=(hi-lo)^2/3.92^2 / mid]
  pstv[,gmk:=mid/gmsc]

  ## return
  list(Fu5=Fu5,prps=prps,pstv=pstv)
}


## doing as PSA parameters dependent on TB
computeCascadePSAPrev <- function(Dx   #Dx solution
                                  ){
  ## C -> D (diagnosis)
  TBP <- Dx[arm=='iph'] #NOTE could have flag here to use idh instead TODO
  TBP[,piC:=tbprev]
  TBP[,piB:=piC * CoB]
  TBP[,piA:=piB]

  ## ICS given TB
  TBP[,pl:=ifelse(location=='PHC',F_ICS,1-F_ICS)]
  TBI <- TBP[,.(tbi=sum(pl*piA)),by=.(arm,age,id)]
  ## odds for DH if TB
  TBIA <- TBP[arm=='iph',.(id,location,age,piA)]
  TBIA <- dcast(TBIA,id + age ~ location,value.var = 'piA')
  TBIA[,OR:= (DH/PHC) / ((1-DH)/(1-PHC))]
  ## TBIA[,summary(OR)]
  TBIAW <- dcast(TBIA[,.(id,age,OR)],id~age,value.var = 'OR')
  TBIW <- dcast(TBI,id~age,value.var = 'tbi')

  ## join for output
  TBics <- merge(TBIW[,.(id,tbu5=`0-4`,tbo5=`5-14`)],
                 TBIAW[,.(id,ORu5=`0-4`,ORo5=`5-14`)],by='id')

  ## return
  list(TBics=TBics,TBP=TBP)
}

## TB-depdendent cascade specificity of presuming
computeCascadePSASpec <- function(D,   #PSA data to alter
                                  TBP  #from last function
                              ){
  ## B -> C (screening)
  ## data inputs: B,C x level
  ## other input: TB prev
  ## parameter inputs: referral=rho, LTFU=lambda
  ## parameter outputs: screening spec = sigma x level
  ## rho= 1-d.idh.rltfu,0
  LTFU <- D[,.(id=1:nrow(D),d.idh.rltfu,F_omega_flat,CoB=F_presume)]
  armdep <- as.data.table(expand.grid(arm=c('idh','iph','soc'),location=c('DH','PHC'),age=c('0-4','5-14')))
  armdep <- armdep[!(location=='PHC' & arm=='idh')]
  ns <- nrow(armdep)
  armdep <- armdep[rep(1:ns,each=nrow(LTFU))]
  armdep[,id:=rep(1:nrow(LTFU),ns)]
  armdep <- merge(armdep,LTFU,by='id',all.x = TRUE)
  armdep[,lambda:=ifelse(arm=='idh',d.idh.rltfu,1e-2)]    #lambda = LTFU
  armdep[,rho:=ifelse(arm=='idh',1.0-1e-2,1e-2)]          #rho = referral TODO check
  armdep[,omega:=F_omega_flat]

  ## add in prevs
  armdep <- merge(armdep,
                  TBP[,.(location,age,id,piA,piB,piC)],
                  by=c('id','location','age'),all.x=TRUE)

  ## specificity of presuming
  Sp <- armdep[location=='PHC',.(arm,id,age,piA,piB,piC,CoB,rho,lambda,omega)]
  Sd <- armdep[location=='DH',.(arm,id,age,piA,piB,piC,CoB,rho,lambda,omega)]
  Sp[,sigmap:=1.0-CoB*(1-piC)/((1-piB)*(1-rho))]
  ## Sp[sigmap<0]                                   #OK
  Sd <- merge(Sd,Sp[,.(id,age,arm,piBp=piB,sigmap)],by=c('id','age','arm'))
  Sd[,sigmad:=1.0-CoB*(1-piC)/(1-piB) + ((1-piBp)/(1-piB))*(1-sigmap)*rho*(1-lambda)/omega]
  ## Sd[,.(id,age,arm,sigmap,sigmad)] #check

  ## return
  dcast(Sd[,.(id,age,arm,sigmap,sigmad)],id~arm+age,value.var = c('sigmap','sigmad'))
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

  keep <- c('costperATT.soc','costperATT.iph','costperATT.idh',
            'DcostperATT.iph','DcostperATT.idh',
            'Ddeaths.iph','Ddeaths.idh',
            'DLYL.iph','DLYL.idh',
            'DLYL0.iph','DLYL0.idh',
            'Dcost.iph','Dcost.idh',
            'attPC.iph','attPC.idh',
            'DcostperLYS0.iph','DcostperLYS0.idh',
            'DcostperLYS.iph','DcostperLYS.idh',
            'Dcostperdeaths.iph','Dcostperdeaths.idh',
            ## D/D
            'DcostperDATT.iph','DcostperDATT.idh',
            'DcostperDLYS0.iph','DcostperDLYS0.idh',
            'DcostperDLYS.iph','DcostperDLYS.idh',
            'DcostperDdeaths.iph','DcostperDdeaths.idh')
  scr <- c(psoc.sc,piph.sc,pidh.sc)
  scrm <- paste0(scr,'.mid')
  keep <- c(keep,scr)

  ## mid/lo/hi
  outa <- MLH(out[,..keep])

  ## more bespoke statistics
  outi <- out[,.(ICER.iph= -mean(Dcost.iph) / mean(DLYL.iph),
                 ICER.idh= -mean(Dcost.idh) / mean(DLYL.idh))]

  ## join
  outs <- do.call(cbind,list(outa$M,outa$L,outa$H,outi)) #combine

  ## pretty version
  pouts <- outs[,.(costperATT.soc = brkt(costperATT.soc.mid,costperATT.soc.lo,costperATT.soc.hi),
                   costperATT.iph = brkt(costperATT.iph.mid,costperATT.iph.lo,costperATT.iph.hi),
                   costperATT.idh = brkt(costperATT.idh.mid,costperATT.idh.lo,costperATT.idh.hi),
                   DcostperATT.iph = brkt(DcostperATT.iph.mid,DcostperATT.iph.lo,DcostperATT.iph.hi),
                   DcostperATT.idh = brkt(DcostperATT.idh.mid,DcostperATT.idh.lo,DcostperATT.idh.hi),
                   DcostperLYS0.iph = brkt(DcostperLYS0.iph.mid,
                                           DcostperLYS0.iph.lo,DcostperLYS0.iph.hi),
                   DcostperLYS0.idh = brkt(DcostperLYS0.idh.mid,DcostperLYS0.idh.lo,DcostperLYS0.idh.hi),
                   DcostperLYS.iph = brkt(DcostperLYS.iph.mid,DcostperLYS.iph.lo,DcostperLYS.iph.hi),
                   DcostperLYS.idh = brkt(DcostperLYS.idh.mid,DcostperLYS.idh.lo,DcostperLYS.idh.hi),
                   Dcostperdeaths.iph = brkt(Dcostperdeaths.iph.mid,
                                            Dcostperdeaths.iph.lo,Dcostperdeaths.iph.hi),
                   Dcostperdeaths.idh = brkt(Dcostperdeaths.idh.mid,
                                             Dcostperdeaths.idh.lo,Dcostperdeaths.idh.hi),
                   ## D/D
                   DcostperDATT.iph = brkt(DcostperDATT.iph.mid,DcostperDATT.iph.lo,DcostperDATT.iph.hi),
                   DcostperDATT.idh = brkt(DcostperDATT.idh.mid,DcostperDATT.idh.lo,DcostperDATT.idh.hi),
                   DcostperDLYS0.iph = brkt(DcostperDLYS0.iph.mid,
                                           DcostperDLYS0.iph.lo,DcostperDLYS0.iph.hi),
                   DcostperDLYS0.idh = brkt(DcostperDLYS0.idh.mid,DcostperDLYS0.idh.lo,DcostperDLYS0.idh.hi),
                   DcostperDLYS.iph = brkt(DcostperDLYS.iph.mid,DcostperDLYS.iph.lo,DcostperDLYS.iph.hi),
                   DcostperDLYS.idh = brkt(DcostperDLYS.idh.mid,DcostperDLYS.idh.lo,DcostperDLYS.idh.hi),
                   DcostperDdeaths.iph = brkt(DcostperDdeaths.iph.mid,
                                             DcostperDdeaths.iph.lo,DcostperDdeaths.iph.hi),
                   DcostperDdeaths.idh = brkt(DcostperDdeaths.idh.mid,
                                             DcostperDdeaths.idh.lo,DcostperDdeaths.idh.hi),
                   ## end D/D
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
                   attPC.iph = brkt(attPC.iph.mid,attPC.iph.lo,attPC.iph.hi),
                   attPC.idh = brkt(attPC.idh.mid,attPC.idh.lo,attPC.idh.hi),
                   DattPC.iph = brkt(attPC.iph.mid-1e2,attPC.iph.lo-1e2,attPC.iph.hi-1e2),
                   DattPC.idh = brkt(attPC.idh.mid-1e2,attPC.idh.lo-1e2,attPC.idh.hi-1e2),
                   ICER.iph=round(ICER.iph,0),ICER.idh=round(ICER.idh,0))]

  ## staged costs
  scouts <- outa$M[,..scrm]

  ## return value
  list(outs=outs,pouts=pouts,scouts=scouts)
}

## make table 1
make.table1 <- function(flout,undiscounted=FALSE){
  ## SOC
  ## DH-focussed
  ## PHC-focussed
  ## X
  ## ATT, deaths, cost, DALYs
  lflout <- copy(flout)
  if(undiscounted){
    lflout[,c("LYL.soc","LYL.idh","LYL.iph"):=.(
              LYL0.soc,LYL0.idh,LYL0.iph)]
    lflout[,c("DLYL.idh","DLYL.iph"):=.(
              DLYL0.idh,DLYL0.iph)]
  }

  keep <- c(
    ## 'full'
    "att.soc","att.idh","att.iph",
    "deaths.soc","deaths.idh","deaths.iph",
    "cost.soc","cost.idh","cost.iph",
    "LYL.soc","LYL.idh","LYL.iph",
    ## incremental
    "Datt.idh","Datt.iph",
    "Ddeaths.idh","Ddeaths.iph",
    "Dcost.idh","Dcost.iph",
    "DLYL.idh","DLYL.iph"
  )

  keepl <- c('iso3',keep)
  tmp <- lflout[,..keepl]

  ## bounds
  tmpl <- list()
  for(cn in lflout[,unique(iso3)]){
    tmpr <- tmp[iso3==cn,..keep]
    tt <- MLH(tmpr)
    tt <- cbind(tt$M,tt$L,tt$H)
    nmz <- names(tt)
    tt[,(nmz):=lapply(.SD,function(x) x * 1e5),.SDcols=nmz]
    tt[,iso3:=cn]
    tt[,ICER.idh:=tmpr[,-mean(Dcost.idh)/mean(DLYL.idh)]]
    tt[,ICER.iph:=tmpr[,-mean(Dcost.iph)/mean(DLYL.iph)]]
    tmpl[[cn]] <- tt
  }
  tmpl <- rbindlist(tmpl)

  ## table
  tmp <- tmpl[,.(iso3,
                 ## ATT
                 att.soc = brkt2(att.soc.mid,att.soc.lo,att.soc.hi),
                 att.idh = brkt2(att.idh.mid,att.idh.lo,att.idh.hi),
                 att.iph = brkt2(att.iph.mid,att.iph.lo,att.iph.hi),
                 ## deaths
                 deaths.soc = brkt2(deaths.soc.mid,deaths.soc.lo,deaths.soc.hi),
                 deaths.idh = brkt2(deaths.idh.mid,deaths.idh.lo,deaths.idh.hi),
                 deaths.iph = brkt2(deaths.iph.mid,deaths.iph.lo,deaths.iph.hi),
                 ## costs
                 cost.soc = brkt2(cost.soc.mid,cost.soc.lo,cost.soc.hi),
                 cost.idh = brkt2(cost.idh.mid,cost.idh.lo,cost.idh.hi),
                 cost.iph = brkt2(cost.iph.mid,cost.iph.lo,cost.iph.hi),
                 ## LYL
                 LYL.soc = brkt2(LYL.soc.mid,LYL.soc.lo,LYL.soc.hi),
                 LYL.idh = brkt2(LYL.idh.mid,LYL.idh.lo,LYL.idh.hi),
                 LYL.iph = brkt2(LYL.iph.mid,LYL.iph.lo,LYL.iph.hi),
                 ## ---- incremental
                 ## ATT
                 Datt.idh = brkt2(Datt.idh.mid,Datt.idh.lo,Datt.idh.hi),
                 Datt.iph = brkt2(Datt.iph.mid,Datt.iph.lo,Datt.iph.hi),
                 ## deaths
                 Ddeaths.idh = brkt2(-Ddeaths.idh.mid,-Ddeaths.idh.hi,-Ddeaths.idh.lo),
                 Ddeaths.iph = brkt2(-Ddeaths.iph.mid,-Ddeaths.iph.hi,-Ddeaths.iph.lo),
                 ## costs
                 Dcost.idh = brkt2(Dcost.idh.mid,Dcost.idh.lo,Dcost.idh.hi),
                 Dcost.iph = brkt2(Dcost.iph.mid,Dcost.iph.lo,Dcost.iph.hi),
                 ## LYL
                 DLYL.idh = brkt2(-DLYL.idh.mid,-DLYL.idh.hi,-DLYL.idh.lo),
                 DLYL.iph = brkt2(-DLYL.iph.mid,-DLYL.iph.hi,-DLYL.iph.lo),
                 ## ICERs
                 ICER.idh = format(round(ICER.idh),big.mark = ','),
                 ICER.iph = format(round(ICER.iph),big.mark = ',')
                 )]

  ## reorder, transpose
  tmp2 <- melt(tmp,id='iso3')
  tmp2 <- dcast(tmp2,variable ~ iso3,value.var = 'value')
  setkey(tmp2,variable)

  tmp2[c(c(
    ## SOC
    "att.soc","deaths.soc","cost.soc","LYL.soc",
    ## IDH
    "att.idh","deaths.idh","cost.idh","LYL.idh",
    "Datt.idh","Ddeaths.idh","Dcost.idh","DLYL.idh",
    "ICER.idh",
    ## IPH
    "att.iph","deaths.iph","cost.iph","LYL.iph",
    "Datt.iph","Ddeaths.iph","Dcost.iph","DLYL.iph",
    "ICER.iph"
    ))]

}


## ---- utilities for making CEACs
make.ceac <- function(CEA,lamz){
    crv <- lamz
    for(i in 1:length(crv)) crv[i] <- CEA[,mean(lamz[i]*Q-P>0)]
    crv
}



