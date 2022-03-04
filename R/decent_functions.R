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
## TODO - remove excess RNG here
CFRtxY <- function(age,hiv=0,art=0){#NB optimized for clarity not speed
  if(length(age)>1 & length(hiv)==1) hiv <- rep(hiv,length(age))
  if(length(age)>1 & length(art)==1) art <- rep(art,length(age))
  tmp <- P$ontxY$r(length(age))
  tmp[age=='5-14'] <- P$ontxO$r(sum(age=='5-14'))  #NB this could be achieved in  the tree model
  ## hivartOR
  Z <- P$hivartOR$r(length(age))
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
  tmp <- P$notxY$r(length(age))          #default a<5 and hiv=art=0
  tmp[age!='5-14' & hiv>0 & art==0] <- P$notxHY$r(sum(age!='5-14' & hiv>0 & art==0)) #u5,HIV+,ART-
  tmp[age!='5-14' & hiv>0 & art>0] <- P$notxHAY$r(sum(age!='5-14' & hiv>0 & art>0)) #u5,HIV+,ART+
  tmp[age=='5-14'] <- P$notxO$r(sum(age=='5-14'))    #o5, HIV-ve
  tmp[age=='5-14' & hiv>0 & art==0] <- P$notxHO$r(sum(age=='5-14' & hiv>0 & art==0)) #o5,HIV+,ART-
  tmp[age=='5-14' & hiv>0 & art>0] <- P$notxHAO$r(sum(age=='5-14' & hiv>0 & art>0)) #o5,HIV+,ART+
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
  D[tb!="noTB",d.cfr.notx:=CFRtxN(age,hiv,art)]
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

makeAttributes <- function(D){
    nrep <- nrow(D)
    D[,id:=1:nrep]
    fx <- list(age=agelevels,hiv=hivlevels,art=artlevels,tb=tblevels)
    cofx <- expand.grid(fx)
    cat('Attribute combinations used:\n')
    print(cofx)
    D <- D[rep(1:nrow(D),each=nrow(cofx))] #expand out data
    D[,names(cofx):=cofx[rep(1:nrow(cofx),nrep),]]
    ## age
    D[,value:=ifelse(age=='5-14',1-d.F.u5,d.F.u5)] #NOTE value first set here
    ## HIV/ART
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
    ## TB
    ## calculate true TB prev among initial care seeking as:
    ## tbi = f x tbp + (1-f) x tbd
    ## where: f=fraction initially seeking care at PHC; tbp=prev @ phc; tbd=prev @ dh
    ## NOTE the 'underlying' TB prev in care-seekers in controlled by soc parms
    D[,tbi:= d.soc.pphc * d.phc.tbinprsmptv + (1-d.soc.pphc) * d.dh.tbinprsmptv]
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
