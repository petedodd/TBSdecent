## TODO questions:
## - NPA and Stool poss
## - check this is or
## - costing for this
## - stool/sputum
## - flag assumption = groups for SA or pending data

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
  d.idh.dh.ptbxns:=TBbac4ST1orST2(ST.poss,ifelse(tb=='TB+',sens.xnpa,1-spec.xnpa),
                                  NPA.poss,ifelse(tb=='TB+',sens.xstool,1-spec.xstool))
  ## IPH, at DH: TB dx bac+ on Xpert Ultra on NPA & stool/sputum, in people identified as having presumptive TB
  d.iph.dh.ptbxns:=TBbac4ST1orST2(ST.poss,ifelse(tb=='TB+',sens.xnpa,1-spec.xnpa),
                                  NPA.poss,ifelse(tb=='TB+',sens.xstool,1-spec.xstool))
  ## IPH, at PHC: TB dx bac+ on Xpert Ultra on NPA & stool/sputum, in people identified as having presumptive TB
  d.iph.phc.ptbxns:=TBbac4ST1orST2(ST.poss,ifelse(tb=='TB+',sens.xnpa,1-spec.xnpa),
                                   NPA.poss,ifelse(tb=='TB+',sens.xstool,1-spec.xstool))
  ## IPD, at DH: TB dx bac+ on Xpert Ultra on NPA & stool/sputum, in people identified as having presumptive TB
  d.ipd.dh.ptbxns:=TBbac4ST1orST2(ST.poss,ifelse(tb=='TB+',sens.xnpa,1-spec.xnpa),
                                  NPA.poss,ifelse(tb=='TB+',sens.xstool,1-spec.xstool))
  ## IPD, at DH: TB dx bac+ on Xpert Ultra on NPA & stool/sputum, in people chosen to be referred from PHC who were bac- on sputum test at PHC
  d.ipd.dhreftest.ptbxns:=TBbac4ST1orST2(ST.poss,ifelse(tb=='TB+',sens.xnpa,1-spec.xnpa),
                                         NPA.poss,ifelse(tb=='TB+',sens.xstool,1-spec.xstool))
  ## IPD, at DH: TB dx bac+ on Xpert Ultra on NPA & stool/sputum, in people chosen to be referred from PHC who were not tested at PHC
  d.ipd.dhrefnotest.ptbxns:=TBbac4ST1orST2(ST.poss,ifelse(tb=='TB+',sens.xnpa,1-spec.xnpa),
                                           NPA.poss,ifelse(tb=='TB+',sens.xstool,1-spec.xstool))

  ## ------- X on sputum ------- #NOTE hardcodes 0 sample availability for u5
  ## IPH, at PHC: TB dx bac+ on Xpert Ultra on sputum, in people identified as having presumptive TB
  d.ipd.phc.ptbxsp:=TBbacsampletest(ifelse(age=='5-14',ES.poss.phc,0),
                                    ifelse(tb=='TB+',sens.xsputum,1-spec.xsputum))
  ## SOC, at DH: TB dx bac+ on Xpert Ultra on sputum, in people identified as having presumptive TB
  d.soc.dh.ptbxsp:=TBbacsampletest(ifelse(age=='5-14',ES.poss.dh,0),
                                   ifelse(tb=='TB+',sens.xsputum,1-spec.xsputum))
  ## SOC, at PHC: TB dx bac+ on Xpert Ultra on sputum, in people identified as having presumptive TB
  d.soc.phc.ptbxsp:=TBbacsampletest(ifelse(age=='5-14',ES.poss.phc,0),
                                    ifelse(tb=='TB+',sens.xsputum,1-spec.xsputum))
  ## SOC, at DH: TB dx bac+ on Xpert Ultra on sputum, in people chosen to be referred from PHC who were bac- on sputum test at PHC
  d.soc.dhreftest.ptbxsp:=TBbacsampletest(ifelse(age=='5-14',ES.poss.dh,0),
                                          ifelse(tb=='TB+',sens.xsputum,1-spec.xsputum))
  ## SOC, at DH: TB dx bac+ on Xpert Ultra on sputum, in people chosen to be referred from PHC who were not tested at PHC
  d.soc.dhrefnotest.ptbxsp:=TBbacsampletest(ifelse(age=='5-14',ES.poss.dh,0),
                                            ifelse(tb=='TB+',sens.xsputum,1-spec.xsputum))

  ## ------- X on GA -------
  ## SOC, at DH: TB dx bac+ on Xpert Ultra on GA, in people identified as having presumptive TB
  d.soc.dh.ptbxga:=TBbacsampletest(ifelse(age=='5-14',GA.poss.dh.o5,GA.poss.dh.u5),
                                   ifelse(tb=='TB+',sens.xga,1-spec.xga))
  ## SOC, at DH: TB dx bac+ on Xpert Ultra on GA, in people chosen to be referred from PHC who were bac- on sputum test at PHC
  d.soc.dhreftest.ptbxga:=TBbacsampletest(ifelse(age=='5-14',GA.poss.dh.o5,GA.poss.dh.u5),
                                          ifelse(tb=='TB+',sens.xga,1-spec.xga))
  ## SOC, at DH: TB dx bac+ on Xpert Ultra on GA, in people chosen to be referred from PHC who were not tested at PHC
  d.soc.dhrefnotest.ptbxga:=TBbacsampletest(ifelse(age=='5-14',GA.poss.dh.o5,GA.poss.dh.u5),
                                            ifelse(tb=='TB+',sens.xga,1-spec.xga))

  ## ------- clinical -------
  ## IDH, at DH: TB dx clinical, in bac- people identified as having presumptive TB
  d.idh.dh.test.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)
  ## IPH, at DH: TB dx clinical, in bac- people identified as having presumptive TB
  d.iph.dh.test.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)
  ## IPH, at PHC: TB dx clinical, in bac- people identified as having presumptive TB
  d.iph.phc.test.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)
  ## IPD, at DH: TB dx clinical, in bac- people identified as having presumptive TB
  d.ipd.dh.test.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)
  ## IPD, at PHC: TB dx clinical, in bac- people identified as having presumptive TB
  d.ipd.phc.test.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)
  ## SOC, at DH: TB dx clinical, in bac- people identified as having presumptive TB
  d.soc.dh.test.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)
  ## SOC, at PHC: TB dx clinical, in bac- people identified as having presumptive TB
  d.soc.phc.test.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)
  ## IPD, at PHC: TB dx clinical, in untested people identified as having presumptive TB
  d.ipd.phc.notest.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)
  ## SOC, at DH: TB dx clinical, in untested people identified as having presumptive TB
  d.soc.dh.notest.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)
  ## SOC, at PHC: TB dx clinical, in untested people identified as having presumptive TB
  d.soc.phc.notest.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)
  ## IPD, at DH: TB dx clinical, in bac- people chosen to be referred from PHC who were bac- on sputum test at PHC
  d.ipd.dhreftest.test.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)
  ## SOC, at DH: TB dx clinical, in bac- people chosen to be referred from PHC who were bac- on sputum test at PHC
  d.soc.dhreftest.test.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)
  ## SOC, at DH: TB dx clinical, in untested people chosen to be referred from PHC who were bac- on sputum test at PHC
  d.soc.dhreftest.notest.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)
  ## IDH, at DH: TB dx clinical, in bac- people chosen to be referred from PHC who were not tested at PHC
  d.ipd.dhrefnotest.test.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)
  ## SOC, at DH: TB dx clinical, in bac- people chosen to be referred from PHC who were not tested at PHC
  d.soc.dhrefnotest.test.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)
  ## SOC, at DH: TB dx clinical, in untested people chosen to be referred from PHC who were not tested at PHC either
  d.soc.dhrefnotest.notest.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)
  ## IDH, at DH: TB dx clinical at 7-day reassessment, in people bac- clin- in previous DH testing
  d.idh.dh.test7.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)
  ## IPH, at DH: TB dx clinical at 7-day reassessment, in people bac- clin- in previous DH testing
  d.iph.dh.test7.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)
  ## IPD, at DH: TB dx clinical at 7-day reassessment, in people bac- clin- in previous DH testing
  d.ipd.dh.test7.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)
  ## IPH, at DH: TB dx clinical following CXR, referred to DH due to ongoing symptoms at 7-day reassessment following bac- clin- test at PHC
  d.iph.dhreftest7.ptbcxr:=ifelse(tb!='noTB',sens.clin,1-spec.clin)

}
