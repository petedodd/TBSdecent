## IDH, at DH: TB dx bac+ on Xpert Ultra on NPA & stool/sputum, in people identified as having presumptive TB
d.idh.dh.ptbxns:=ifelse(tb=="TB+",d.idh.dh.ptbxns.se,1-d.idh.dh.ptbxns.sp)

## IPH, at DH: TB dx bac+ on Xpert Ultra on NPA & stool/sputum, in people identified as having presumptive TB                                
d.iph.dh.ptbxns:=ifelse(,d.iph.dh.ptbxns.se,1-d.iph.dh.ptbxns.sp)

## IPH, at PHC: TB dx bac+ on Xpert Ultra on NPA & stool/sputum, in people identified as having presumptive TB                               
d.iph.phc.ptbxns:=ifelse(,d.iph.phc.ptbxns.se,1-d.iph.phc.ptbxns.sp)

## IPD, at DH: TB dx bac+ on Xpert Ultra on NPA & stool/sputum, in people identified as having presumptive TB                                
d.ipd.dh.ptbxns:=ifelse(,d.ipd.dh.ptbxns.se,1-d.ipd.dh.ptbxns.sp)

## IPD, at DH: TB dx bac+ on Xpert Ultra on NPA & stool/sputum, in people chosen to be referred from PHC who were bac- on sputum test at PHC 
d.ipd.dhreftest.ptbxns:=ifelse(,d.ipd.dhreftest.ptbxns.se,1-d.ipd.dhreftest.ptbxns.sp)

## IPD, at DH: TB dx bac+ on Xpert Ultra on NPA & stool/sputum, in people chosen to be referred from PHC who were not tested at PHC          
d.ipd.dhrefnotest.ptbxns:=ifelse(,d.ipd.dhrefnotest.ptbxns.se,1-d.ipd.dhrefnotest.ptbxns.sp)

## IPH, at PHC: TB dx bac+ on Xpert Ultra on sputum, in people identified as having presumptive TB                                           
d.ipd.phc.ptbxsp:=ifelse(,d.ipd.phc.ptbxsp.se,1-d.ipd.phc.ptbxsp.sp)

## SOC, at DH: TB dx bac+ on Xpert Ultra on sputum, in people identified as having presumptive TB                                            
d.soc.dh.ptbxsp:=ifelse(,d.soc.dh.ptbxsp.se,1-d.soc.dh.ptbxsp.sp)

# SOC, at PHC: TB dx bac+ on Xpert Ultra on sputum, in people identified as having presumptive TB                                           
d.soc.phc.ptbxsp:=ifelse(,d.soc.phc.ptbxsp.se,1-d.soc.phc.ptbxsp.sp)

 ## SOC, at DH: TB dx bac+ on Xpert Ultra on sputum, in people chosen to be referred from PHC who were bac- on sputum test at PHC             
d.soc.dhreftest.ptbxsp:=ifelse(,d.soc.dhreftest.ptbxsp.se,1-d.soc.dhreftest.ptbxsp.sp)

## SOC, at DH: TB dx bac+ on Xpert Ultra on sputum, in people chosen to be referred from PHC who were not tested at PHC                      
d.soc.dhrefnotest.ptbxsp:=ifelse(,d.soc.dhrefnotest.ptbxsp.se,1-d.soc.dhrefnotest.ptbxsp.sp)

## SOC, at DH: TB dx bac+ on Xpert Ultra on GA, in people identified as having presumptive TB                                                
d.soc.dh.ptbxga:=ifelse(,d.soc.dh.ptbxga.se,1-d.soc.dh.ptbxga.sp)

## SOC, at DH: TB dx bac+ on Xpert Ultra on GA, in people chosen to be referred from PHC who were bac- on sputum test at PHC                 
d.soc.dhreftest.ptbxga:=ifelse(,d.soc.dhreftest.ptbxga.se,1-d.soc.dhreftest.ptbxga.sp)

 ## SOC, at DH: TB dx bac+ on Xpert Ultra on GA, in people chosen to be referred from PHC who were not tested at PHC                          
d.soc.dhrefnotest.ptbxga:=ifelse(,d.soc.dhrefnotest.ptbxga.se,1-d.soc.dhrefnotest.ptbxga.sp)

## IDH, at DH: TB dx clinical, in bac- people identified as having presumptive TB                                                            
d.idh.dh.test.ptbc:=ifelse(,d.idh.dh.test.ptbc.se,1-d.idh.dh.test.ptbc.sp)

## IPH, at DH: TB dx clinical, in bac- people identified as having presumptive TB                                                            
d.iph.dh.test.ptbc:=ifelse(,d.iph.dh.test.ptbc.se,1-d.iph.dh.test.ptbc.sp)

## IPH, at PHC: TB dx clinical, in bac- people identified as having presumptive TB                                                           
d.iph.phc.test.ptbc:=ifelse(,d.iph.phc.test.ptbc.se,1-d.iph.phc.test.ptbc.sp)

 ## IPD, at DH: TB dx clinical, in bac- people identified as having presumptive TB                                                            
d.ipd.dh.test.ptbc:=ifelse(,d.ipd.dh.test.ptbc.se,1-d.ipd.dh.test.ptbc.sp)

## IPD, at PHC: TB dx clinical, in bac- people identified as having presumptive TB                                                           
d.ipd.phc.test.ptbc:=ifelse(,.se,1-.sp)

## SOC, at DH: TB dx clinical, in bac- people identified as having presumptive TB                                                            
d.soc.dh.test.ptbc:=ifelse(,d.soc.dh.test.ptbc.se,1-d.soc.dh.test.ptbc.sp)

## SOC, at PHC: TB dx clinical, in bac- people identified as having presumptive TB                                                           
d.soc.phc.test.ptbc:=ifelse(,d.soc.phc.test.ptbc.se,1-d.soc.phc.test.ptbc.sp)

## IPD, at PHC: TB dx clinical, in untested people identified as having presumptive TB                                                       
d.ipd.phc.notest.ptbc:=ifelse(,d.ipd.phc.notest.ptbc.se,1-d.ipd.phc.notest.ptbc.sp)

## SOC, at DH: TB dx clinical, in untested people identified as having presumptive TB                                                        
d.soc.dh.notest.ptbc:=ifelse(,d.soc.dh.notest.ptbc.se,1-d.soc.dh.notest.ptbc.sp)

## SOC, at PHC: TB dx clinical, in untested people identified as having presumptive TB                                                       
d.soc.phc.notest.ptbc:=ifelse(,d.soc.phc.notest.ptbc.se,1-d.soc.phc.notest.ptbc.sp)

## IPD, at DH: TB dx clinical, in bac- people chosen to be referred from PHC who were bac- on sputum test at PHC                             
d.ipd.dhreftest.test.ptbc:=ifelse(,d.ipd.dhreftest.test.ptbc.se,1-d.ipd.dhreftest.test.ptbc.sp)

d.soc.dhreftest.test.ptbc.sp     ## SOC, at DH: TB dx clinical, in bac- people chosen to be referred from PHC who were bac- on sputum test at PHC                             
(,d.soc.dhreftest.test.ptbc.se,1-d.soc.dhreftest.test.ptbc.sp)

## SOC, at DH: TB dx clinical, in untested people chosen to be referred from PHC who were bac- on sputum test at PHC                         
d.soc.dhreftest.notest.ptbc:=ifelse(,d.soc.dhreftest.notest.ptbc.se,1-d.soc.dhreftest.notest.ptbc.sp)

## IDH, at DH: TB dx clinical, in bac- people chosen to be referred from PHC who were not tested at PHC                                      
d.ipd.dhrefnotest.test.ptbc:=ifelse(,d.ipd.dhrefnotest.test.ptbc.se,1-d.ipd.dhrefnotest.test.ptbc.sp)

## SOC, at DH: TB dx clinical, in bac- people chosen to be referred from PHC who were not tested at PHC                                      
d.soc.dhrefnotest.test.ptbc:=ifelse(,d.soc.dhrefnotest.test.ptbc.se,1-d.soc.dhrefnotest.test.ptbc.sp)

## SOC, at DH: TB dx clinical, in untested people chosen to be referred from PHC who were not tested at PHC either                           
d.soc.dhrefnotest.notest.ptbc:=ifelse(,d.soc.dhrefnotest.notest.ptbc.se,1-d.soc.dhrefnotest.notest.ptbc.sp)

## IDH, at DH: TB dx clinical at 7-day reassessment, in people bac- clin- in previous DH testing                                             
d.idh.dh.test7.ptbc:=ifelse(,d.idh.dh.test7.ptbc.se,1-d.idh.dh.test7.ptbc.sp)

## IPH, at DH: TB dx clinical at 7-day reassessment, in people bac- clin- in previous DH testing                                             
d.iph.dh.test7.ptbc:=ifelse(,d.iph.dh.test7.ptbc.se,1-d.iph.dh.test7.ptbc.sp)

## IPD, at DH: TB dx clinical at 7-day reassessment, in people bac- clin- in previous DH testing                                             
d.ipd.dh.test7.ptbc:=ifelse(,d.ipd.dh.test7.ptbc.se,1-d.ipd.dh.test7.ptbc.sp)

## IPH, at DH: TB dx clinical following CXR, referred to DH due to ongoing symptoms at 7-day reassessment following bac- clin- test at PHC   
d.iph.dhreftest7.ptbcxr:=ifelse(,d.iph.dhreftest7.ptbcxr.se,1-d.iph.dhreftest7.ptbcxr.sp)
