## looking at the decentralization cascade data
library(here)
library(data.table)
library(ggplot2)
library(scales)
library(brms)
library(glue)

## helpers
rot45 <- theme(axis.text.x = element_text(angle = 45, hjust = 1))
## add CIs
## helper functions
getCI1 <- function(x) binom.test(x[1],x[2],p=.025)$conf.int #k,N
getCI <- function(k,N) t(apply(cbind(k,N),1,getCI1))
getCI(c(5,5,5),c(10,10,10))
## function to add binomial CIs
MLH <- function(k,N) {
  if(length(k)!=length(N)) stop('k and N have different lengths!')
  k <- as.integer(k); N <- as.integer(N)
  mid <- lo <- hi <- rep(NA,length(k))
  who <- which(!is.na(k) & !is.na(N) & N>0)
  if(any(k[who]<0)){ stop('k<0!')}
  if(any(k[who]>N[who])){ stop('k>N!')}
  HL <- getCI(k[who],N[who])
  mid[who] <- k[who]/N[who]
  lo[who] <- HL[,1]; hi[who] <- HL[,2]
  list(mid*1,lo*1,hi*1)
}
MLH(c(5,5,5),c(10,10,10))
odds <- function(x) x/(1-x)
logit <- function(x)log(x/(1-x))
ilogit <- function(x)1/(1+exp(-x))

## extractor function
bayessmy <- function(mdl){
  fe <- fixef(mdl)
  smy <- data.table(type='model',location='All',
                    mid=inv_logit_scaled(fe[1]),
                    lo=inv_logit_scaled(fe[3]),
                    hi=inv_logit_scaled(fe[4]))
  ss <- coef(mdl)$location[,,]
  S <- as.data.table(coef(mdl)$location[,,])
  S[,type:='model']
  S[,location:=rownames(ss)]
  S[,c('mid','lo','hi'):=.(inv_logit_scaled(Estimate),
                           inv_logit_scaled(Q2.5),
                           inv_logit_scaled(Q97.5))]
  rbind(S[,.(type,location,mid,lo,hi)],
        smy)
}


## read data
D <- fread(here('indata/cascades/MainCascade.csv')) #cascade
R <- fread(here('indata/cascades/refd.csv')) #referrals
U <- fread(here('indata/cascades/refu.csv')) #refusals NOTE could combine with above
P <- fread(here('indata/cascades/pphcdata.csv')) #initial care-seeking
key <- data.table(
  location=D[,unique(location)],
  cn=R$cn[c(1,4,7,2,3,5,6)])

## changes to D
D$OPD <- gsub(',','',D$OPD)
D$OPD <- as.integer(D$OPD)
D$Screened <- gsub(',','',D$Screened)
D$Screened <- as.integer(D$Screened)

## levels
arms <- c('SOC','IDH','IPH')
acats <- c('0-4','5-14')
stage <- c('OPD','Screened','Presumptive','Treated')
lvl <- c('PHC','DH')

## reformat
DM <- melt(D,id=c('location','Arm','Facility','Age'))
DM$Arm <- factor(DM$Arm,levels=arms,ordered = TRUE)
DM$Facility <- factor(DM$Facility,levels=lvl,ordered = TRUE)
DM$Age <- factor(DM$Age,levels=acats,ordered = TRUE)
DM$variable <- factor(DM$variable,levels=stage,ordered = TRUE)
names(DM)[names(DM)=='variable'] <- 'stage'
DM

## aggregate over facility
DMA <- DM[,.(value=sum(value)),by=.(location,Arm,Age,stage)]

## wide again
DMW <- dcast(DM,location+Arm+Facility+Age~stage,value.var = 'value')
DMW[,c('Screened/OPD',
       'Presumptive/Screened',
       'Treated/Presumptive'):=
       .(Screened/OPD,Presumptive/Screened,Treated/Presumptive)]
DC <- melt(DMW[,.(location,Arm,Facility,Age,
                  `Screened/OPD`,
                  `Presumptive/Screened`,
                  `Treated/Presumptive`)],
           id=c('location','Arm','Facility','Age'))

stp <- c('Screened/OPD','Presumptive/Screened','Treated/Presumptive')
names(DC)[5] <- 'step'
DC$step <- factor(DC$step,levels=stp,ordered = TRUE)


## working with version aggreagted over location
## aggregate over facility
DMA <- DM[,.(value=sum(value)),by=.(location,Arm,Age,stage)]
DMWA <- dcast(DMA,location+Arm+Age~stage,value.var = 'value')
stpl <- c('Screened/OPD','Presumptive/Screened',
          'Treated/Presumptive','Treated/OPD')
DMWA[,(stpl):=.(Screened/OPD,Presumptive/Screened,
                Treated/Presumptive,Treated/OPD)]

DA <- melt(DMWA[,.(location,Arm,Age,
                   `Screened/OPD`,`Presumptive/Screened`,
                   `Treated/Presumptive`,`Treated/OPD`)],
           id=c('location','Arm','Age'))

names(DA)[4] <- 'step'
DA$step <- factor(DA$step,levels=stpl,ordered = TRUE)

DA[Arm!='SOC' & location=='All'][order(step,Age,Arm)]
DA[Arm!='SOC' & location=='All' & Arm=='IDH'][order(step,Age,Arm)]
DA[Arm!='SOC' & location=='All' & Arm=='IPH'][order(step,Age,Arm)]

DA[step=='Treated/OPD',.(location,Arm,Age,1/value)]

## make directories if missing
if(!file.exists(here('graphs/cascades'))) dir.create(here('graphs/cascades'))
if(!file.exists(here('graphs/cascades/data'))) dir.create(here('graphs/cascades/data'))

ggplot(DA,aes(Arm,value,
              col=location,
              shape=factor(location=='All'),
              group=paste(location,Age)))+
  geom_point()+
  geom_line()+
  theme_light()+
  facet_grid(step~Age,scales = 'free_y')

ggsave(here('graphs/cascades/ArmCountryData.png'),w=7,h=7)

save(DMWA,file=here('graphs/cascades/data/DMWA.Rdata'))


## changes to R & U & P
R <- merge(R,key,by='cn')
R$location <- factor(R$location,levels=R$location,ordered = TRUE)
R[,c('mid','lo','hi'):=MLH(lost,total)]
R[,type:='data']

U <- merge(U,key,by='cn')
U$location <- factor(U$location,levels=U$location,ordered = TRUE)
U[,c('mid','lo','hi'):=MLH(no,asked)]
U[,type:='data']

P <- merge(P,key,by='cn')
P$location <- factor(P$location,levels=P$location,ordered = TRUE)
ditch <- grep('A',names(P),value = TRUE)
ditch <- c(ditch,'d.idh.pphc','d.iph.pphc')
P[,(ditch):=NULL]
P <- melt(P[,.(location,c,d.idh.pphc.k, d.idh.pphc.N, d.iph.pphc.k, d.iph.pphc.N)],
     id=c('location','c'))
P[,Arm:=ifelse(grepl('iph',variable),'IPH','IDH')]
P[,qty:=ifelse(grepl('k',variable),'k','N')]
P <- dcast(P,location + c + Arm ~ qty,value.var = 'value')
P[,c('mid','lo','hi'):=MLH(k,N)]
P[,type:='data']


## ==== referrals
## forest plots
ggplot(R,aes(location,y=mid,ymin=lo,ymax=hi))+
  scale_y_continuous(label=percent,limits=c(0,1))+
  geom_pointrange()+
  geom_pointrange(data=U,col='red')+
  coord_flip()+
  ylab('LTFU')


## modelling: don't include UGA & MOZ as different health system structures
mdl2 <- brm(formula=lost | trials(total) ~1+(1|location),
           data=R[!location %in% c('All','Uganda','Mozambique')],family=binomial(link='logit'))

## extract data
SM <- bayessmy(mdl2)

## joing and plot both
RS <- rbind(R[,.(type,location,mid,lo,hi)],SM)
ggplot(RS,aes(location,y=mid,ymin=lo,ymax=hi,group=type,col=type))+
  scale_y_continuous(label=percent,limits=c(0,1))+
  geom_pointrange(position=position_dodge(0.1),shape=1)+coord_flip()+
  ylab('Referral') + theme_classic()+ggpubr::grids()+theme(legend.position = 'top')

ggsave(here('graphs/cascades/F_refs.png'),w=5,h=6)

fwrite(RS[type=='model'],file=here('graphs/cascades/data/F_refs.csv'))


## ==== refusals

## modelling
mdl2b <- brm(formula=no | trials(asked) ~1+(1|location),
            data=U[!location %in% c('All')],family=binomial(link='logit'))

## extract data
SMb <- bayessmy(mdl2b)

## joing and plot both
RSb <- rbind(U[,.(type,location,mid,lo,hi)],SMb)
ggplot(RSb,aes(location,y=mid,ymin=lo,ymax=hi,group=type,col=type))+
  scale_y_continuous(label=percent,limits=c(0,1))+
  geom_pointrange(position=position_dodge(0.1),shape=1)+coord_flip()+
  ylab('Refusal') + theme_classic()+ggpubr::grids()+theme(legend.position = 'top')

ggsave(here('graphs/cascades/F_refu.png'),w=5,h=6)

fwrite(RSb[type=='model'],file=here('graphs/cascades/data/F_refu.csv'))



## ===== rest of cascade


## ----- Screening
## -------- B/A
## w/o age/arm
tmp <- DMWA[location!='All' & Arm!='SOC',.(Screened=sum(Screened),OPD=sum(OPD)),by=location]
## tmp <- DMWA[location!='All' & Arm!='SOC',.(location,Screened,OPD)]
mdl3 <- brm(formula=Screened | trials(OPD) ~1+(1|location),
            data=tmp,
            family=binomial(link='logit'))


## extract data
SM1 <- bayessmy(mdl3)
R1 <- DMWA[Arm!='SOC',.(type='data',location,Screened,OPD)]
R1[,c('mid','lo','hi'):=MLH(Screened,OPD)]

## joing and plot both
RS1 <- rbind(R1[,.(type,location,mid,lo,hi)],SM1)
ggplot(RS1,aes(location,y=mid,ymin=lo,ymax=hi,group=type,col=type))+
  scale_y_continuous(label=percent)+
  geom_pointrange(position=position_dodge(0.1),shape=1)+coord_flip()+
  ylab('Screening/OPD') + theme_classic()+ggpubr::grids()+theme(legend.position = 'top')

ggsave(here('graphs/cascades/F_screen.png'),w=5,h=6)

fwrite(RS1[type=='model'],file=here('graphs/cascades/data/F_screenOPD.csv'))




## ---- Presuming
## --------- C/B
## w/o age/arm
tmp <- DMWA[location!='All' & Arm!='SOC',
            .(Screened=sum(Screened),Presumptive=sum(Presumptive)),by=location]

mdl4 <- brm(formula=Presumptive | trials(Screened) ~1+(1|location),
            data=tmp,
            family=binomial(link='logit'))


## extract data
SM2 <- bayessmy(mdl4)
R2 <- DMWA[Arm!='SOC',.(type='data',location,Screened,Presumptive)]
R2[,c('mid','lo','hi'):=MLH(Presumptive,Screened)]

## joing and plot both
RS2 <- rbind(R2[,.(type,location,mid,lo,hi)],SM2)
ggplot(RS2,aes(location,y=mid,ymin=lo,ymax=hi,group=type,col=type))+
  scale_y_continuous(label=percent)+
  geom_pointrange(position=position_dodge(0.1),shape=1)+coord_flip()+
  ylab('Presumptive/Screening') + theme_classic()+ggpubr::grids()+theme(legend.position = 'top')

ggsave(here('graphs/cascades/F_presume.png'),w=5,h=6)

fwrite(RS2[type=='model'],file=here('graphs/cascades/data/F_presume.csv'))



## --------- initial level of care-seeking
## now w/o age/arm
tmp <- P[location!='All',.(k=sum(k),N=sum(N)),by=location]

mdl6 <- brm(formula=k | trials(N) ~1+(1|location),
            data=tmp,
            family=binomial(link='logit'))


## extract data
SM6 <- bayessmy(mdl6)

R3 <- DMWA[Arm!='SOC',.(type='data',location,Treated,Presumptive)]
R3[,c('mid','lo','hi'):=MLH(Treated,Presumptive)]

## joing and plot both
CS6 <- rbind(P[,.(type,location,mid,lo,hi)],SM6)

ggplot(CS6,aes(location,y=mid,ymin=lo,ymax=hi,group=type,col=type))+
  scale_y_continuous(label=percent)+
  geom_pointrange(position=position_dodge(0.1),shape=1)+
  coord_flip()+
  ylab(expression(A[p]/A[d])) + theme_classic()+ggpubr::grids()+theme(legend.position = 'top')

ggsave(here('graphs/cascades/F_ICS_raw.png'),w=5,h=6) #uncorrected

## correct
ch <- function(x,a) x/(x+a)
work <- c('mid','lo','hi')
old <- paste0(work,'_raw')
CS6[,(old):=.(mid,lo,hi)]
CS6 <- merge(CS6,unique(P[,.(c,location)]),by='location')
CS6[,c(work):=.(ch(mid_raw,c),ch(lo_raw,c),ch(hi_raw,c))]



ggplot(CS6,aes(location,y=mid,ymin=lo,ymax=hi,group=type,col=type))+
  scale_y_continuous(label=percent)+
  geom_pointrange(position=position_dodge(0.1),shape=1)+
  coord_flip()+
  ylab('Proportion initially seeking care at PHC') +
  theme_classic()+ggpubr::grids()+theme(legend.position = 'top')

ggsave(here('graphs/cascades/F_ICS_cor.png'),w=5,h=6)


fwrite(CS6[type=='model',.(location,mid,lo,hi)],
       file=here('graphs/cascades/data/F_ICS.csv'))



## --- omega

## explore
RO <- D[,.(location,Arm,Facility,Age,Screened)]
ROW <- dcast(RO,location+Arm+Age~Facility,value.var = 'Screened')

ggplot(ROW[Arm!='SOC'],aes(Arm,DH/PHC,shape=Age,lty=Age,col=location,group=paste(location,Age)))+
  geom_point()+geom_line()+
  theme_classic()+ggpubr::grids()

## aggregate by age, keep arm consistent as intervention changes screening
ROW <- ROW[Arm!='SOC',.(DH=sum(DH),PHC=sum(PHC)),by=.(location,Arm)]
ROW[,total:=DH+PHC]
ROW[,type:='data']
ROW[,c('mid','lo','hi'):=MLH(DH,total)]

## -- flat version (ie no ARM dependence)
ROWf <- ROW[Arm!='SOC',.(DH=sum(DH),PHC=sum(PHC)),by=.(location)]
ROWf[,total:=DH+PHC]
ROWf[,type:='data']
ROWf[,c('mid','lo','hi'):=MLH(DH,total)]

mdl7 <- brm(formula=DH | trials(total) ~1+(1|location),
               data=ROWf[location !='All'],
               family=binomial(link='logit'))

## extract data
SM7f <- bayessmy(mdl7)


## joing and plot both
CS7f <- rbind(ROWf[,.(type,location,mid,lo,hi)],SM7f)

ggplot(CS7f,aes(location,y=mid,ymin=lo,ymax=hi,group=paste(type),col=type))+
  scale_y_continuous(label=percent)+
  geom_pointrange(position=position_dodge(0.1),shape=1)+
  coord_flip()+
  ylab('Fraction of screening at DH') + theme_classic()+ggpubr::grids()+theme(legend.position = 'top')

ggsave(here('graphs/cascades/F_Omega_model_flat.png'),w=5,h=6)



## correct
work <- c('mid','lo','hi')
old <- paste0(work,'_raw')
CS7f[,(old):=.(mid,lo,hi)]
CS7f[,c(work):=.(odds(mid_raw),odds(lo_raw),odds(hi_raw))]

fwrite(CS7f[type=='model',.(location,mid,lo,hi)],
       file=here('graphs/cascades/data/F_omega_flat.csv'))


## --- alpha
## explore
RA <- D[,.(location,Arm,Facility,Age,Screened,OPD)]


ggplot(RA[Arm!='SOC'],aes(Arm,Screened/OPD,shape=Age,lty=Facility,col=location,
                          group=paste(Facility,location,Age)))+
  geom_point()+geom_line()+
  theme_classic()+ggpubr::grids()

ggsave(here('graphs/cascades/F_Alpha_ScreenCov.png'),w=7,h=5)

## aggregate by age, arm, 
RAA <- RA[Arm!='SOC',.(Screened=sum(Screened),OPD=sum(OPD)),by=.(location)]
RAA[,type:='data']
RAA[,c('mid','lo','hi'):=MLH(Screened,OPD)]

mdl8 <- brm(formula=Screened | trials(OPD) ~1+(1|location),
               data=RAA[location !='All'],
               family=binomial(link='logit'))


## extract data
SM8 <- bayessmy(mdl8)

## joing and plot both
CS8 <- rbind(RAA[,.(type,location,mid,lo,hi)],SM8)

ggplot(CS8,aes(location,y=mid,ymin=lo,ymax=hi,group=paste(type),col=type))+
  scale_y_continuous(label=percent)+
  geom_pointrange(position=position_dodge(0.1),shape=1)+
  coord_flip()+
  ylab('Screening coverage') + theme_classic()+ggpubr::grids()+theme(legend.position = 'top')

ggsave(here('graphs/cascades/F_alpha_model.png'),w=5,h=6)


fwrite(CS8[type=='model',.(location,mid,lo,hi)],
       file=here('graphs/cascades/data/F_alpha.csv'))

## ---- Treating/Presumed
## ---------- D/C
## now w/o age/arm
tmp <- DMWA[location!='All' & Arm!='SOC',
            .(Treated=sum(Treated),Presumptive=sum(Presumptive)),by=location]

mdl5 <- brm(formula=Treated | trials(Presumptive) ~1+(1|location),
            data=tmp,
            family=binomial(link='logit'))


## extract data
SM3 <- bayessmy(mdl5)
R3 <- DMWA[Arm!='SOC',.(type='data',location,Treated,Presumptive)]
R3[,c('mid','lo','hi'):=MLH(Treated,Presumptive)]

## joing and plot both
RS3 <- rbind(R3[,.(type,location,mid,lo,hi)],SM3)
ggplot(RS3,aes(location,y=mid,ymin=lo,ymax=hi,group=type,col=type))+
  scale_y_continuous(label=percent)+
  geom_pointrange(position=position_dodge(0.1),shape=1)+coord_flip()+
  ylab('Treated/Presumptive') + theme_classic()+ggpubr::grids()+theme(legend.position = 'top')

ggsave(here('graphs/cascades/F_treat.png'),w=5,h=6)

fwrite(RS3[type=='model'],file=here('graphs/cascades/data/F_treat.csv'))



## NOTE remainder are stratified by age and arm
## =================== loop
TPD <- TPM <- TPS <- list()
TOD <- TOM <- TOS <- list()
for(arm in c('SOC','IPH','IDH')){
  for(age in c('0-4','5-14')){
    cat('working on ',arm,' & ',age,'...\n')
    if(arm!='SOC'){
      ## --- treated/presumptive
      ## data
      TPD[[glue("{arm}{age}")]] <- DMWA[location!='All' & Arm==arm & Age==age,
                                        .(Treated=sum(Treated),Presumptive=sum(Presumptive)),
                                        by=location]
      ## model
      TPM[[glue("{arm}{age}")]] <- brm(formula=Treated | trials(Presumptive) ~1+(1|location),
                                       data=TPD[[glue("{arm}{age}")]],
                                       family=binomial(link='logit'))

      ## extracts
      SMx <- bayessmy(TPM[[glue("{arm}{age}")]])
      Rx <- DMWA[Arm==arm & Age==age,.(type='data',location,Treated,Presumptive)]
      Rx[,c('mid','lo','hi'):=MLH(Treated,Presumptive)]
      tmp <- rbind(Rx[,.(type,location,mid,lo,hi)],SMx)
      tmp[,c('Arm','Age'):=.(arm,age)]
      TPS[[glue("{arm}{age}")]] <- tmp
    }

    ## --- treated/OPD (all arms, including SOC)
    ## data
    TOD[[glue("{arm}{age}")]] <- DMWA[location!='All' & Arm==arm & Age==age,
                                      .(Treated=sum(Treated),OPD=sum(OPD)),
                                      by=location]
    ## model
    TOM[[glue("{arm}{age}")]] <- brm(formula=Treated | trials(OPD) ~1+(1|location),
                                     data=TOD[[glue("{arm}{age}")]],
                                     family=binomial(link='logit'))

    ## extracts
    SMx <- bayessmy(TOM[[glue("{arm}{age}")]])
    Rx <- DMWA[Arm==arm & Age==age,.(type='data',location,Treated,OPD)]
    Rx[,c('mid','lo','hi'):=MLH(Treated,OPD)]
    tmp <- rbind(Rx[,.(type,location,mid,lo,hi)],SMx)
    tmp[,c('Arm','Age'):=.(arm,age)]
    TOS[[glue("{arm}{age}")]] <- tmp

  }
}
TPS <- rbindlist(TPS)
TOS <- rbindlist(TOS)

fwrite(TPS[type=='model'],file=here('graphs/cascades/data/TPS.csv'))
fwrite(TOS[type=='model'],file=here('graphs/cascades/data/TOS.csv'))


ggplot(TPS,aes(location,y=mid,ymin=lo,ymax=hi,group=type,col=type))+
  scale_y_continuous(label=percent)+
  geom_pointrange(position=position_dodge(0.1),shape=1)+coord_flip()+
  facet_grid(Arm ~ Age)+
  ylab('Treated/Presumptive') +
  theme_light()+
  theme(legend.position = 'top')

ggsave(here('graphs/cascades/F_treat_AA.png'),w=10,h=10)


ggplot(TOS,aes(location,y=mid,ymin=lo,ymax=hi,group=type,col=type))+
  scale_y_continuous(label=percent)+
  geom_pointrange(position=position_dodge(0.1),shape=1)+coord_flip()+
  facet_grid(Arm ~ Age)+
  ylab('Treated/OPD') +
  theme_light()+
  theme(legend.position = 'top')

ggsave(here('graphs/cascades/F_treatOPD_AA.png'),w=10,h=15)



## ================= SOC screening prediction
## missing B/A for SOC
## do have D/A for SOC
## SOC:B/A =
## SOC:(D/A) / (D/B)

## explore
SP <- D[,.(location,Arm,Facility,Age,Treated,OPD)]


ggplot(SP,aes(Arm,Treated/OPD,shape=Age,lty=Facility,col=location,
                          group=paste(Facility,location,Age)))+
  geom_point()+geom_line()+
  facet_wrap(~Facility,scales='free')+
  theme_classic()+ggpubr::grids()

ggplot(SP,aes(Age,Treated/OPD,shape=Arm,lty=Facility,col=location,
              group=paste(Facility,location,Arm)))+
  geom_point()+geom_line()+
  facet_wrap(~Facility,scales='free')+
  theme_classic()+ggpubr::grids()


## aggregate by Age
RSP <- SP[Arm=='SOC',.(Treated=sum(Treated),OPD=sum(OPD)),by=.(location)]
RSP[,type:='data']
RSP[,c('mid','lo','hi'):=MLH(Treated,OPD)]

mdl9 <- brm(formula=Treated | trials(OPD) ~1+(1|location),
               data=RSP[location !='All'],
               family=binomial(link='logit'))


## extract data
SM9 <- bayessmy(mdl9)

## joing and plot both
CS9 <- rbind(RSP[,.(type,location,mid,lo,hi)],SM9)

ggplot(CS9,aes(location,y=mid,ymin=lo,ymax=hi,group=paste(type),col=type))+
  scale_y_continuous(label=percent)+
  geom_pointrange(position=position_dodge(0.1),shape=1)+
  coord_flip()+
  ylab('SOC D/A') + theme_classic()+ggpubr::grids()+theme(legend.position = 'top')

ggsave(here('graphs/cascades/F_SOC_DoA.png'),w=5,h=6)


fwrite(CS9[type=='model',.(location,mid,lo,hi)],
       file=here('graphs/cascades/data/F_SOC_DoA.csv'))


## prediction of alpha


## explore
SP2 <- D[,.(location,Arm,Facility,Age,Screened,Treated)]


ggplot(SP2,aes(Arm,Treated/Screened,shape=Age,lty=Facility,col=location,
                          group=paste(Facility,location,Age)))+
  geom_point()+geom_line()+
  facet_wrap(~Facility,scales='free')+
  theme_classic()+ggpubr::grids()


## aggregate by Age
RSP2 <- SP2[Arm!='SOC',.(Treated=sum(Treated),Screened=sum(Screened)),by=.(location)]
RSP2[,type:='data']
RSP2[,c('mid','lo','hi'):=MLH(Treated,Screened)]

mdl10 <- brm(formula=Treated | trials(Screened) ~1+(1|location),
               data=RSP2[location !='All'],
               family=binomial(link='logit'))


## extract data
SM10 <- bayessmy(mdl10)

## joing and plot both
CS10 <- rbind(RSP2[,.(type,location,mid,lo,hi)],SM10)

ggplot(CS10,aes(location,y=mid,ymin=lo,ymax=hi,group=paste(type),col=type))+
  scale_y_continuous(label=percent)+
  geom_pointrange(position=position_dodge(0.1),shape=1)+
  coord_flip()+
  ylab('Intervention D/B') + theme_classic()+ggpubr::grids()+theme(legend.position = 'top')

ggsave(here('graphs/cascades/F_IND_BoA.png'),w=5,h=6)


fwrite(CS10[type=='model',.(location,mid,lo,hi)],
       file=here('graphs/cascades/data/F_INT_BoA.csv'))


## implied alpha!

A <- merge(CS9[type=='model',.(location,midN=mid,VN=(hi-lo)^2/3.92^2)],
           CS10[type=='model',.(location,midD=mid,VD=(hi-lo)^2/3.92^2)],
           by='location')

A <- merge(CS9[type=='model',.(location,midN=mid,loN=lo,hiN=hi)],
           CS10[type=='model',.(location,midD=mid,loD=lo,hiD=hi)],
           by='location')

A[,c('Nmu','Nsig'):=.(logit(midN),(logit(hiN)-logit(loN))/3.92)]
A[,c('Dmu','Dsig'):=.(logit(midD),(logit(hiD)-logit(loD))/3.92)]

A[,mu:=Nmu-Dmu]
A[,sig:=sqrt(Nsig^2+Dsig^2)]

A[,alpha:=ilogit(mu)]
A[,alpha.lo:=ilogit(mu - 1.96*sig)]
A[,alpha.hi:=ilogit(mu + 1.96*sig)]


ggplot(A,aes(location,y=alpha,ymin=alpha.lo,ymax=alpha.hi))+
  scale_y_continuous(label=percent,limits=c(0,1))+
  geom_pointrange(shape=1)+
  coord_flip()+
  ylab('Modelled SOC B/A') + theme_classic()+ggpubr::grids()

ggsave(here('graphs/cascades/F_SOC_BoA_modelled.png'),w=5,h=6)


fwrite(A[,.(location,alpha,alpha.lo,alpha.hi)],
       file=here('graphs/cascades/data/F_SOC_BoA_modelled.csv'))
