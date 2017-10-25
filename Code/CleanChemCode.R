########################################################################
# This script is the chemistry and community analysis for the MASCOT Project:
# Here, I calculate net ecosystem calcification and net community production rates from 
# water samples collected in the summer of 2016.
# 
#Created by Dr. Nyssa Silbiger
#Created on 7/19/2016
#edited on 10/23/2017
#########################################################################



## clear workspace--------------
rm(list=ls())

#load libraries-----------
library(seacarb)
library(vegan)
library(plyr)
library(tidyr)
library(lmerTest)
library(lme4)
library(nls2)
library(usdm)
library(MuMIn)
library(MESS)
library(nlme)
library(r2glmm)
library(effects)


# set working directory -------------------------------
#setwd(~'/Biolophysical_feedbacks_in_coastal_ecosystems/Code')

## functions------------------------------------------
#Function for calculating fugosity of CO2
FCo2Calc<-function(Temp=25,Sal=35,u=0, CO2.water, CO2.air=400, rho=1025){
  
  # function for calculating fCO2
  #FCO2<-k*s*rho*(CO2.water-CO2.air)
  # k is gas transfer velocity k(600) = (0.266 ± 0.019)u10^2 (u10 is m s-1 and k(600 is cm h-1))
  # s is solubility of CO2 in seawater (K0 in seacarb)
  # CO2 water is CO2 in water
  # CO2 air is CO2 in air
  
  # calculate gas transfer velocity from Ho et al. 2006 parameterization
  k<-(0.266*u^2)/100 # m hr-1
  # calculate solubility of CO2 from temp and salinity from seacarb (Weiss et al. 1974)
  s<-as.numeric(K0(S=Sal, T=Temp,P=0, Patm=1))/1000 #mmol kg-1 uatm-1
  #calculate FCO2
  FCO2<-(k*s*rho*(CO2.water-CO2.air)) #fCO2 is mmol m-2 hr-1
  
  return(FCO2)
}

# function to rotate a matrix
rotate <- function(x) t(apply(x, 2, rev))

#load data------------------------------------------
#read in the light data

par(mfrow=c(4,2))
source('ImportTempCode.R')
par(mrow=c(1,1))

#Chemistry data
CData<-read.csv('../Data/RawData/ChemDataAll.csv')

#pool descriptions
PData<-read.csv('../Data/RawData/PoolDescriptionsAll.csv')

# Data Analysis------------------------------------------

#Run CO2Sys to get all carbonate parameters-------------------

#First calculate salinity from conductivity
CData$SalCal<-swSCTp(as.numeric(CData$Conductivity)/1000, CData$Temp.pool, pressure=rep(0, each=nrow(CData)),
                     conductivityUnit="mS/cm")

#plot salinity straight from instrument versus calculated from conuctivity and temp
plot(CData$SalCal, CData$Salinity, xlab='Salinity from conductivity', ylab="Salinity")
#looks very similar...good


##CO2 Calculations--------------------------------------------

## there is missing data because of monterey pools (only 13 pools with data instead of 15). replace NA for salinity and temp with 35 and 25 just so the code will run.
#This will not effect the data or any of the analysis because I am not using any data from those pools
CData$SalCal[is.na(CData$SalCal)]<-35
CData$Temp.pool[is.na(CData$Temp.pool)]<-25
CData$Temp.in[is.na(CData$Temp.in)]<-25
CData$Salinity[is.na(CData$Salinity)]<-35

#calculate pH at in situ temperature
CData$pH<-pHinsi(pH=CData$pH.pre, ALK=CData$TA, Tinsi=CData$Temp.pool, Tlab=CData$Temp.in, S=CData$Salinity,
                 Pt=CData$PO/1000000, Sit=0, k1k2 = "x",  kf = "x", ks="d", pHscale = "T", b="u74")

#correct the TA for the calibration error
CData$TA<-CData$TA-2235*(CData$CRM.off/100)

#pool 7 for CDM  ia s crazy outlier for TA....  We think this pool has an underground "river" and was not sealed off. We have therefore removed it from the dataset
CData$TA.Norm[CData$Pool==7 & CData$Site=='CDM'] <-NA
CData$TA[CData$Pool==7 & CData$Site=='CDM'] <-NA

#calculate all the CO2Sys params
PoolCO2<-carb(flag=8, CData$pH, CData$TA/1000000, S=CData$Salinity, T=CData$Temp.pool, Patm=1, P=0, Pt=CData$PO/1000000, Sit=0,
              k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential")
#TA is divided by 1000 because all calculations are in mol/kg in the seacarb package

#convert CO2, HCO3, CO3, DIC, and Alk back to micromol for easier interpretation
PoolCO2[,c("CO2","HCO3","CO3","DIC","ALK")]<-PoolCO2[,c("CO2","HCO3","CO3","DIC","ALK")]*1000000

#combine all the data
CData<-cbind(CData,PoolCO2[,-1:-6])

# calculate the fugosity of CO2 from the wind data
CData$FugCO2<-FCo2Calc(Temp = CData$Temp.pool, Sal = CData$Salinity, u = CData$Wind ,CO2.water = CData$pCO2insitu)

#Normalize the TA to nutrients and salinity------------------------------------
CData$TA.N.Norm<-CData$TA+CData$NN+CData$PO-CData$NH4 #nutrient normalized

#Normalize TA salinity after accounting for nutrients
CData$TA.Norm<-CData$TA.N.Norm*CData$Salinity/36
#TA normalized to constant salinity, but not nutrients here
CData$TA.NormSal<-CData$TA*CData$Salinity/36

#Normalize DIC constant salinity
CData$DIC.Norm<-CData$DIC*CData$Salinity/36

#convert data to date time
CData$DateTime<-strptime(paste(CData$Date, CData$Time), "%m/%d/%Y %H:%M")
CData$Time<-strptime( CData$Time, "%H:%M")

#round the times to the nearest minute
CData$DateTime<-format(round(CData$DateTime, units="mins"), format="%Y-%m-%d %H:%M")

##for some reason the dates are showing up as characters... this converts back to datetime
CData$DateTime<-strptime(CData$DateTime, '%Y-%m-%d %H:%M')

#make an array with the unique pool names for the for loops
Pools<-unique(CData$Pool)

# for loops for each site
CData$Site<-factor(CData$Site, levels=c("CDM","Monterey","Bodega","BobCreek")) #put the factors in order of latitude
Sites<-unique(CData$Site) # site names
DN<-c('Day','Night') # day versys night 

# Add in the light data at each time point----------------------------
#find the light environment during the actual sampling time
lightTimesBodega<-BodegaPAR[BodegaPAR$DateTime%in%CData$DateTime,] 
lightTimesBC<-BobCreekPAR[BobCreekPAR$DateTime%in%CData$DateTime,]
lightTimesCDM<-CDMPAR[CDMPAR$DateTime%in%CData$DateTime,]
lightTimesMonterey<-MontereyPAR[MontereyPAR$DateTime%in%CData$DateTime,]
#CDM was missing light data from pool 1 and 5 (Poseidon stole the loggers)... replace with data from pool 4, the closest pool
lightTimesCDM$TP1_light<-lightTimesCDM$TP4_light
lightTimesCDM$TP5_light<-lightTimesCDM$TP4_light
lightTimesCDM<-lightTimesCDM[,c(1,15,2:4,16,5:14)] #put back in order
#PAR
CDMPAR$TP1_light<-CDMPAR$TP4_light
CDMPAR$TP5_light<-CDMPAR$TP4_light
CDMPAR<-CDMPAR[,c(1,15,2:4,16,5:14)] #put back in order

#temp
CDMTemp$TP1<-CDMTemp$TP4
CDMTemp$TP5<-CDMTemp$TP4
CDMTemp<-CDMTemp[,c(1,15,2:4,16,5:14)] #put back in order

# calculate integrated PAR measurement for each pool from full timeseries data
# from start to endtime during the day

#first extract just the data from the daytime sampling period
BodegaDayPAR.all<-BodegaPAR[BodegaPAR$DateTime> lightTimesBodega[1,1] & BodegaPAR$DateTime< lightTimesBodega[96,1],]
BobCreekDayPAR.all<-BobCreekPAR[BobCreekPAR$DateTime> lightTimesBC[1,1] & BobCreekPAR$DateTime< lightTimesBC[96,1],]
CDMDayPAR.all<-CDMPAR[CDMPAR$DateTime> lightTimesCDM[1,1] & CDMPAR$DateTime< lightTimesCDM[96,1],]
MontereyDayPAR.all<-MontereyPAR[MontereyPAR$DateTime> lightTimesMonterey[1,1] & MontereyPAR$DateTime< lightTimesMonterey[96,1],]

# using auc (area under curve) in MESS to calculate an integrative measure of PAR
PAR.Int.Bodega<-NA
PAR.Int.BobCreek<-NA
PAR.Int.CDM<-NA
PAR.Int.Monterey<-NA

#Bodega
for (i in 1:15){# loop through all the pools
  PAR.Int.Bodega[i]<-auc(as.numeric(BodegaDayPAR.all[,1]),BodegaDayPAR.all[,i+1])
}
#BobCreek
for (i in 1:15){
  PAR.Int.BobCreek[i]<-auc(as.numeric(BobCreekDayPAR.all[,1]),BobCreekDayPAR.all[,i+1])
}
#CDM
for (i in 1:15){
  PAR.Int.CDM[i]<-auc(as.numeric(CDMDayPAR.all[,1]),CDMDayPAR.all[,i+1])
}
#Monterey
for (i in 1:15){
  PAR.Int.Monterey[i]<-auc(as.numeric(MontereyDayPAR.all[,1]),MontereyDayPAR.all[,i+1])
}

## add integrated PAR to PData
PData$PAR.Int<-c(PAR.Int.Bodega, PAR.Int.BobCreek, PAR.Int.CDM, PAR.Int.Monterey)

# put the light data for each time point in CData
CData$PAR<-rep('NA',nrow(CData))
#add a new column with the PAR data
for (i in 1:15){
  CData$PAR[CData$Pool==i & CData$Site=='Bodega']<-lightTimesBodega[CData$Pool==i,i+1]
  CData$PAR[CData$Pool==i & CData$Site=='BobCreek']<-lightTimesBC[CData$Pool==i,i+1]
  CData$PAR[CData$Pool==i & CData$Site=='CDM']<-lightTimesCDM[CData$Pool==i,i+1]
  CData$PAR[CData$Pool==i & CData$Site=='Monterey']<-lightTimesMonterey[CData$Pool==i,i+1]
}

CData$PAR<-as.numeric(CData$PAR) #converts to numeric because of the NAs

### For loop to calculate the change in each parameter between each time point--------------

# make all nighttime data 0 for light
CData$PAR[CData$Day.Night=='Night']<-0

#create a dataframe that calculates the change in each parameter between two consecutive time points----------------
time<-c(0:5)

## pH, O2, and temp metrics--------------

#calculate differet pH metrics for each pool
maxpH<-NA
rangepH<-NA
meanpH<-NA
maxdiffpH<-NA # max difference between ocean and tide pool pH over the 24 hours cycle
s<-c(0,15,30,45) #add to i for each site... number of pools
for (j in 1:length(Sites)){
  for (i in 1:15){
    rangepH[i+s[j]]<-max(CData$pH[CData$Pool==i& CData$Site==Sites[j]], na.rm = TRUE)-min(CData$pH[CData$Pool==i& CData$Site==Sites[j]], na.rm = TRUE)
    maxpH[i+s[j]]<-max(CData$pH[CData$Pool==i& CData$Site==Sites[j]], na.rm = TRUE)
    meanpH[i+s[j]]<-mean(CData$pH[CData$Pool==i& CData$Site==Sites[j]], na.rm = TRUE)
    maxdiffpH[i+s[j]]<-max(CData$pH[CData$Pool==i& CData$Site==Sites[j]] - CData$pH[CData$Pool=='Ocean' & CData$Site==Sites[j]], na.rm = TRUE)
  }
}

#calculate change in pH for day and night 
pH.night.mean<-NA
pH.day.mean<-NA
pH.day.max<-NA # max difference between ocean and tide pool during day
pH.night.max<-NA # max difference between ocean and tide pool during the night
for (j in 1:length(Sites)){
  for (i in 1:15){

    for (k in 1:5){
      pH.day.mean[i+s[j]]<-mean(CData$pH[CData$Pool==i& CData$Site==Sites[j] & CData$Time.pt==k & CData$Day.Night=='Day'], na.rm=TRUE)
      pH.day.max[i+s[j]]<-max(CData$pH[CData$Pool==i& CData$Site==Sites[j] & CData$Time.pt==k & CData$Day.Night=='Day'] -
                                CData$pH[CData$Pool=='Ocean'& CData$Site==Sites[j] & CData$Time.pt==k & CData$Day.Night=='Day'], na.rm=TRUE)
      
    }
  }
}
# change in pH in a pool during the night 
j<-3 # site 3 #CDM
for (i in 1:15){ # loop through pools
  #CDM night only had 4 time points instead of 5 so I need to do this for loop twice   
  
  for (k in 1:4){ # loop through time points
    pH.night.mean[i+s[j]]<-mean(CData$pH[CData$Pool==i& CData$Site==Sites[j] & CData$Time.pt==k & CData$Day.Night=='Night'], na.rm=TRUE)
    pH.night.max[i+s[j]]<-max(CData$pH[CData$Pool==i& CData$Site==Sites[j] & CData$Time.pt==k & CData$Day.Night=='Night'] -
                                CData$pH[CData$Pool=='Ocean'& CData$Site==Sites[j] & CData$Time.pt==k & CData$Day.Night=='Night'], na.rm=TRUE)
    
    
  }
}
for (j in c(1,2,4)){ # loop through sites
  for (i in 1:15){ # loop through pools
      for (k in 1:5){ # loop through time points
      pH.night.mean[i+s[j]]<-mean(CData$pH[CData$Pool==i& CData$Site==Sites[j] & CData$Time.pt==k & CData$Day.Night=='Night'], na.rm=TRUE)
      pH.night.max[i+s[j]]<-max(CData$pH[CData$Pool==i& CData$Site==Sites[j] & CData$Time.pt==k & CData$Day.Night=='Night'] -
                                  CData$pH[CData$Pool=='Ocean'& CData$Site==Sites[j] & CData$Time.pt==k & CData$Day.Night=='Night'], na.rm=TRUE)
      
    } 
  }
}

#create a dataframe with all the pH metrics
pHmetrics<-data.frame(PData$Site, rangepH, maxpH, meanpH, pH.day.mean, pH.night.mean, maxdiffpH)
colnames(pHmetrics)[1]<-'Site'

#because monterey is missing two pools I need to convert the inf values to NA
is.na(pHmetrics) <- sapply(pHmetrics, is.infinite)

## Temperature metrics
rangeT<-NA
meanT<-NA
s<-c(0,15,30,45) #add to i for each site... number of pools
for (j in 1:length(Sites)){
  for (i in 1:15){
    rangeT[i+s[j]]<-max(CData$Temp.pool[CData$Pool==i& CData$Site==Sites[j]], na.rm = TRUE)-min(CData$Temp.pool[CData$Pool==i& CData$Site==Sites[j]], na.rm = TRUE)
    meanT[i+s[j]]<-mean(CData$Temp.pool[CData$Pool==i& CData$Site==Sites[j]], na.rm = TRUE)
    
  }
}
#calculate mean temp for day and night 
Temp.night.mean<-NA
Temp.day.mean<-NA
for (j in 1:length(Sites)){
  for (i in 1:15){
    for (k in 1:5){
      
      Temp.day.mean[i+s[j]]<-mean(CData$Temp.pool[CData$Pool==i& CData$Site==Sites[j] & CData$Time.pt==k & CData$Day.Night=='Day'], na.rm=TRUE)
      
    }
  }
}

#  mean change in Temp in a pool during the night
j<-3
for (i in 1:15){
  #CDM night only had 4 time points so I need to do this   
  for (k in 1:4){
    Temp.night.mean[i+s[j]]<-mean(CData$Temp.pool[CData$Pool==i& CData$Site==Sites[j] & CData$Time.pt==k & CData$Day.Night=='Night'], na.rm=TRUE)
    
  }
}
for (j in c(1,2,4)){
  for (i in 1:15){
    
    for (k in 1:5){
      Temp.night.mean[i+s[j]]<-mean(CData$Temp.pool[CData$Pool==i& CData$Site==Sites[j] & CData$Time.pt==k & CData$Day.Night=='Night'], na.rm=TRUE)
      
    } 
  }
}


#create a dataframe with all the Temp metrics
Tempmetrics<-data.frame(PData$Site, rangeT, meanT,Temp.night.mean,Temp.day.mean)
colnames(Tempmetrics)[1]<-'Site'

#because monterey is missing two pools I need to convert the inf values to NA
is.na(Tempmetrics) <- sapply(Tempmetrics, is.infinite)

# calculate change in all parameters between two time ponts
param<-colnames(CData)[c(7,9,10,16:18,23:28, 33,34,37,38,42,44,43,39)] # all parameters I am interested in
delta<-as.data.frame(x=matrix(nrow=640,ncol=length(param)+3)) #pre-allocate space for the new data frame
colnames(delta)<-c('Pool', 'Time.pt','Day.Night') #pre-name the first 3 columns for each pool, time point, and whether it was taken during the day or night
delta$Pool<-rep(Pools,40) #Pool ID
delta$Time.pt<-rep(c(rep(1,16),rep(2,16),rep(3,16),rep(4,16),rep(5,16)),8) #time points
delta$Day.Night<-rep(c(rep('Day',80),rep('Night',80)),4) # day or night
delta$Site<-c(rep('Bodega',160),rep('BobCreek',160),rep('CDM',160),rep('Monterey',160)) #site names
delta$Site<-factor(delta$Site, levels=c("CDM","Monterey","Bodega","BobCreek"))  #make the site names factors

for (m in 1:length(Sites)){ #loops through the sites
  for (l in 1:length(param)){ #loop through parameters
    for (k in 1:2){ #loop through day night
      for (i in 1:16){ #loop through pools
        for (j in 1:(length(unique(CData$Time.pt))-1)){ #loop through sites
          
# time point 2 minus 1 and so on
          delta[delta$Time.pt==time[j+1] & delta$Pool==Pools[i] & delta$Day.Night==DN[k] & delta$Site==Sites[m], l+3]<-
            CData[CData$Pool==Pools[i] & CData$Time.pt==j & CData$Day.Night==DN[k] & CData$Site==Sites[m], param[l]]-
            CData[CData$Pool==Pools[i] & CData$Time.pt==j-1 & CData$Day.Night==DN[k]& CData$Site==Sites[m], param[l]]
          colnames(delta)[l+3]<-param[l] #change the column names appropriately
        }
      }
    }
  }
}

#calculate the average FCO2 between two sampling points.  This is what I will use for the Fugosity of CO2 in the NCP calculations
#replace NA for FUgCO2 with a 0
CData$FugCO2[is.na(CData$FugCO2)]<-0
for (m in 1:length(Sites)){ #loops through the sites
  for (k in 1:2){ #loop through day night
    for (i in 1:16){ #loop through pools
      for (j in 1:(length(unique(CData$Time.pt))-1)){ #loop through sites
        
        delta$FugCO2[delta$Time.pt==time[j+1] & delta$Pool==Pools[i] & delta$Day.Night==DN[k] & delta$Site==Sites[m]]<-mean(
          CData$FugCO2[CData$Pool==Pools[i] & CData$Time.pt==j & CData$Day.Night==DN[k] & CData$Site==Sites[m]],
          CData$FugCO2[CData$Pool==Pools[i] & CData$Time.pt==j-1 & CData$Day.Night==DN[k]& CData$Site==Sites[m]])
      }
    }
  }
}

# All params without datetime
param2<-colnames(CData)[c(7,9,10,16:18,23:28, 33,34,37,38,42,44,39)]

# make a dataframe that calculates the mean between consecutive time points so that I can plot rates versus actual values----
MeanData<-as.data.frame(x=matrix(nrow=640,ncol=length(param2)+3)) #pre-allocate space for the new data frame
colnames(MeanData)<-c('Pool', 'Time.pt','Day.Night') #pre-name the first 3 columns for each pool, time point, and whether it was taken during the day or night
MeanData$Pool<-rep(Pools,40) #Pool ID
MeanData$Time.pt<-rep(c(rep(1,16),rep(2,16),rep(3,16),rep(4,16),rep(5,16)),8) #time points
MeanData$Day.Night<-rep(c(rep('Day',80),rep('Night',80)),4) # day or night
MeanData$Site<-c(rep('Bodega',160),rep('BobCreek',160),rep('CDM',160),rep('Monterey',160)) #site names
MeanData$PoolID<-CData$PoolID[CData$Time.pt != 5]
MeanData$SetID<-CData$SetID[CData$Time.pt != 5]

for (m in 1:length(Sites)){ #loops through the sites
  for (l in 1:length(param2)){ #loop through parameters
    for (k in 1:2){ #loop through day night
      for (i in 1:16){ #loop through pools
        for (j in 1:(length(unique(CData$Time.pt))-1)){ #loop through sites
          

          MeanData[MeanData$Time.pt==time[j+1] & MeanData$Pool==Pools[i] & MeanData$Day.Night==DN[k] & MeanData$Site==Sites[m], l+3]<-
            mean(CData[CData$Pool==Pools[i] & CData$Time.pt==j & CData$Day.Night==DN[k] & CData$Site==Sites[m], param2[l]],
                 CData[CData$Pool==Pools[i] & CData$Time.pt==j-1 & CData$Day.Night==DN[k]& CData$Site==Sites[m], param2[l]], na.rm=TRUE)
          colnames(MeanData)[l+3]<-param2[l] #change the column names appropriately
        }
      }
    }
  }
}

#remove Ocean from data
MeanData<-data.frame(MeanData[-c(which(MeanData$Pool=='Ocean')),])

# order the factors of site for Mean Data
MeanData$Site<-factor(MeanData$Site, levels=c("CDM","Monterey","Bodega","BobCreek")) #put the factors in order of latitude

# order the sites in the PData dataframe
PData$Site<-factor(PData$Site, levels=c("CDM","Monterey","Bodega","BobCreek"))
#add a column for change in time in the appropriate units
delta$TimeDiff<-NULL
for (m in 1:length(Sites)){
  for (k in 1:2){
    for (i in 1:16){
      for (j in 1:(length(unique(CData$Time.pt))-1)){
        
        
        delta$TimeDiff[delta$Time.pt==time[j+1] & delta$Pool==Pools[i] & delta$Day.Night==DN[k]& delta$Site==Sites[m]]<-
          difftime(CData$DateTime[CData$Pool==Pools[i] & CData$Time.pt==j & CData$Day.Night==DN[k]& CData$Site==Sites[m]],CData$DateTime[CData$Pool==Pools[i] & CData$Time.pt==j-1 & CData$Day.Night==DN[k]& CData$Site==Sites[m]], units="hours")
      }
    }
  }
}


# Metabolism calculations----------------------------
#NEC rates in mmol m-2 hr-1
#convert rates by multiplying by density and volume and dividing by surface area and time. 

#remove Ocean from delta
delta<-delta[-c(which(delta$Pool=='Ocean')),]

#convert volume from L to m3 
PData$Vol<-PData$Vol*0.001

VolRemoved<-c(.4,0.8,1.2,1.6, 2.0)*0.001 #how much water was removed by each sampling point... need to substract this to get the actual volume
#at each sampling time.  We removed 400ml of water each time. 

delta$Vol<-NULL
Volume<-PData$Vol
# calculate volume of pools after water was removed during sampling points
for (j in 1:length(Sites)){
  for (i in 1:5){
    delta$Vol[delta$Time.pt==i & delta$Site==Sites[j]]<-Volume[PData$Site==Sites[j]] - VolRemoved[i]
  }
}

# repeated the surface area of the pools
SA<-c(rep(PData$SA[1:15],10),rep(PData$SA[16:30],10),rep(PData$SA[31:45],10),rep(PData$SA[46:60],10))

# calculate all the rates 
#NEC delta TA/2 * seawater density * volume/ surface area / time.  Divided by 100 so mmols m2 hr and multiplied by negative 1 so pos values are calcification
delta$NEC<--1*(delta$TA.Norm/2)*1023*(delta$Vol/SA)*(1/delta$TimeDiff)/100
#NCP delta DIC * SW density * volume/SA / time  -  NEC + Fugosity of CO2
delta$NCP<-(((-delta$DIC.Norm)*1023*(delta$Vol/SA)*(1/delta$TimeDiff)/100)-delta$NEC +delta$FugCO2) #positive values are net photosynthesis (the equation is - fugosity, but because it is time 2-1 it is plus here)
delta$Pool<-as.numeric(levels(delta$Pool))[delta$Pool] # make pool numbers numeric

# average pH of the ocean sample for each site/day night
OceanData<-ddply(CData, c("Site", "Day.Night"), summarize,
                 OceanpH = mean(pH[CData$Pool=='Ocean'], na.rm=TRUE),
                 OceanTemp = mean(Temp.pool[CData$Pool=='Ocean'], na.rm=TRUE),
                 OceanNH4 = mean(NH4[CData$Pool=='Ocean'], na.rm=TRUE),
                 OceanNN = mean(NN[CData$Pool=='Ocean'], na.rm=TRUE)
)

# average ocean just by site
OceanDataMean<-ddply(CData, c("Site"), summarize,
                     OceanpH = mean(pH[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanMinpH = min(pH[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanMaxpH = max(pH[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanTemp = mean(Temp.pool[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanNH4 = mean(NH4[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanNN = mean(NN[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanPO = mean(PO[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanDO = mean(DO[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanArag = mean(OmegaAragonite[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanpCO2 = mean(pCO2insitu[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanDIC = mean(DIC.Norm[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanTA = mean(TA.NormSal[CData$Pool=='Ocean'], na.rm=TRUE)
)

# make all the data for pool 4 and 9 for monterey = NA and ocean time point 5 for CSM... the temperature data is still showing up even though it is not real (no temp data was taken)
CData[CData$Site=='Monterey' & CData$Pool==4 | CData$Pool==9,'Temp.pool' ]<-NA
CData[CData$Site=='CDM' & CData$Time.pt==5 & CData$Day.Night=='Night','Temp.pool' ]<-NA

# ocean data range
OceanDataRange<-ddply(CData, c("Site","Day.Night"), summarize,
                     
                     OceanMinpH = min(pH[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanMaxpH = max(pH[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanMinTemp = min(Temp.pool[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanMaxTemp = max(Temp.pool[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanMinNH4 = min(NH4[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanMaxNH4 = max(NH4[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanMinNN = min(NN[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanMaxNN = max(NN[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanMinPO = min(PO[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanMaxPO = max(PO[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanMinDO = min(DO[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanMaxDO = max(DO[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanMinArag = min(OmegaAragonite[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanMaxArag = max(OmegaAragonite[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanMinpCO2 = min(pCO2insitu[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanMaxpCO2 = max(pCO2insitu[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanMinDIC = min(DIC.Norm[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanMaxDIC = max(DIC.Norm[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanMinTA = min(TA.NormSal[CData$Pool=='Ocean'], na.rm=TRUE),
                     OceanMaxTA = max(TA.NormSal[CData$Pool=='Ocean'], na.rm=TRUE)
)


# Pool ranges for table (excluding ocean sample)
PoolDataRange<-ddply(CData, c("Site","Day.Night"), summarize,
                     pHmin = min(pH[CData$Pool!='Ocean'], na.rm=TRUE),
                     pHmax = max(pH[CData$Pool!='Ocean'], na.rm=TRUE),
                     Tempmin = min(Temp.pool[CData$Pool!='Ocean'], na.rm=TRUE),
                     Tempmax = max(Temp.pool[CData$Pool!='Ocean'], na.rm=TRUE),
                     DOmin = min(DO[CData$Pool!='Ocean'], na.rm=TRUE),
                     DOmax = max(DO[CData$Pool!='Ocean'], na.rm=TRUE),
                     Aragmin = min(OmegaAragonite[CData$Pool!='Ocean'], na.rm=TRUE),
                     Aragmax = max(OmegaAragonite[CData$Pool!='Ocean'], na.rm=TRUE),
                     pCO2min = min(pCO2insitu[CData$Pool!='Ocean'], na.rm=TRUE),
                     pCO2max = max(pCO2insitu[CData$Pool!='Ocean'], na.rm=TRUE),
                     DICmin = min(DIC.Norm[CData$Pool!='Ocean'], na.rm=TRUE),
                     DICmax = max(DIC.Norm[CData$Pool!='Ocean'], na.rm=TRUE),
                     TAmin = min(TA.NormSal[CData$Pool!='Ocean'], na.rm=TRUE),
                     TAmax = max(TA.NormSal[CData$Pool!='Ocean'], na.rm=TRUE)
)
#ranges for the physical data
PoolPhysDataRange<-ddply(PData, c("Site"), summarize,
                         Permmin = min(Perimeter, na.rm=TRUE),
                         Permmax = max(Perimeter, na.rm=TRUE),
                         Depthmin = min(MaxDepth, na.rm=TRUE),
                         Depthmax = max(MaxDepth, na.rm=TRUE),
                         SAmin = min(SA, na.rm=TRUE),
                         SAmax = max(SA, na.rm=TRUE),
                         Volmin = min(Vol, na.rm=TRUE),
                         Volmax = max(Vol, na.rm=TRUE)
                         
)

#calculate Ta/DIC slope for each pool-----------------
cols<-rainbow(15)
PData$slopes.TADIC<-NA
s<-c(0,15,30,45)
par(mar=c(5, 4, 4, 2) + 0.1, mfrow=c(1,1)) #return the margins back to default
for (j in 1:length(Sites)){
  #because there is an uneven number of pools at each site.. this removes the pools that have no data
  n<-unique(delta$Pool[delta$Site==Sites[j] & !is.na(delta$NEC)])
  
  for (i in n){
    #slopes for  normalized TA vs DIC
    x<-CData$DIC.Norm[delta$Pool==i & delta$Site==Sites[j]]
    y<-CData$TA.Norm[delta$Pool==i & delta$Site==Sites[j]]
    model4<-lm(y~x, na.action=na.exclude)
    PData$slopes.TADIC[i+s[j]]<-model4$coefficients[2] #save the slopes
    
  }
  
}

#plot TA vs DIC and color by either aragonie or pH-----------
source('ContourFunction.R') #source my contour plot function
# this is to make the background of the property-property plot aragonite saturation state
AT <- seq(1300e-6, 2750e-6, length.out=200) # create TA data
DIC <- seq(850e-6, 2900e-6, length.out=200) # create DIC data

dat <- expand.grid(AT, DIC)
carb <- carb(flag=15, var1=dat$Var1, var2=dat$Var2, S=35, T=25, P=0, Pt=0, Sit=0, k1k2="l", kf="pf", pHscale="T")
arag<-matrix(carb$OmegaAragonite, 200)
lev=20 #levels for the contour

pdf(file = '../Output/TAvsDICPropertyPlots.pdf', height = 5, width = 7)
par(mar=c(5, 4, 4, 2) + 0.1, mfrow=c(1,2)) #return the margins back to default
for (j in 1:length(Sites)){
  
  #because there is an uneven number of pools at each site.. this removes the pools that have no data
  n<-(unique(CData$Pool[CData$Site==Sites[j] & !is.na(CData$TA.Norm)]))
  n<- as.numeric(as.character(n[n!='Ocean'])) # remove the ocean sample
  my.filled.contour(DIC*1e6,AT*1e6, rotate(rotate(arag)), 
                    xlab="DICn",
                    ylab="TAn",
                    xlim= c(min(CData$DIC.Norm,na.rm=TRUE)-20, max(CData$DIC.Norm,na.rm=TRUE)+10),
                    ylim= c(min(CData$TA.Norm,na.rm=TRUE)-20, max(CData$TA.Norm,na.rm=TRUE)+10),
                    zlim = c(0,16),
                    key.title = title(main = expression(paste(Omega[arg]))),
                    levels=c(seq(0,3.5,0.1), seq(4,15,1)),
                    lwd=2,
                    lty="solid",
                    plot.title ={for (i in n){ 
                      x<-CData$DIC.Norm[CData$Pool==i & CData$Site==Sites[j]]
                      y<-CData$TA.Norm[CData$Pool==i & CData$Site==Sites[j]]
                      model1<-lm(y~x, na.action=na.exclude)
                      CI<-predict(model1, interval="confidence")
                      CI.order<-order(CI[,1])
                      lines(x[CI.order],CI[CI.order,1], col='black')
                      points(x, y, col='black', pch=19)};axis(1); axis(2)},
                    
                    plot.axes={ contour(DIC*1e6,AT*1e6, rotate(rotate(arag)), nlevels = 15,
                                        levels = c(0,0.25,0.5,0.75,1,2,3,4,5,6,7),          
                                        drawlabels = TRUE, axes = FALSE, 
                                        frame.plot = FALSE, add = TRUE, method='flattest',
                                        col='grey'

                    )}

                    
  )
  title(main = Sites[j], cex=1.5, adj=0.5, outer = TRUE, line=-1.5)
}

dev.off()
#steeper slopes are more calcification driven and flatter slopes are more photosynthesis driven

#Calculate relative tide heights
PData$TidePercent<-NA
PData$TidePercent[PData$Site=='BobCreek']<-100*PData$TideHeight[PData$Site=='BobCreek']/max(PData$TideHeight[PData$Site=='BobCreek'])
PData$TidePercent[PData$Site=='Monterey']<-100*PData$TideHeight[PData$Site=='Monterey']/max(PData$TideHeight[PData$Site=='Monterey'], na.rm=T)
PData$TidePercent[PData$Site=='Bodega']<-100*PData$TideHeight[PData$Site=='Bodega']/max(PData$TideHeight[PData$Site=='Bodega'])
PData$TidePercent[PData$Site=='CDM']<-100*PData$TideHeight[PData$Site=='CDM']/max(PData$TideHeight[PData$Site=='CDM'])


#Add surface area to volume metric 
PData$SAV<-PData$SA/PData$Vol


#Mean NCP and NEC by pool

MeanNCPNEC<-ddply(delta, c('Site','Pool', 'Day.Night'), summarize,
                  NCP = mean(NCP, na.rm = TRUE),
                  NEC = mean(NEC, na.rm = TRUE))
is.na(MeanNCPNEC) <- sapply(MeanNCPNEC, is.infinite)

# Figure of pH versus DO
pdf('../Output/DOvspH.pdf', height = 6, width = 6, useDingbats=FALSE)
par(mfrow=c(1,1))

# get the colors right
cols.all2<-rep(NA, length(CData$Site))
for (i in 1:length(CData$Site)){
  if(CData$Site[i]=='CDM'){
    cols.all2[i]<- 'lightblue'
  }
  if (CData$Site[i]=='Monterey'){
    cols.all2[i]<-'#0606f7'}
  if (CData$Site[i]=='Bodega'){ 
    cols.all2[i]<-'magenta'}
  if (CData$Site[i]=='BobCreek'){
    cols.all2[i]<-'red'
  }
}

col.u<-c('magenta','red','lightblue','#0606f7')
plot(CData$DO, CData$pH, col =cols.all2, pch = 19, xlab = expression(paste('DO mg L'^ {-1})), ylab = 'pH' )

#MME for pH vs DO by site
pHDO<-lme(pH~DO, random = ~1|Site, data=CData, na.action = na.omit)
pHDO.f<-fitted(pHDO, level = 1)
# add a legend
legend("topleft",legend= c('Bob Creek','Bodega Bay','Monterey Bay','Corona del Mar'), col=c('red','magenta','#0606f7','lightblue'), pch=19, bty="n")
dev.off()

# calculate average NEC and NCP rates per site for supplemental table
aggregate(NCP~Site + Day.Night,data =  MeanNCPNEC, FUN = function(x) mean(x, na.rm=T))
aggregate(NEC~Site + Day.Night,data =  MeanNCPNEC, FUN = function(x) mean(x, na.rm=T))

##### Community data cleaning and analysis #####
#load data------------------------------------------
#read in the community data
Sessile<-read.csv('../Data/RawData/SessileAll_edited.csv', stringsAsFactors=FALSE) #Sessile are algae and sessile inverts in number of squares filled (to be convered to percent cover)
Mobile<-read.csv('../Data/RawData/MobileAll3.csv', stringsAsFactors=FALSE) #Mobile inverts are counts

# Pull out the info for the subgroups
SessileGroups<- Sessile[1,] #pull out the first row
InvertGroups <- Mobile[1,]

#remove the first row and convert everything back to numeric
Sessile<-data.frame(Sessile[-1,])
Mobile<-data.frame(Mobile[-1,])

#convert everything back to numeric
for (i in 3: ncol(Mobile)){
Mobile[,i]<-as.numeric(Mobile[,i])
}

for (i in 4: ncol(Sessile)){
  Sessile[, i]<-as.numeric(Sessile[,i])
}


#replace NA with 0
Sessile[is.na(Sessile)]<-0 
Mobile[is.na(Mobile)]<-0

# Make all the community data a relative percent
PercentSessile<-100*Sessile[4:ncol(Sessile)]/Sessile$Squares
#create % cover for the mobile inverts.
MobileCover<-Mobile

## Rescale the mobile inverts after measuring them to convert to percent cover
# these were calculated by outlining lots of individuals in ImageJ.  See supplemental methods from Silbiger and Sorte for details
MobileCover$Chlorostoma<-(Mobile$Chlorostoma*2.6)/100 # divide by 100 cm2 because that is how much one square is
MobileCover$Littorina <- (MobileCover$Littorina*0.28)/100 # Litorrina
MobileCover$Lottia.sp<- (MobileCover$Lottia.sp*1.84)/100 # limpets
MobileCover$Fissurella<- (MobileCover$Fissurella*1.84)/100 # limpets
MobileCover$Cyanoplax<- (MobileCover$Cyanoplax*1.19)/100 #chitons
MobileCover$Nuttalina<-(MobileCover$Nuttalina*1.19)/100 #chitons
MobileCover$Mopilia<-(MobileCover$Mopilia*1.19)/100 #chitons

MobileCover$Pagurus<- (MobileCover$Pagurus*2.6)/100 # most pagurus are in chlorostoma shells
# I did not measure Nucella, but they are slightly smaller than chlorostoma
MobileCover$Nucella<-(Mobile$Nucella*2)/100 # divide by 100 cm2 because that is how much one square is

# for the remaining groups
MobileCover[,which(InvertGroups=='Med')]<-(MobileCover[,which(InvertGroups=='Med')]*2)/100 # 2cm2
MobileCover[,which(InvertGroups=='Fish')]<-(MobileCover[,which(InvertGroups=='Fish')]*2)/100 # these aren't actually included in the analyses because there were so little
MobileCover[,which(InvertGroups=='Large')]<-(MobileCover[,which(InvertGroups=='Large')]*10)/100 # 10 cm2

#Add the Mobile cover to the sessile
PercentTotal<-cbind(PercentSessile,MobileCover[3:ncol(MobileCover)])
#normalize to the sum of the total cover (since it can be greater than 100%)
PercentTotal<- 100*PercentTotal/rowSums(PercentTotal)

## plot percent of pools by group for each site
##remove the site and square info from the groups and join them into one vector so its the same length as percent total
SessileGroups<-SessileGroups[-c(1:3)]
InvertGroups<- InvertGroups[-c(1:2)]
Groups<- cbind(SessileGroups, InvertGroups)

#divide by groups
AllAlgae<-rowSums(PercentTotal[,c(which(Groups=='Fleshy' | Groups=='Crust' | Groups=='Coralline'))]) #sum across the algae rows for total % cover
SurfGrass<- PercentTotal$Phyllospadix
FleshyAlgae<-rowSums(PercentTotal[,c(which(Groups=='Fleshy'))])
CorallineAlgae<-rowSums(PercentTotal[,c(which(Groups=='Coralline'))])
Crust<- PercentTotal$NonCorallinecCust
RockSand<- rowSums(PercentTotal[,c(which(Groups=='Rock'))])
Inverts<-rowSums(PercentTotal[,c(which(Groups=='Invert' | Groups=='Small' | Groups=='Med'| Groups=='Star'| Groups=='Large'))])#sum across the invert rows for total % cover
Fish<-rowSums(PercentTotal[,c(which(Groups=='Fish'))])
Mussel<-PercentTotal$Mytilus
OtherInverts<-Inverts-Mussel
# group them into one data frame
CoverbyGroups<-data.frame(Sessile$Site,RockSand,FleshyAlgae,CorallineAlgae,Crust,SurfGrass,Inverts,Fish, Mussel, OtherInverts)

colnames(CoverbyGroups)[1]<-'Site'

#average cover by site
Cover.mean <- ddply(CoverbyGroups, c("Site"), summarize,
                    RockSand = mean(RockSand, na.rm = T),
                    FleshyAlgae = mean(FleshyAlgae, na.rm = T),
                    CorallineAlgae = mean(CorallineAlgae, na.rm = T),
                    Crust = mean(Crust, na.rm = T),
                    SurfGrass = mean(SurfGrass, na.rm = T),
                    Inverts = mean(Inverts, na.rm = T),
                    Fish = mean(Fish, na.rm = T),
                    Mussel = mean(Mussel, na.rm = T),
                    OtherInvers = mean(OtherInverts, na.rm = T)
                    )
rownames(Cover.mean)<- Cover.mean$Site
Cover.mean<-Cover.mean[,-1]
site<-unique(Sessile$Site)

#sum across all tide pools at each site.
TotalCover<-rbind(colSums(PercentTotal[1:15,]),colSums(PercentTotal[16:30,]),colSums(PercentTotal[31:45,]),colSums(PercentTotal[46:60,]))
rownames(TotalCover)<-c('Bodega','BobCreek','CDM','Monterey')
TotalCover<-TotalCover[,-c(1:2)] #remove sand and rock

# col numbers for each group
Algae<-1:36
Inverts<-37:ncol(TotalCover)

#colors
cols.all<-c(rep('red',15),rep('blue',15),rep('green',15),rep('magenta',15))
# the community composition metric
b<-CoverbyGroups$FleshyAlgae+CoverbyGroups$SurfGrass-CoverbyGroups$Inverts

#PCA for physical variables
x<-data.frame(scale(PData[,c('Perimeter','MaxDepth','SA','Vol','TidePercent','SAV')], scale= apply(PData[,c('Perimeter','MaxDepth','SA','Vol','TidePercent','SAV')], 2, sd, na.rm = TRUE), 
                    center = TRUE))
PCA<- princomp( na.omit(x ))

#add back the NA so that I can use the PC scores in tests
PCA$scores<-rbind(PCA$scores[1:48,], rep(NA,ncol(PCA$scores)),PCA$scores[c(49:52),], rep(NA,ncol(PCA$scores)), PCA$scores[c(53:58),])

# to get the variance explained by each axis
s<-matrix(summary(PCA))

########## statistical analysis###############

## Figure 1
### plot pH time-series #####
pdf(file = '../Output/pHTimeSeries.pdf', width = 4.5, height = 7, useDingbats=FALSE)
par(mfrow=c(1,1), bg=NA)
par(mar=c(2.1,4.1,4.1,2.1))
DN<-c('Day','Night')
par(mfrow=c(4,2))

for (j in c(2,1,4,3)){
  for (k in 1:2){
    
    #because monterey has some missing pools I need to put this if statement in to pass over the missing ones
    pool<-if(Sites[j]=='Monterey') c(1:3,5:8,10:15) else 1:15
    plot(CData$DateTime[CData$Pool==Pools[16] & CData$Site==Sites[j] & CData$Day.Night==DN[k]], CData$pH[CData$Pool==Pools[16]& CData$Site==Sites[j] & CData$Day.Night==DN[k]], 
         xlab='Time', ylab = 'pH', ylim = c(7, 9), pch=19, cex.lab=1.5,
         col='grey', type='b', lwd=3, axes=F, cex=0.75)
    for (i in pool){
      points(CData$DateTime[CData$Pool==Pools[i]& CData$Site==Sites[j]& CData$Day.Night==DN[k]], CData$pH[CData$Pool==Pools[i]& CData$Site==Sites[j]& CData$Day.Night==DN[k]], 
             main = paste('Tidepool',Pools[i]), ylab = 'pH', ylim = c(7.1, 9.0), pch=19, type='b', cex=0.75)
    }
    axis(side=2, cex.axis=1)
    box(which = "plot", lty = "solid")
    axis.POSIXct(side = 1,x = CData$DateTime[CData$Pool==Pools[16] & CData$Site==Sites[j] & CData$Day.Night==DN[k]], cex.axis=1)
  }
  
}
dev.off()


#### Mixed effects models for mean and range pH ###
### Figure 2 ####
pdf(file = '../Output/pHEffectSizePlot.pdf', height = 7, width = 10, useDingbats=FALSE)

## Range pH data frame with scaled data
scaledData<-data.frame(scale(b, center=TRUE, scale = TRUE),scale(Tempmetrics$rangeT, center=TRUE, scale = TRUE),
                       scale(log(PData$PAR.Int), center=TRUE, scale = TRUE),scale(PData$slopes.TADIC, center=TRUE, scale = TRUE))

# rename the columns
colnames(scaledData)<-c('b','Temp','PAR','TADIC')

#remove missing data
no<-!is.na(PData$slopes.TADIC)

# look at the variance inflation factors to check for multicollinearity
vif(cbind(scaledData,PCA$scores[,1:2]))

# all are less than 2 , as long as < 2 then OK!

# scale the dependent variable as well
y<-scale(pHmetrics$rangepH, center=TRUE, scale = TRUE) # scale the pH data so that I can compare mean and range

## Mean DataFrame
#pHMean x DayNight
# separate by day and night for pH mean model
DataDayNight<-data.frame(rep(pHmetrics$Site,2),c(rep('Day', length(pHmetrics$Site)),rep('Night', length(pHmetrics$Site))),
                         c(pHmetrics$pH.day.mean, pHmetrics$pH.night.mean),
                         rep(PData$slopes.TADIC,2),
                         rep(CoverbyGroups$FleshyAlgae+CoverbyGroups$SurfGrass - CoverbyGroups$Inverts,2),
                         c(Temp.day.mean, Temp.night.mean),
                         c(PData$PAR.Int,PData$PAR.Int))
colnames(DataDayNight)<-c('Site','DN','pH','TADIC','b','Temp', 'PAR')

# scale the day and night data
scaledDataDN<-data.frame(scale(DataDayNight$b, center=TRUE, scale = TRUE),scale(DataDayNight$Temp, center=TRUE, scale = TRUE),
                         scale(log(DataDayNight$PAR), center=TRUE, scale = TRUE), scale(DataDayNight$pH, center=TRUE, scale = TRUE),
                         DataDayNight$Site, DataDayNight$DN,
                         rep(PCA$scores[,1],2),
                         rep(PCA$scores[,2],2),
                         scale(DataDayNight$TADIC, center=TRUE, scale = TRUE))


colnames(scaledDataDN)<-c('b','Temp','PAR','pH', 'Site','DN','PC1','PC2','TADIC')

#remove the NA
scaledDataDN<-scaledDataDN[no,]


##-------------------------------Dealing with Heteroscedastic data between day and night--------------------------------------
#range model
lmeData<-data.frame(scaledData,y,PData$Site,PCA$scores[,1],PCA$scores[,2])
lmeData<-lmeData[no,]
colnames(lmeData)[5:8]<-c('pH','Site','PC1','PC2')

FullModel<-lme(pH ~ b + Temp + PAR + PC1+ PC2 +TADIC, random = ~1|Site, data = lmeData) # include TADIC as a predictor
CI<-intervals(FullModel, level = 0.95, which = 'fixed')

# plot the results in an effects plot
par(mfrow=c(1,1), pty='m', mar= c(5.1,18.1,4.1,2.1))

plot(0,type='n', xlim=c(-0.7,1), ylim=c(1,8.2),xlab='Standardized Effect Size', ylab='',
     pch=19, yaxt='n', cex.lab = 1.5, cex.axis = 1.5)
abline(v=0)
axis(2, at=c(1,2,3,4,5,6,7,8), labels=FALSE, tck = -0.01)
text(y = c(1,2,3,4,5,6,7,8), par("usr")[1],labels=c("Temperature",expression(integral('PAR')),"PC 2 (Pool Location)","PC 1 (Pool Size)",
                                                    "TA/DIC Slopes x Day Night","TA/DIC Slopes", '% Cover X Day Night', '% Cover'), srt = 0, pos = 2, xpd = TRUE, cex = 1.5)

arrows(CI$fixed[,1][c(3,4,6,5,7,2)],c(1:4,6,8),CI$fixed[,3][c(3,4,6,5,7,2)],c(1:4,6,8), col = 'black',angle = 90, code = 3, length = 0, lty = 1)
points(CI$fixed[,2][c(3,4,6,5,7,2)],c(1:4,6,8), cex=2, pch=19, col = 'black')


#mean
FullModel<-lme(pH~b*DN + b + Temp + PAR + PC1+ PC2 + TADIC*DN, random=  ~1|Site/DN, data=scaledDataDN, weights=varIdent(form=~1|DN))

# calculate 95% CI from lme model
CI<-intervals(FullModel, level = 0.95, which = 'fixed')

arrows(CI$fixed[,1][c(4,5,7,6,10,8,9,2)],1:8+0.2,CI$fixed[,3][c(4,5,7,6,10,8,9,2)],1:8+0.2, col = 'grey',angle = 90, code = 3, length = 0, lty = 1)
points(CI$fixed[,2][c(4,5,7,6,10,8,9,2)],1:8+0.2, cex=2, pch=19, col = 'grey')

legend('bottomright', c('pH Mean', 'pH Range'), col = c('grey','black'), pch = 19, bty = 'n', cex = 1.5)

dev.off()

#checking model assumptions for mean pH
pdf('../Output/AssumptionsCheck.pdf', height = 10, width = 10)
par(mfrow=c(4,3))
#before accounting for heteroscedasticity
FullModel1<-lme(pH~b*DN + Temp + PAR + PC1+ PC2 +TADIC*DN, random=  ~1|Site/DN, data=scaledDataDN)

E2<-resid(FullModel1, type = 'pearson')
F2<-fitted(FullModel1)
plot(F2,E2, xlab = 'Fitted values', ylab ='Pearson Residuals', col = scaledDataDN$DN)
abline(h=0)
legend('topleft', legend = c('Day','Night'), col = c('black','red'), pch = 16, bty='n')
boxplot(E2~DN, data=scaledDataDN,ylab = 'Residuals')
mtext('Not accounting for unequal variances between Day/Night (Mean Model)', outer = FALSE, line = 2)
qqnorm(E2)
qqline(E2)

#after accounting for unequal variance
FullModel<-lme(pH~b*DN + Temp + PAR + PC1+ PC2 +TADIC*DN, random=  ~1|Site/DN, data=scaledDataDN, weights=varIdent(form=~1|DN))

E2<-resid(FullModel, type = 'pearson')
F2<-fitted(FullModel)
plot(F2,E2, xlab = 'Fitted values', ylab ='Pearson Residuals', col = scaledDataDN$DN)
abline(h=0)
boxplot(E2~DN, data=scaledDataDN)
mtext('Accounting for unequal variances between Day/Night (Mean Model)', outer = FALSE, line = 2)
qqnorm(E2)
qqline(E2)

# now check for equal variance across site
plot(F2,E2, xlab = 'Fitted values', ylab ='Pearson Residuals', col = scaledDataDN$Site)
legend('topleft', legend = unique(scaledDataDN$Site), col = unique(scaledDataDN$Site), pch = 16, bty = 'n')
abline(h=0)
boxplot(E2~Site, data=scaledDataDN,ylab = 'Residuals')
mtext('Checking for equal variance across sites (Mean Model)', outer = FALSE, line = 2)
qqnorm(E2)
qqline(E2)

### assumptions check for range
FullModel1<-lme(pH ~ b + Temp + PAR + PC1+ PC2 +TADIC, random = ~1|Site, data = lmeData)
E2<-resid(FullModel1, type = 'pearson')
F2<-fitted(FullModel1)
plot(F2,E2, xlab = 'Fitted values', ylab ='Pearson Residuals', col = lmeData$Site)
abline(h=0)
boxplot(E2~Site, data=lmeData,ylab = 'Residuals')
mtext('Checking for equal variance across sites (Range Model)', outer = FALSE, line = 2)
qqnorm(E2)
qqline(E2)
# looks good no need to do anything

dev.off()

### visualization of significant factors in mixed effects model for mean and range pH
## Figure 3 ###
pdf(file = '../Output/MMESigSlopes.pdf', height = 11, width = 9, useDingbats=FALSE)
par(mfrow=c(3,2), pty = "s", mai=c(.35,0.01,0.01,0.01))

#pH----------------------------
#make transparent greys
c<-col2rgb("grey")
myGrey<-rgb(c[1],c[2],c[3], alpha=150,  maxColorValue = 255)

c<-col2rgb("grey35")
myGrey35<-rgb(c[1],c[2],c[3], alpha=150,  maxColorValue = 255)

# make the day night pH data long format so that I can look at the day/night interaction for mean pH
pHDayNight<-data.frame(rep(pHmetrics$Site,2),c(rep('Day', length(pHmetrics$Site)),rep('Night', length(pHmetrics$Site))),c(pHmetrics$pH.day.mean, pHmetrics$pH.night.mean),
                       rep(PData$slopes.TADIC,2),rep(CoverbyGroups$FleshyAlgae+CoverbyGroups$SurfGrass - CoverbyGroups$Inverts,2))
colnames(pHDayNight)<-c('Site','DN','pH','TADIC','b')

#slopes versus mean pH with a random effect of site
pHDayNight<-pHDayNight[!is.na(pHDayNight$TADIC),] # remove the missing data

x<-pHDayNight$TADIC 
y<-pHDayNight$pH

# plot mean pH model results
plot(x[pHDayNight$DN=='Day'],y[pHDayNight$DN=='Day'], pch=19, ylab='Mean pH', xlab='', cex=1.5, cex.axis=2, cex.lab=2,
     xlim=c(min(x, na.rm=T),max(x, na.rm=T)),ylim=c(min(y, na.rm=T),max(y, na.rm=T)), col='grey', xaxt ='n')
points(x[pHDayNight$DN=='Night'],y[pHDayNight$DN=='Night'], pch=19,  cex=1.5,col='grey35')

axis(side=1, labels = FALSE) #add back tick marks

#Need to account for variance structure between day and night
model3<-lme(pH~TADIC*DN , random=  ~1|Site/DN, data=pHDayNight, weights=varIdent(form=~1|DN))

# calculate 95% CI for the interaction terms
ef<-effect(term = "TADIC:DN", mod = model3, xlevels = list(TADIC = seq(0.3,0.9,.01)))
ef2<-as.data.frame(ef)

# put color vectors together
gcols<-c(myGrey, myGrey35)

for (i in 1:2){
  polygon(c(ef2$TADIC[ef2$DN==DN[i]], rev(ef2$TADIC[ef2$DN==DN[i]])),
          c(ef2$upper[ef2$DN==DN[i]], rev(ef2$lower[ef2$DN==DN[i]])),
          col = gcols[i], border = NA)
  lines(ef2$TADIC[ef2$DN==DN[i]],ef2$fit[ef2$DN==DN[i]], col='black', lwd=2)
}
legend('topright',c('Day','Night'), col = c('grey','grey35'), pch=19, bty = 'n', cex = 1.5)
r.squaredGLMM(model3) 

# next community versus pH----------------------
#mean pH
x<-pHDayNight$b 
y<-pHDayNight$pH

plot(x[pHDayNight$DN=='Day'],y[pHDayNight$DN=='Day'], pch=19, xlab='', ylab='', cex=1.5, cex.axis=2, cex.lab=2,
     xlim=c(min(x, na.rm=T),max(x, na.rm=T)), ylim=c(min(y, na.rm=T),max(y, na.rm=T)), col="grey", axes=FALSE)
box()

axis(side=1, labels = FALSE) #add back tick marks
axis(side=2, labels = FALSE) #add back tick marks
abline(v=0, lty=2)
points(x[pHDayNight$DN=='Night'],y[pHDayNight$DN=='Night'], pch=19,  cex=1.5,col='grey35')

model3<-lme(pH~b*DN , random=  ~1|Site/DN, data=pHDayNight, weights=varIdent(form=~1|DN))

# calculate 95% CI for the interaction terms
ef<-effect(term = "b:DN", mod = model3, xlevels = list(b = seq(min(pHDayNight$b),max(pHDayNight$b),.1)))
ef2<-as.data.frame(ef)

for (i in 1:2){
  polygon(c(ef2$b[ef2$DN==DN[i]], rev(ef2$b[ef2$DN==DN[i]])),
          c(ef2$upper[ef2$DN==DN[i]], rev(ef2$lower[ef2$DN==DN[i]])),
          col = gcols[i], border = NA)
  lines(ef2$b[ef2$DN==DN[i]],ef2$fit[ef2$DN==DN[i]], col='black', lwd=2)
}


r.squaredGLMM(model3) 

#range pH
#slopes versus range pH
# for the effects function to work, everything needs to be in a clean dataframe
NewD<-data.frame(PData$slopes.TADIC[!is.na(PData$slopes.TADIC)],pHmetrics$rangepH[!is.na(PData$slopes.TADIC)],
                 pHmetrics$Site[!is.na(PData$slopes.TADIC)], b[!is.na(PData$slopes.TADIC)], pHmetrics$maxdiffpH[!is.na(PData$slopes.TADIC)])

colnames(NewD)<-c('TADIC','pHRange','Site','b','maxdiffpH')
x<-PData$slopes.TADIC[!is.na(PData$slopes.TADIC)] #remove the missing data
y<-pHmetrics$rangepH[!is.na(PData$slopes.TADIC)]
plot(x,y, pch=19, ylab='Range pH', xlab='TA vs DIC slopes', cex=1.5, cex.axis=2, cex.lab=2,
     xlim=c(min(x, na.rm=T),max(x, na.rm=T)), col='grey')

model3<-lme(pHRange~TADIC, random=~1|Site, data=NewD, na.action = na.exclude)

# calculate 95% CI for the interaction terms
ef<-effect(term = "TADIC", mod = model3, xlevels = list(TADIC = seq(min(NewD$TADIC),max(NewD$TADIC),.01)))
ef2<-as.data.frame(ef)

polygon(c(ef2$TADIC, rev(ef2$TADIC)),
        c(ef2$upper, rev(ef2$lower)),
        col = myGrey, border = NA)
lines(ef2$TADIC,ef2$fit, col='black', lwd=2)


#community versus range
y<- pHmetrics$rangepH[!is.na(PData$slopes.TADIC)] #remove the missing data
x<- CoverbyGroups$FleshyAlgae[!is.na(PData$slopes.TADIC)]+CoverbyGroups$SurfGrass[!is.na(PData$slopes.TADIC)] - 
  CoverbyGroups$Inverts[!is.na(PData$slopes.TADIC)]

plot(x,y, pch=19, xlab='', ylab='', cex=1.5, cex.axis=2, cex.lab=2,
     xlim=c(min(x, na.rm=T),max(x, na.rm=T)), col="grey", yaxt='n')

axis(side=2, labels = FALSE) #add back tick marks
abline(v=0, lty=2)

model3<-lme(pHRange~b, random = ~1|Site, data=NewD, na.action = na.exclude)

# calculate 95% CI for the interaction terms
ef<-effect(term = "b", mod = model3, xlevels = list(b = seq(min(NewD$b),max(NewD$b),.1)))
ef2<-as.data.frame(ef)

polygon(c(ef2$b, rev(ef2$b)),
        c(ef2$upper, rev(ef2$lower)),
        col = myGrey, border = NA)
lines(ef2$b,ef2$fit, col='black', lwd=2)

title(xlab = expression(paste(Delta, ' Benthic Cover')), cex.lab = 2, line = 4)
title(xlab = 'Producers - Invertebrates', cex.lab = 2, line = 2.5)


#next cover versus slopes-----------------
y<-PData$slopes.TADIC[!is.na(PData$slopes.TADIC)] #remove the missing data
x<- CoverbyGroups$FleshyAlgae[!is.na(PData$slopes.TADIC)]+CoverbyGroups$SurfGrass[!is.na(PData$slopes.TADIC)] - 
  CoverbyGroups$Inverts[!is.na(PData$slopes.TADIC)]

plot(x,y, pch=19, xlab='', ylab='TA vs DIC slopes', cex=1.5, cex.axis=2, cex.lab=2,
     xlim=c(min(x, na.rm=T),max(x, na.rm=T)), col="grey")
abline(v=0, lty=2)

model3<-lme(TADIC~b, random=~1|Site, data=NewD, na.action = na.exclude)

# calculate 95% CI for the interaction terms
ef<-effect(term = "b", mod = model3, xlevels = list(b = seq(min(NewD$b),max(NewD$b),.1)))
ef2<-as.data.frame(ef)

polygon(c(ef2$b, rev(ef2$b)),
        c(ef2$upper, rev(ef2$lower)),
        col = myGrey, border = NA)
lines(ef2$b,ef2$fit, col='black', lwd=2)
title(xlab = expression(paste(Delta, ' Benthic Cover')), cex.lab = 1.5, line = 4)
title(xlab = 'Producers - Invertebrates', cex.lab = 1.5, line = 2.5)

dev.off()


### Community and PCA plot ###
## Figure 4 ###
pdf(file = '../Output/CommunityFig.pdf', width = 13, height = 13, useDingbats=FALSE)

#colors for sites
cols.all<-c(rep('magenta',15),rep('red',15),rep('lightblue',15),rep('#0606f7',15))

par(mar=c(7,4.5,4.1,1), bg=NA)
layout(matrix(c(1,1,2,2,3,4,5,5,6,6), nrow = 2, ncol = 5, byrow = TRUE))

#NMDS
ord<-metaMDS(PercentTotal,k=2, distance='bray', trymax = 30) #create the ordination output using bray curtis
plot(1, type='n', xlim=c(-2,2), ylim=c(-1.5,1.5), xlab='nMDS1', ylab='nMDS2',cex.axis =2, cex.lab=2,
     xaxt='n', yaxt='n', bty='n')
box(which = "plot", lty = "solid")
axis(side=2, cex.axis=2)
axis(1,  cex.axis=2)


LightBlueT=rgb(173,216,230,max=255,alpha=0.75) 
MagentaT=rgb(255,0,255,max=255,alpha=0.75) 
RedT=rgb(255,0,0,max=255,alpha=0.75) 
DBlueT=rgb(6,6,247,max=255,alpha=0.75) 

# CI around the points in NMDS
ordiellipse(ord, groups=Sessile$Site, label=F, kind='ehull', border='white', col=c(RedT,MagentaT,LightBlueT,DBlueT ), lwd=2, draw ='polygon')

# plot the points
points(ord$points[Sessile$Site=='Bodega',1],ord$points[Sessile$Site=='Bodega',2], pch=19, col='magenta', cex=2)
points(ord$points[Sessile$Site=='BobCreek',1],ord$points[Sessile$Site=='BobCreek',2], pch=19, col='red', cex=2)
points(ord$points[Sessile$Site=='CDM',1],ord$points[Sessile$Site=='CDM',2], pch=19, col='lightblue', cex=2)
points(ord$points[Sessile$Site=='Monterey',1],ord$points[Sessile$Site=='Monterey',2], pch=19, col='#0606f7', cex=2)
legend('topright',legend=c('Bob Creek','Bodega Bay','Monterey Bay','Corona del Mar'),
       col=c('red','magenta','#0606f7','lightblue'), pch=19, bty='n', y.intersp = 1, x.intersp = 0.5, cex = 2)

## plot percent of pools by group for each site
##remove the site and square info from the groups and join them into one vector so its the same length as percent total

bp<-barplot(t(Cover.mean)[c(1:5,8,9),c(3,4,2,1)], # remove the fish because they are barely there 
            col=c('black','darkgreen','pink','brown','lightgreen','gray38','grey'),
            ylab='',cex.axis =2,cex.names =2 , cex.lab=2, xaxt='n',
            xaxt='n', yaxt='n', bty='n', ylim=c(0,105))
title(ylab = 'Precent Cover', line = 2.75, cex.lab=2)
axis(1, bp, labels = FALSE, cex.axis=2)
axis(2, c(0,50,100), cex.axis=2)
text(bp-.2, par("usr")[3] - 8,  labels = c('Corona del Mar','Monterey Bay','Bodega Bay','Bob Creek'), srt = 45, cex=2, pos = 1, xpd = TRUE)
# add the legend
par(mar=c(5.1,0,4.1,1))
plot(0, 0, type = "n", ann = F, axes = F)
legend('left', legend = c('Invertebrates','Mussels Only','Surf Grass','Non-coralline Crust', 'Coralline Crust','Fleshy Algae','Bare Rock'),
       bty='n', pch=15, col=c('grey','grey38','lightgreen','brown','pink','darkgreen','black'),cex=2)


#PCA for physical variables
x<-data.frame(scale(PData[,c('Perimeter','MaxDepth','SA','Vol','TidePercent','SAV')], scale= apply(PData[,c('Perimeter','MaxDepth','SA','Vol','TidePercent','SAV')], 2, sd, na.rm = TRUE), 
                    center = TRUE))
PCA<- princomp( na.omit(x ))

#add back the NA so that I can use the PC scores in tests
PCA$scores<-rbind(PCA$scores[1:48,], rep(NA,ncol(PCA$scores)),PCA$scores[c(49:52),], rep(NA,ncol(PCA$scores)), PCA$scores[c(53:58),])

# to get the variance explained by each axis
s<-matrix(summary(PCA))

#create a biplot
par(mar=c(7,4.5,4.1,1), bg=NA)
plot(0, 0, type = "n", ann = F, axes = F)
plot(PCA$scores[,1]/10, PCA$scores[,2]/10, col=cols.all, xlim=c(-0.5,0.4), ylim=c(-0.4,0.4), pch=19, axes=F, xlab="", ylab="", cex=1.5, cex.axis = 2)
par(new=TRUE)
biplot(PCA, xlabs=rep("", dim(PCA$scores)[1]), cex=2, xlab = 'Comp.1 (55%)', ylab='Comp.2 (19%)', col = 'black', cex.lab = 2, cex.axis = 2)
dev.off() 


##### NCP vs pH vs NCP models ####################
## Figure 5 ###
pdf(file = '../Output/NECNCP.pdf', height = 10, width = 15, useDingbats=FALSE)
par(mfrow=c(1,2), pty = 's',mar=c(5.1,6.1,4.1,2.1))
# put all data into a data.frame
MetabModelData<-data.frame(delta$NEC,MeanData$pH, MeanData$Temp.pool,MeanData$PoolID,MeanData$Site, delta$NCP)

# make the sites a specific color
color_pallete_function <- colorRampPalette(
  colors = c('lightblue','magenta','#0606f7','red'),
  space = "Lab" # Option used when colors do not represent a quantitative scale
)


scols<-c('#0606f7','red','lightblue','magenta')
#transparent cols
c<-col2rgb("red")
myred<-rgb(c[1],c[2],c[3], alpha=50,  maxColorValue = 255)
c<-col2rgb("lightblue")
myblue<-rgb(c[1],c[2],c[3], alpha=50,  maxColorValue = 255)
c<-col2rgb("#0606f7")
mydblue<-rgb(c[1],c[2],c[3], alpha=50,  maxColorValue = 255)
c<-col2rgb("magenta")
mymagenta<-rgb(c[1],c[2],c[3], alpha=50,  maxColorValue = 255)


tcols<-c(mydblue, myred, myblue,mymagenta)

# remove missing data
MetabModelData<-MetabModelData[!is.na(delta$NEC),]
colnames(MetabModelData)<-c('NEC','pH','Temp','PoolID','Site','NCP')

##### pH versus NCP ###############
# model of NEC as a function of pH with a site interaction and random effects of site and pool
NCPModel<-lme(pH ~ NCP*Site, random = ~1|PoolID ,na.action = na.omit, data = MetabModelData)

# put colors in the right order
num_colors <- nlevels(MetabModelData$Site)
site_colors <- color_pallete_function(num_colors)

# plot the results
plot(MetabModelData$NCP,MetabModelData$pH, col = site_colors[MetabModelData$Site],
     ylab = 'pH', xlab ='' , pch = 19, cex.lab = 2, cex.axis = 2)

title(xlab = expression(paste('NCP (mmol m'^{-2}, ' hr'^{-1},')')), cex.lab = 2, line = 4)
# calulate prediction intervals
ef<-effect(term = "NCP:Site", mod = NCPModel, xlevels = list(NCP = seq(min(MetabModelData$NCP),max(MetabModelData$NCP),.1)))
ef2<-as.data.frame(ef)
# change the order of the levels
ef2$Site<-factor(ef2$Site, levels=c("CDM","Monterey","Bodega","BobCreek"))

# plot CI 
for (i in c(3,4,1,2)){
  polygon(c(ef2$NCP[ef2$Site==Sites[i]], rev(ef2$NCP[ef2$Site==Sites[i]])),
          c(ef2$upper[ef2$Site==Sites[i]], rev(ef2$lower[ef2$Site==Sites[i]])),
          col = tcols[i], border = NA)
  lines(ef2$NCP[ef2$Site==Sites[i]],ef2$fit[ef2$Site==Sites[i]], col=scols[i], lwd=2)
}
abline(v=0, lty=2)

legend('bottomright',legend =c("Bob Creek","Bodega","Monterey","Corona del Mar"), col = c('red','#0606f7','magenta','lightblue'), pch = 19, bty = 'n' )
r.squaredGLMM(NCPModel) 
##### NEC versus pH ###############
# model of NEC as a function of pH with a site interaction and random effects of site and pool
NECModel<-lme(NEC ~ pH*Site, random = ~1|PoolID ,na.action = na.omit, data = MetabModelData)

# plot the results
plot(MetabModelData$pH,MetabModelData$NEC, col = site_colors[MetabModelData$Site],
     xlab = '', ylab = expression(paste('NEC (mmol m'^{-2}, ' hr'^{-1},')')), pch = 19, cex.lab = 2, cex.axis = 2)
title(xlab = 'pH', cex.lab = 2, line = 4)
# calulate prediction intervals
ef<-effect(term = "pH:Site", mod = NECModel, xlevels = list(pH = seq(min(MetabModelData$pH),max(MetabModelData$pH),.01)))
ef2<-as.data.frame(ef)
# change the order of the levels
ef2$Site<-factor(ef2$Site, levels=c("CDM","Monterey","Bodega","BobCreek"))

for (i in c(3,4,1,2)){
  polygon(c(ef2$pH[ef2$Site==Sites[i]], rev(ef2$pH[ef2$Site==Sites[i]])),
          c(ef2$upper[ef2$Site==Sites[i]], rev(ef2$lower[ef2$Site==Sites[i]])),
          col = tcols[i], border = NA)
  lines(ef2$pH[ef2$Site==Sites[i]],ef2$fit[ef2$Site==Sites[i]], col=scols[i], lwd=2)
}
r.squaredGLMM(NECModel) 
abline(h=0, lty=2)

### NCP versus NEC
NCPNECModel<-lme(NEC ~ NCP*Site, random = ~1|PoolID ,na.action = na.omit, data = MetabModelData)

# plot the results
plot(MetabModelData$NCP,MetabModelData$NEC, col = site_colors[MetabModelData$Site],
     ylab = '', xlab ='' , pch = 19, cex.lab = 2, cex.axis = 2)

title(ylab = expression(paste('NEC (mmol m'^{-2}, ' hr'^{-1},')')), cex.lab = 2, line = 4)
title(xlab = expression(paste('NCP (mmol m'^{-2}, ' hr'^{-1},')')), cex.lab = 2, line = 4)
# calulate prediction intervals
ef<-effect(term = "NCP:Site", mod = NCPNECModel, xlevels = list(NCP = seq(min(MetabModelData$NCP, na.rm=T),max(MetabModelData$NCP, na.rm=T),.1)))
ef2<-as.data.frame(ef)
# change the order of the levels
ef2$Site<-factor(ef2$Site, levels=c("CDM","Monterey","Bodega","BobCreek"))

for (i in c(3,4,1,2)){
  polygon(c(ef2$NCP[ef2$Site==Sites[i]], rev(ef2$NCP[ef2$Site==Sites[i]])),
          c(ef2$upper[ef2$Site==Sites[i]], rev(ef2$lower[ef2$Site==Sites[i]])),
          col = tcols[i], border = NA)
  lines(ef2$NCP[ef2$Site==Sites[i]],ef2$fit[ef2$Site==Sites[i]], col=scols[i], lwd=2)
}
abline(v=0, lty=2)
abline(h=0, lty=2)
legend('bottomright',legend =c("Bob Creek","Bodega","Monterey","Corona del Mar"), col = c('red','#0606f7','magenta','lightblue'), pch = 19, bty = 'n' )
r.squaredGLMM(NCPModel) 
#

dev.off()

## community versus the max difference between tide pools and the ocean for supplement
#community versus range
pdf(file = '../Output/MaxDiff.pdf', height = 7.5, width = 7, useDingbats=FALSE)

par(mar=c(5.1,7.1,4.1,2.1))
y<- pHmetrics$maxdiffpH[!is.na(PData$slopes.TADIC)] #remove the missing data
x<- CoverbyGroups$FleshyAlgae[!is.na(PData$slopes.TADIC)]+CoverbyGroups$SurfGrass[!is.na(PData$slopes.TADIC)] - 
  CoverbyGroups$Inverts[!is.na(PData$slopes.TADIC)]

plot(x,y, pch=19, xlab='', cex=1.5, cex.axis=1.5, cex.lab=1.5,
     xlim=c(min(x, na.rm=T),max(x, na.rm=T)), col="grey", ylab = 'Max (Tide Pool - Ocean pH)')


#axis(side=2, labels = FALSE) #add back tick marks
abline(v=0, lty=2)
abline(h=0, lty=2)

model3<-lme(maxdiffpH~b, random = ~1|Site, data=NewD, na.action = na.exclude)

# calculate 95% CI for the interaction terms
ef<-effect(term = "b", mod = model3, xlevels = list(b = seq(min(NewD$b),max(NewD$b),.1)))
ef2<-as.data.frame(ef)

polygon(c(ef2$b, rev(ef2$b)),
        c(ef2$upper, rev(ef2$lower)),
        col = myGrey, border = NA)
lines(ef2$b,ef2$fit, col='black', lwd=2)

title(xlab = 'Producer Dominance', cex.lab = 1.5, line = 2.5)
title(xlab = '% Producers - % Consumers', cex.lab = 1.5, line = 4)
dev.off()

#chemistry ranges by day and night

