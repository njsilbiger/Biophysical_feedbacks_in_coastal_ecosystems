## script to read in the temperature data  from the different MASCOT sites
## created Nyssa Silbiger 9/13/2016


# set working directory
dname<-'../Data/RawData'

## read in the data ----------

# first source the calibration code data which is based off of common garden light data from Monterey
# This code calculates the average deviation from the mean of each sensor during a common garden
source('LightCalibrationCode.R')
# The output here is a slope and intercept for each pool.  Add this to the light data for every site

#Bodega Temp-----------------------------------------------------
#create an empty dataframe to read in all the data
BodegaTemp<-data.frame(matrix(, nrow = 4618, ncol = 16))

#files to read in
files<-paste(rep('TP',15),1:15,sep='')

#for loop looping through all the temperature files
for (i in 1:15){
d<-read.csv(paste(dname,'/BodegaTemp/',files[i],'.csv',sep=''),
         stringsAsFactors=FALSE,skip=1, col.names = c('NA','DateTime','Temp','NA','NA','NA','NA','NA'))[ ,2:3]

#convert the dates to datetime
d$DateTime<-as.POSIXct(d$DateTime, format="%m/%d/%y %I:%M:%S %p")

#round the times to the nearest minute
d$DateTime<-format(round(d$DateTime, units="mins"), format="%Y-%m-%d %H:%M")

#find the start and stop time of all the dates
start<- as.POSIXct('2016-7-25 00:01', format="%Y-%m-%d %H:%M")
end<- as.POSIXct('2016-7-28 05:00', format="%Y-%m-%d %H:%M")

#only pull out the data that I need during the actual sampling period

d<-d[d$DateTime>start & d$DateTime<end,]
BodegaTemp[,1]<-d$DateTime
BodegaTemp[,i+1]<-d$Temp
  
}
#change the column names so that it is labeled by tidepool
colnames(BodegaTemp)<-c('DateTime',files)
##for some reason the dates are showing up as characters... this converts back to datetime
BodegaTemp$DateTime<-strptime(BodegaTemp$DateTime, '%Y-%m-%d %H:%M')

#plot all the data on top of each other
plot(BodegaTemp$DateTime, BodegaTemp$TP1, type = "l")
for (i in 1:15){
  lines(BodegaTemp$DateTime, BodegaTemp[,i+1], col=i)
}

legend('topleft',legend = files,col = 1:15, lty=1, bty = 'n')

#Bodega light-----------------------------------------------------
#create an empty dataframe to read in all the data
BodegaLight<-data.frame(matrix(, nrow = 4618, ncol = 16))

#files to read in
files<-paste(rep('TP',15),1:15,'_light',sep='')

#for loop looping through all the temperature files
for (i in 1:15){
  d<-read.csv(paste(dname,'/BodegaTemp/',files[i],'.csv',sep=''),
              stringsAsFactors=FALSE,skip=1, col.names = c('NA','DateTime','Temp','Light','NA','NA','NA','NA','NA'))[ ,2:4]
  
  #convert the dates to datetime
  d$DateTime<-as.POSIXct(d$DateTime, format="%m/%d/%y %I:%M:%S %p")
  
  #round the times to the nearest minute
  d$DateTime<-format(round(d$DateTime, units="mins"), format="%Y-%m-%d %H:%M")
  
  #only pull out the data that I need
  
  # filter a moving average for every 30 minutes
  d$Light2<-filter(d$Light, sides=2, rep(1/30,30))  
  
  d<-d[d$DateTime>start & d$DateTime<end,]
  BodegaLight[,1]<-d$DateTime
  BodegaLight[,i+1]<-d$Light2
  
}
#change the column names so that it is labeled by tidepool
colnames(BodegaLight)<-c('DateTime',files)
##for some reason the dates are showing up as characters... this converts back to datetime
BodegaLight$DateTime<-strptime(BodegaLight$DateTime, '%Y-%m-%d %H:%M')

# add the correction factor (slope and int are output from the calibration code)
for (i in 2:14){
  BodegaLight[,i]<-BodegaLight[,i]+ (BodegaLight[,i]*(slope[i]))  
}



# convert the light data from lumes/ft2 to m2
BodegaLight[,2:16]<-BodegaLight[,2:16]*0.092903

#Now convert to PAR based on Long 2012
#parameters from Long 2012
#A1 = -8165.9
#t1 = 1776.4
#y0 =  8398.2

A1 = -4924
t1 = 20992.9
y0 =  4929.0


BodegaPAR<-BodegaLight
BodegaPAR[,2:16] = A1*exp(-BodegaLight[,2:16]/t1) + y0

#plot all the data on top of each other
plot(BodegaPAR$DateTime, BodegaPAR$TP1_light, type = "l", ylim=c(0,max(BodegaPAR[,2:15])))
for (i in 1:15){
  lines(BodegaPAR$DateTime, BodegaPAR[,i+1], col=i)
}

legend('topleft',legend = files,col = 1:15, lty=1, bty = 'n')

##Bob Creek------------------
#Temp-----------------------------------------------------
#create an empty dataframe to read in all the data
BobCreekTemp<-data.frame(matrix(, nrow = 4538, ncol = 16))

#files to read in
files<-paste(rep('TP',15),1:15,sep='')

#for loop looping through all the temperature files
for (i in 1:15){
  d<-read.csv(paste(dname,'/BobCreekTemp/',files[i],'.csv',sep=''),
              stringsAsFactors=FALSE,skip=1, col.names = c('NA','DateTime','Temp','NA','NA','NA','NA','NA'))[ ,2:3]
  
  #convert the dates to datetime
  d$DateTime<-as.POSIXct(d$DateTime, format="%m/%d/%y %I:%M:%S %p")
  #round the times to the nearest minute
  d$DateTime<-format(round(d$DateTime, units="mins"), format="%Y-%m-%d %H:%M")
  
  #find the start and stop time of all the dates
  start<- as.POSIXct('2016-8-8 00:01', format="%Y-%m-%d %H:%M")
  end<- as.POSIXct('2016-8-11 03:40', format="%Y-%m-%d %H:%M")
  
  #only pull out the data that I need
  d<-d[d$DateTime>start & d$DateTime<end,]
  BobCreekTemp[,1]<-d$DateTime
  BobCreekTemp[,i+1]<-d$Temp
  
}
#change the column names so that it is labeled by tidepool
colnames(BobCreekTemp)<-c('DateTime',files)
##for some reason the dates are showing up as characters... this converts back to datetime
BobCreekTemp$DateTime<-strptime(BobCreekTemp$DateTime, '%Y-%m-%d %H:%M')


#plot all the data on top of each other
plot(BobCreekTemp$DateTime, BobCreekTemp$TP1, type = "l", ylim=c(8,30))
for (i in 1:15){
  lines(BobCreekTemp$DateTime, BobCreekTemp[,i+1], col=i)
}

legend('topleft',legend = files,col = 1:15, lty=1, bty = 'n')

#BobCreek light-----------------------------------------------------
#create an empty dataframe to read in all the data
BobCreekLight<-data.frame(matrix(, nrow = 4538, ncol = 16))

#files to read in
files<-paste(rep('TP',15),1:15,'_light',sep='')

#for loop looping through all the temperature files
for (i in 1:15){
  d<-read.csv(paste(dname,'/BobCreekTemp/',files[i],'.csv',sep=''),
              stringsAsFactors=FALSE,skip=1, col.names = c('NA','DateTime','Temp','Light','NA','NA','NA','NA','NA'))[ ,2:4]
  
  #convert the dates to datetime
  d$DateTime<-as.POSIXct(d$DateTime, format="%m/%d/%y %I:%M:%S %p")
  
  #round the times to the nearest minute
  d$DateTime<-format(round(d$DateTime, units="mins"), format="%Y-%m-%d %H:%M")
  
  # filter a moving average for every 30 minutes
  d$Light<-filter(d$Light, sides=2, rep(1/30,30))  
  
  #only pull out the data that I need
  d<-d[d$DateTime>start & d$DateTime<end,]
  BobCreekLight[,1]<-d$DateTime
  BobCreekLight[,i+1]<-d$Light
  
}
#change the column names so that it is labeled by tidepool
colnames(BobCreekLight)<-c('DateTime',files)
##for some reason the dates are showing up as characters... this converts back to datetime
BobCreekLight$DateTime<-strptime(BobCreekLight$DateTime, '%Y-%m-%d %H:%M')

# add the correction factor (slope and int are output from the calibration code)
for (i in 2:14){
  BobCreekLight[,i]<-BobCreekLight[,i]+ (BobCreekLight[,i]*(slope[i])+int[i])  
}


# convert the light data from lumes/ft2 to m2
BobCreekLight[,2:16]<-BobCreekLight[,2:16]*0.092903

#Now convert to PAR based on Long 2012
BobCreekPAR<-BobCreekLight
BobCreekPAR[,2:16] = A1*exp(-BobCreekLight[,2:16]/t1) + y0


#plot all the data on top of each other
plot(BobCreekPAR$DateTime, BobCreekLight$TP1_light, type = "l", ylim=c(0,max(BobCreekPAR[,2:15])))
for (i in 1:15){
  lines(BobCreekPAR$DateTime, BobCreekPAR[,i+1], col=i)
}

legend('topleft',legend = files,col = 1:15, lty=1, bty = 'n')

## Monterey Data--------------------
#Temp-----------------------------------------------------
#create an empty dataframe to read in all the data
MontereyTemp<-data.frame(matrix(, nrow = 4658, ncol = 16))

#files to read in
files<-paste(rep('TP',15),1:15,sep='')

#for loop looping through all the temperature files
for (i in 1:15){
  d<-read.csv(paste(dname,'/MontereyTemp/',files[i],'.csv',sep=''),
              stringsAsFactors=FALSE,skip=1, col.names = c('NA','DateTime','Temp','NA','NA','NA','NA','NA'))[ ,2:3]
  
  #convert the dates to datetime
  d$DateTime<-as.POSIXct(d$DateTime, format="%m/%d/%y %I:%M:%S %p")
  #round the times to the nearest minute
  d$DateTime<-format(round(d$DateTime, units="mins"), format="%Y-%m-%d %H:%M")
  
  #find the start and stop time of all the dates
  start<- as.POSIXct('2016-7-10 00:01', format="%Y-%m-%d %H:%M")
  end<- as.POSIXct('2016-7-13 05:40', format="%Y-%m-%d %H:%M")
  
  #only pull out the data that I need
  d<-d[d$DateTime>start & d$DateTime<end,]
  MontereyTemp[,1]<-d$DateTime
  MontereyTemp[,i+1]<-d$Temp
  
}
#change the column names so that it is labeled by tidepool
colnames(MontereyTemp)<-c('DateTime',files)
##for some reason the dates are showing up as characters... this converts back to datetime
MontereyTemp$DateTime<-strptime(MontereyTemp$DateTime, '%Y-%m-%d %H:%M')


#plot all the data on top of each other
plot(MontereyTemp$DateTime, MontereyTemp$TP1, type = "l", ylim=c(11,30))
for (i in 1:15){
  lines(MontereyTemp$DateTime, MontereyTemp[,i+1], col=i)
}

legend('topleft',legend = files,col = 1:15, lty=1, bty = 'n')

#Monterey light-----------------------------------------------------
#create an empty dataframe to read in all the data
MontereyLight<-data.frame(matrix(, nrow = 4658, ncol = 16))

#files to read in
files<-paste(rep('TP',15),1:15,'_light',sep='')

#for loop looping through all the temperature files
for (i in 1:15){
  d<-read.csv(paste(dname,'/MontereyTemp/',files[i],'.csv',sep=''),
              stringsAsFactors=FALSE,skip=1, col.names = c('NA','DateTime','Temp','Light','NA','NA','NA','NA','NA'))[ ,2:4]
  
  #convert the dates to datetime
  d$DateTime<-as.POSIXct(d$DateTime, format="%m/%d/%y %I:%M:%S %p")
  
  #round the times to the nearest minute
  d$DateTime<-format(round(d$DateTime, units="mins"), format="%Y-%m-%d %H:%M")
  
  # filter a moving average for every 30 minutes
  d$Light<-filter(d$Light, sides=2, rep(1/30,30))  
  
  #only pull out the data that I need
  d<-d[d$DateTime>start & d$DateTime<end,]
  MontereyLight[,1]<-d$DateTime
  MontereyLight[,i+1]<-d$Light
  
}
#change the column names so that it is labeled by tidepool
colnames(MontereyLight)<-c('DateTime',files)
##for some reason the dates are showing up as characters... this converts back to datetime
MontereyLight$DateTime<-strptime(MontereyLight$DateTime, '%Y-%m-%d %H:%M')

# add the correction factor
#MontereyLight[,2:16]<-MontereyLight[,2:16]+AverageDev

for (i in 2:14){
  MontereyLight[,i]<- MontereyLight[,i] + (MontereyLight[,i]*(slope[i])+int[i])  
}


# convert the light data from lumes/ft2 to m2
MontereyLight[,2:16]<-MontereyLight[,2:16]*0.092903

#Now convert to PAR based on Long 2012
MontereyPAR<-MontereyLight
MontereyPAR[,2:16] = A1*exp(-MontereyLight[,2:16]/t1) + y0

#plot all the data on top of each other
plot(MontereyPAR$DateTime, MontereyPAR$TP1_light, type = "l", ylim=c(0,max(MontereyPAR[,2:15])))
for (i in 1:15){
  lines(MontereyPAR$DateTime, MontereyPAR[,i+1], col=i)
}

legend('topleft',legend = files,col = 1:15, lty=1, bty = 'n')

#clear unneeded variables
rm(d)

## CDM Data--------------------
#Temp-----------------------------------------------------
#create an empty dataframe to read in all the data
CDMTemp<-data.frame(matrix(, nrow = 16078, ncol = 14))

#files to read in
files<-paste(rep('TP',13),c(2:4,6:15),sep='')

#for loop looping through all the temperature files
for (i in 1:13){ #two of the sensors got pulled off :(
  d<-read.csv(paste(dname,'/CDMTemp/',files[i],'.csv',sep=''),
              stringsAsFactors=FALSE,skip=1, col.names = c('NA','DateTime','Temp','NA','NA','NA','NA','NA'))[ ,2:3]
  
  #convert the dates to datetime
  d$DateTime<-as.POSIXct(d$DateTime, format="%m/%d/%y %I:%M:%S %p")
  #round the times to the nearest minute
  d$DateTime<-format(round(d$DateTime, units="mins"), format="%Y-%m-%d %H:%M")
  
  #find the start and stop time of all the dates
  start<- as.POSIXct('2016-9-16 00:01', format="%Y-%m-%d %H:%M")
  end<- as.POSIXct('2016-9-27 04:00', format="%Y-%m-%d %H:%M")
  
  #only pull out the data that I need
  d<-d[d$DateTime>start & d$DateTime<end,]
  CDMTemp[,1]<-d$DateTime
  CDMTemp[,i+1]<-d$Temp
  
}
#change the column names so that it is labeled by tidepool
colnames(CDMTemp)<-c('DateTime',files)
##for some reason the dates are showing up as characters... this converts back to datetime
CDMTemp$DateTime<-strptime(CDMTemp$DateTime, '%Y-%m-%d %H:%M')


#plot all the data on top of each other
plot(CDMTemp$DateTime, CDMTemp$TP2, type = "l", ylim=c(15,27))
for (i in 1:13){
  lines(CDMTemp$DateTime, CDMTemp[,i+1], col=i)
}

legend('topleft',legend = files,col = 1:13, lty=1, bty = 'n')

#CDM light-----------------------------------------------------
#create an empty dataframe to read in all the data
CDMLight<-data.frame(matrix(, nrow = 16078, ncol = 14))

#files to read in
files<-paste(rep('TP',13),c(2:4,6:15),'_light',sep='')

#for loop looping through all the temperature files
for (i in 1:13){
  d<-read.csv(paste(dname,'/CDMTemp/',files[i],'.csv',sep=''),
              stringsAsFactors=FALSE,skip=1, col.names = c('NA','DateTime','Temp','Light','NA','NA','NA','NA','NA'))[ ,2:4]
  
  #convert the dates to datetime
  d$DateTime<-as.POSIXct(d$DateTime, format="%m/%d/%y %I:%M:%S %p")
  
  #round the times to the nearest minute
  d$DateTime<-format(round(d$DateTime, units="mins"), format="%Y-%m-%d %H:%M")
  
  # filter a moving average for every 30 minutes
  d$Light<-filter(d$Light, sides=2, rep(1/30,30))  
  
  #only pull out the data that I need
  d<-d[d$DateTime>start & d$DateTime<end,]
  CDMLight[,1]<-d$DateTime
  CDMLight[,i+1]<-d$Light
  
}
#change the column names so that it is labeled by tidepool
colnames(CDMLight)<-c('DateTime',files)
##for some reason the dates are showing up as characters... this converts back to datetime
CDMLight$DateTime<-strptime(CDMLight$DateTime, '%Y-%m-%d %H:%M')

# add the correction factor (slope and int are output from the calibration code)
for (i in 2:13){
CDMLight[,i]<-CDMLight[,i]+ (CDMLight[,i]*(slope[i])+int[i])  
}

# convert the light data from lumes/ft2 to m2
CDMLight[,2:14]<-CDMLight[,2:14]*0.092903

#Now convert to PAR based on Long 2012
CDMPAR<-CDMLight
CDMPAR[,2:14] = A1*exp(-CDMLight[,2:14]/t1) + y0

#plot all the data on top of each other
plot(CDMPAR$DateTime, CDMPAR$TP2_light, type = "l", ylim=c(0,max(CDMPAR[,2:13])))
for (i in 1:13){
  lines(CDMPAR$DateTime, CDMPAR[,i+1], col=i)
}

legend('topleft',legend = files,col = 1:13, lty=1, bty = 'n')

#clear unneeded variables
rm(d)

