## QA/QC the light data

MontereyLight<-data.frame(matrix(, nrow = 10000, ncol = 16))

#files to read in
files<-paste(rep('TP',15),1:15,'_light',sep='')

#directory
dname<-'../Data/RawData'

#for loop looping through all the temperature files
for (i in 1:15){
  d<-read.csv(paste(dname,'/MontereyTemp/',files[i],'.csv',sep=''),
              stringsAsFactors=FALSE,skip=1, col.names = c('NA','DateTime','Temp','Light','NA','NA','NA','NA','NA'))[ ,2:4]
  
  #convert the dates to datetime
  d$DateTime<-as.POSIXct(d$DateTime, format="%m/%d/%y %I:%M:%S %p")
  
  #round the times to the nearest minute
  d$DateTime<-format(round(d$DateTime, units="mins"), format="%Y-%m-%d %H:%M")
  
  #only # filter a moving average for every 30 minutes
  d$Light2<-filter(d$Light, sides=2, rep(1/30,30))  
  
    #d<-d[d$DateTime>start & d$DateTime<end,]
  MontereyLight[,1]<-d$DateTime[4:10003]
  MontereyLight[,i+1]<-d$Light2[4:10003]
  
}

#change the column names so that it is labeled by tidepool
colnames(MontereyLight)<-c('DateTime',files)
##for some reason the dates are showing up as characters... this converts back to datetime
MontereyLight$DateTime<-strptime(MontereyLight$DateTime, '%Y-%m-%d %H:%M')

png('../Output/LightCalibration.png', height = 1500, width = 1500, res=200)
par(mfrow=c(3,2))
# plot all data
plot(MontereyLight$DateTime, MontereyLight$TP1_light, type="l", main = 'raw light for all sensors', ylab = 'Light', xlab='Time')
for (i in 1:14){
  lines(MontereyLight$DateTime, MontereyLight[,i+2], col=i)
}

# plot common garden
plot(MontereyLight$DateTime[100:350], MontereyLight$TP1_light[100:350], type="l", 
     main = 'raw light for all sensors during common garden', ylab = 'Light', xlab='Time')
for (i in 1:14){
  lines(MontereyLight$DateTime[100:350], MontereyLight[100:350,i+2], col=i)
}

#now, take the average light at each time point and look at the deviations away when there is light
# i.e. all points that are not 0

mean.light<-rowMeans(MontereyLight[100:350,2:16])
MontereyCal<-MontereyLight[100:350,]

#plot the mean light on top
lines(MontereyLight$DateTime[100:350], mean.light, lwd=4)

#calculate the deviations of each point from the mean when there is light


Deviations<-mean.light-MontereyCal[,2:16]
#plot the deviations

#You will see that as light increases individual loggers 
# get increasingly higher or lower from the mean.  So I can't just normalize to one value, but I need to use a regression
slope<-NA
int<-NA
plot(MontereyCal[,2],Deviations$TP1_light, type = 'l', ylim=c(-250,250),
     main = 'Deviations from mean light across all sensors', xlab = 'Deviations from mean', ylab = 'Light')
for (i in 1:14){
  lines(MontereyCal[,i+1], Deviations[,i+1], col=i)
mod<-lm(Deviations[,i+1]~MontereyCal[,i+1])
newy<-predict(mod)
lines(MontereyCal[,i+1],newy, col=i, lwd='3')

# save slopes (intercepts are all nearly 0)
slope[i]<-mod$coefficients[2]
int[i]<-mod$coefficients[1]
}

# in the import temp code the new light is light + slope*light
MontereyLight2<-MontereyLight
for (i in 2:15){
MontereyLight2[,i]<- MontereyLight[,i] + (MontereyLight[,i]*(slope[i])+int[i])  
}

# plot common garden after calibration
plot(MontereyLight2$DateTime[100:350], MontereyLight2$TP1_light[100:350], type="l", 
     main = 'calibrated light for all sensors during common garden', ylab = 'Light', xlab='Time')
for (i in 1:14){
  lines(MontereyLight2$DateTime[100:350], MontereyLight2[100:350,i+2], col=i)
}

#plot the mean light on top
lines(MontereyLight2$DateTime[100:350], mean.light, lwd=4)


# plot all data again after calibration
plot(MontereyLight2$DateTime, MontereyLight2$TP1_light, type="l", main = 'Calibrated light for all sensors', ylab = 'Light', xlab='Time')
for (i in 1:14){
  lines(MontereyLight2$DateTime, MontereyLight2[,i+2], col=i)
}

dev.off()

# remove everything but the average deviations
rm(list=c('MontereyCal','MontereyLight','Deviations','mean.light','i','files','MontereyLight2'))