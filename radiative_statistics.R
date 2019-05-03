###############Statistics on radiative budget things...
#Written by A Cooper
#Last edited April 30, 2015.
#####################################################
###

library(ggplot2)
library(car)
library(compute.es)
library(effects)
library(multcomp)
library(pastecs)
library(bayesm)

#Import csv files and put into a proper dataframe.

regionList <- c('R1','R2','R3','R4','R5','R6')

for (i in 1:7){
  data <- data.frame(yearDist = 0, yearAfterDist = 0, region = 0, dist = 0, lst = 0 ,et = 0, evi = 0, x = 0, y = 0)
  
  region <- regionList[i]
  
#Insects.
  lstdir_insect <- paste("G:/BFAST_v2/LST/values_at_detected_areas/",region, "/insects/" , sep="")
  etdir_insect <- paste("G:/BFAST_v2/ET/values_at_detected_areas/",region,"/insects/", sep="")
  evidir_insect <- paste("G:/BFAST_v2/percEVI/values_at_detected_areas/",region,"/insects/",sep="")

  lstList_insect <- list.files(lstdir_insect, pattern=".txt", full.names=TRUE,recursive=FALSE,include.dirs=FALSE)
  etList_insect <- list.files(etdir_insect, pattern=".txt", full.names=TRUE,recursive=FALSE,include.dirs=FALSE)
  eviList_insect <- list.files(evidir_insect, pattern=".txt", full.names=TRUE,recursive=FALSE,include.dirs=FALSE)

#Fires.
  lstdir_fire <- paste("G:/BFAST_v2/LST/values_at_detected_areas/",region,"/fires/",sep="")
  etdir_fire <- paste("G:/BFAST_v2/ET/values_at_detected_areas/",region,"/fires/", sep="")
  evidir_fire <- paste("G:/BFAST_v2/percEVI/values_at_detected_areas/",region,"/fires/", sep="")

  lstList_fire <- list.files(lstdir_fire, pattern=".txt", full.names=TRUE,recursive=FALSE,include.dirs=FALSE)
  etList_fire <- list.files(etdir_fire, pattern=".txt", full.names=TRUE,recursive=FALSE,include.dirs=FALSE)
  eviList_fire <- list.files(evidir_fire, pattern=".txt", full.names=TRUE,recursive=FALSE,include.dirs=FALSE)
  for (i in 1:length(lstList_fire)){
    lst_insect <- read.csv(lstList_insect[i],sep=",",header=T)
    et_insect <- read.csv(etList_insect[i],sep=",",header=T)
    lst_fire <- read.csv(lstList_fire[i],sep=",",header=T)
    et_fire <- read.csv(etList_fire[i],sep=",",header=T)
    
    temp_lst_insect <- lst_insect$RASTERVALU
    temp_et_insect <- et_insect$RASTERVALU
    temp_lst_fire <- lst_fire$RASTERVALU
    temp_et_fire <- et_fire$RASTERVALU
    temp_insect_x <- lst_insect$POINT_X
    temp_insect_y <- lst_insect$POINT_Y
    temp_fire_x <- lst_fire$POINT_X
    temp_fire_y <- lst_fire$POINT_Y
    
    if (substr(lstList_insect[i],59,62)==substr(lstList_insect[i],68,71)){
        insect_file <- eviList_insect[substr(eviList_insect,57,60)==substr(lstList_insect[i],59,62)][1]
        fire_file <- eviList_fire[substr(eviList_fire,55,58)==substr(lstList_insect[i],59,62)][1]
        evi_insect <- read.csv(insect_file,sep=",",header=T) ###Fix this.
        evi_fire <- read.csv(fire_file,sep=",",header=T)
        temp_evi_fire <- evi_fire$RASTERVALU
        temp_evi_insect <- evi_insect$RASTERVALU
        
        temp <- data.frame(yearDist = rep(substr(lstList_insect[i],59,62),(length(temp_lst_fire)+length(temp_lst_insect))),
                       yearAfterDist = rep(substr(lstList_insect[i],68,71),(length(temp_lst_fire)+length(temp_lst_insect))),
                       region = rep(region,(length(temp_lst_fire)+length(temp_lst_insect))),
                       dist = c(rep('fire',length(temp_lst_fire)),rep('insect',length(temp_lst_insect))),
                       lst = c(temp_lst_fire,temp_lst_insect),
                       et = c(temp_et_fire,temp_et_insect),
                       evi = c(temp_evi_fire,temp_evi_insect),
                       x = c(temp_fire_x,temp_insect_x),
                       y = c(temp_fire_y,temp_insect_y))
    }
    if (substr(lstList_insect[i],59,62)!=substr(lstList_insect[i],68,71)){
              temp <- data.frame(yearDist = rep(substr(lstList_insect[i],59,62),(length(temp_lst_fire)+length(temp_lst_insect))),
                       yearAfterDist = rep(substr(lstList_insect[i],68,71),(length(temp_lst_fire)+length(temp_lst_insect))),
                       region = rep(region,(length(temp_lst_fire)+length(temp_lst_insect))),
                       dist = c(rep('fire',length(temp_lst_fire)),rep('insect',length(temp_lst_insect))),
                       lst = c(temp_lst_fire,temp_lst_insect),
                       et = c(temp_et_fire,temp_et_insect),
                       evi = rep(0,length(temp_lst_fire)+length(temp_lst_insect)),
                       x = c(temp_fire_x,temp_insect_x),
                       y = c(temp_fire_y,temp_insect_y))
    }
    
    data <- rbind(data,temp)
  }
  data_name <- paste('data',region,sep="")
  assign(data_name,data)
}

dataMT <- dataMT[2:nrow(dataMT),]
dataID <- dataID[2:nrow(dataID),]
dataR2 <- dataR2[2:nrow(dataR2),]
dataR3 <- dataR3[2:nrow(dataR3),]
dataR4 <- dataR4[2:nrow(dataR4),]
dataR5 <- dataR5[2:nrow(dataR5),]
dataR6 <- dataR6[2:nrow(dataR6),]

data <- rbind(dataMT,dataID,dataR2,dataR3,dataR4,dataR5,dataR6)

data <- data[data$lst!=-9999,]
data <- data[data$et!=-9999,]
data <- data[data$evi!=-9999,]
data <- data[!is.na(data$lst),]
data <- data[!is.na(data$et),]
data <- data[!is.na(data$evi),]

###Read in as csv from now on.
data <- read.csv('G:/BFAST_v2/RadiativeBalance//disturbance_info_actual_lstet.csv',header=T)

data$region <- as.character(data$region)
data$region[data$region=='MT' | data$region=='ID'] <- 'R1'

data$region <- as.factor(data$region)
data$dist <- as.factor(data$dist)

#Get data frame of lst and et for nondisturbance pixels.

for (i in 1:7){
  data <- data.frame(year = 0, region = 0, lst = 0 ,et = 0)
  
  region <- regionList[i]
  
  lstdir <- paste("G:/BFAST_v2/LST/values_nondisturbance/",region, "/" , sep="")
  etdir <- paste("G:/BFAST_v2/ET/values_nondisturbance/",region,"/", sep="")

  lstList <- list.files(lstdir, pattern=".txt", full.names=TRUE,recursive=FALSE,include.dirs=FALSE)
  etList <- list.files(etdir, pattern=".txt", full.names=TRUE,recursive=FALSE,include.dirs=FALSE)

  for (i in 1:11){
    lst <- read.csv(lstList[i],sep=",",header=T)
    et <- read.csv(etList[i],sep=",",header=T)
    
    temp_lst <- lst$RASTERVALU
    temp_et <- et$RASTERVALU
    
    temp <- data.frame(year = rep((2001+i),length(temp_lst)),
                       region = rep(region,length(temp_lst)),
                       lst = temp_lst,
                       et = temp_et)
    
    data <- rbind(data,temp)
  }
  data_name <- paste('dataNoDist',region,sep="")
  assign(data_name,data)
}

dataNoDistMT <- dataNoDistMT[2:nrow(dataNoDistMT),]
dataNoDistID <- dataNoDistID[2:nrow(dataNoDistID),]
dataNoDistR2 <- dataNoDistR2[2:nrow(dataNoDistR2),]
dataNoDistR3 <- dataNoDistR3[2:nrow(dataNoDistR3),]
dataNoDistR4 <- dataNoDistR4[2:nrow(dataNoDistR4),]
dataNoDistR5 <- dataNoDistR5[2:nrow(dataNoDistR5),]
dataNoDistR6 <- dataNoDistR6[2:nrow(dataNoDistR6),]

dataNoDist <- rbind(dataNoDistMT,dataNoDistID,dataNoDistR2,dataNoDistR3,dataNoDistR4,dataNoDistR5,dataNoDistR6)

dataNoDist <- dataNoDist[dataNoDist$lst!=-9999,]
dataNoDist <- dataNoDist[dataNoDist$et!=-9999,]
dataNoDist <- dataNoDist[!is.na(dataNoDist$lst),]
dataNoDist <- dataNoDist[!is.na(dataNoDist$et),]

###read in as csv.
dataNoDist <- read.csv('G:/BFAST_v2/RadiativeBalance/undisturbed_info_actual_lstet.csv',header=T)

dataNoDist$region <- as.character(dataNoDist$region)
dataNoDist$region[dataNoDist$region=='MT' | dataNoDist$region=='ID'] <- 'R1'

dataNoDist$region <- as.factor(dataNoDist$region)

#Shrink dataNoDist to same size as disturbance data.
dataNoDist <- dataNoDist[sample(1:nrow(dataNoDist), nrow(data),replace=FALSE),]

#Create actual difference values. ###You can skip this from now on.
#LST.
regionList <- c('R1','R2','R3','R4','R5','R6')

data$lst_diff <- 0
for (i in 1:6){
  for (j in 2002:2012){
    data$lst_diff[data$region==regionList[i] & data$yearAfterDist==j] <- data$lst[data$region==regionList[i] & data$yearAfterDist==j] - mean(dataNoDist$lst[dataNoDist$region==regionList[i] & dataNoDist$year==j])
  }
}

#ET.
data$et <- 0
for (i in 1:6){
  for (j in 2002:2012){
    data$et[data$region==regionList[i] & data$yearAfterDist==j] <- data$et[data$region==regionList[i] & data$yearAfterDist==j] - mean(dataNoDist$et[dataNoDist$region==regionList[i] & dataNoDist$year==j])
  }
}

#Create indexed time after disturbance variable.
data$yearAfterIndex <- as.numeric(data$yearAfterDist) - as.numeric(data$yearDist)

###############Plot lst and et differences by disturbance.
#Merge datasets.
#Create other columns for nondisturbance data.
dataNoDist$dist <- 'none'
dataNoDist$x <- 'na'
dataNoDist$y <- 'na'
dataNoDist$yearAfterDist <- 'na'
dataNoDist$yearAfterIndex <- 'na'
dataNoDist$evi <- 'na'

names(dataNoDist)[2] <- 'yearDist'
dataNoDist <- dataNoDist[,2:11]

data<-data[,2:11]
comboData <- rbind(data[data$yearAfterIndex==0,],dataNoDist)

#Get average lst and et differences at each time for each region.
lst_averages <- aggregate(comboData$lst,by=list(comboData$region, comboData$dist),FUN="mean")
et_averages <- aggregate(comboData$et,by=list(comboData$region, comboData$dist),FUN="mean")
names(lst_averages) <- c('region','dist','lst')
names(et_averages) <- c('region','dist','et')

#Plot.
#FIGURES!!!
ggplot(data=comboData, aes(x=region , y=lst, fill=dist)) + geom_boxplot() +
  xlab('Region') +
  ylab('LST (K)') +
  theme(axis.title.x = element_text(face="bold", size=20),
          axis.text.x  = element_text(size=16),
          axis.title.y = element_text(face="bold",size=20),
          axis.text.y = element_text(size=16),
          legend.title = element_text(face="bold",size=16),
          legend.text = element_text(size=16)) +
  scale_x_discrete(breaks=c("R1", "R2",'R3','R4','R5','R6'),
                      labels=c("MT/ID", "WY/CO", "AZ/NM","UT/NV","CA","WA/OR")) +
  scale_fill_discrete(name = "Disturbance",
                      breaks = c("fire","insect","none"),
                      labels = c("Fire","Insects","None"))
ggplot(data=comboData, aes(x=region , y=et, fill=dist)) + geom_boxplot() +
  xlab('Region') +
  ylab('ET (mm/mo)') +
    theme(axis.title.x = element_text(face="bold", size=20),
          axis.text.x  = element_text(size=16),
          axis.title.y = element_text(face="bold",size=20),
          axis.text.y = element_text(size=16),
          legend.title = element_text(face="bold",size=16),
          legend.text = element_text(size=16)) +
  scale_x_discrete(breaks=c("R1", "R2",'R3','R4','R5','R6'),
                      labels=c("MT/ID", "WY/CO", "AZ/NM","UT/NV","CA","WA/OR")) +
  scale_fill_discrete(name = "Disturbance",
                      breaks = c("fire","insect","none"),
                      labels = c("Fire","Insects","None"))

#################Test (ANOVA) differences between disturbances by type and region.

#By disturbance. Significant difference (p<0.001)
#STATS!!!
dist_lst_aov <- aov(data$lst[data$yearAfterIndex==0]~data$dist[data$yearAfterIndex==0])
dist_et_aov <- aov(data$et[data$yearAfterIndex==0]~data$dist[data$yearAfterIndex==0])
dist_evi_aov <- aov(data$evi[data$yearAfterIndex==0]~data$dist[data$yearAfterIndex==0])

#By region. Significant difference (p<0.001)
reg_lst_aov <- aov(data$lst[data$yearAfterIndex==0]~data$region[data$yearAfterIndex==0])
reg_et_aov <- aov(data$et[data$yearAfterIndex==0]~data$region[data$yearAfterIndex==0])
reg_evi_aov <- aov(data$evi[data$yearAfterIndex==0]~data$region[data$yearAfterIndex==0])

#Boxplot summary of differences in lst, evi, et by region and disturbance.
#FIGURES!!!
ggplot(data[data$yearAfterIndex==0,], aes(x=region, y=lst, fill=dist)) + geom_boxplot() +
  xlab('Region') +
  ylab('Difference between Disturbance and Non-Disturbance LST')
ggplot(data[data$yearAfterIndex==0,], aes(x=region, y=et, fill=dist)) + geom_boxplot() +
  xlab('Region') +
  ylab('Difference between Disturbance and Non-Disturbance ET')
ggplot(data[data$yearAfterIndex==0,], aes(x=region, y=evi, fill=dist)) + geom_boxplot() +
  xlab('Region') +
  ylab('Difference between Disturbance and Non-Disturbance %EVI Decline')

#Differences between region lst,et, keeping severity the same. 
#Summary of the glht function provides t-tests for each combination.
##STATS!!!!!!
dist_lst_cov <- aov(lst~dist + evi, data=data[data$yearAfterIndex==0,])
Anova(dist_lst_cov, type='III')
dist_lst_posthocs <- glht(dist_lst_cov,linfct = mcp(dist = "Tukey"))
summary(dist_lst_posthocs)
confint(dist_lst_posthocs)

reg_lst_cov <- aov(lst~region + evi, data=data[data$yearAfterIndex==0,])
Anova(reg_lst_cov, type='III')
reg_lst_posthocs <- glht(reg_lst_cov,linfct = mcp(region = "Tukey"))
summary(reg_lst_posthocs)
confint(reg_lst_posthocs)

dist_et_cov <- aov(et~dist + evi, data=data[data$yearAfterIndex==0,])
Anova(dist_et_cov, type='III')
dist_et_posthocs <- glht(dist_et_cov,linfct = mcp(dist = "Tukey"))
summary(dist_et_posthocs)
confint(dist_et_posthocs)

reg_et_cov <- aov(et~region + evi, data=data[data$yearAfterIndex==0,])
Anova(reg_et_cov, type='III')
reg_et_posthocs <- glht(reg_et_cov,linfct = mcp(region = "Tukey"))
summary(reg_et_posthocs)
confint(reg_et_posthocs)


#################Recovery analysis.

###First, get p values of t tests.
#Blank variables.
data$lst_tValue <- 0
data$et_tValue <- 0

regionList <-c('R1','R2','R3','R4','R5','R6')

#Get p values for each pixel combo.
for (i in 1:6){
  for (j in 2002:2012){
    for (k in 1:2){
      for (l in 0:10){
        region = regionList[i]
        yearAfterDist = j
        yearAfterIndex = l
        dist = c('fire','insect')[k]
        tempDist <- data[data$region==region & data$yearAfterIndex==yearAfterIndex & data$dist==dist & data$yearAfterDist==yearAfterDist,]
        tempUnDist <- dataNoDist[dataNoDist$yearDist==j & dataNoDist$region==region,]
        if (nrow(tempDist)<=10){
          lst_t_temp <- 1
          et_t_temp <- 1
        }
        if (nrow(tempDist)>10){
          lst_t_temp <- t.test(tempDist$lst,tempUnDist$lst,"greater")$p.value
          et_t_temp <- t.test(tempDist$et,tempUnDist$et,"less")$p.value
        }  
        data$lst_tValue[data$region==region & data$yearAfterDist==yearAfterDist & data$dist==dist & data$yearAfterIndex==yearAfterIndex] <- lst_t_temp
        data$et_tValue[data$region==region & data$yearAfterDist==yearAfterDist & data$dist==dist & data$yearAfterIndex==yearAfterIndex] <- et_t_temp
      }
    }
  }
}

#Create new dataset with only significant p values. This is years dist is still significant,
#need to add 1 to the recovery time number from this.
dataRecov <- data

#Get maximum year after disturbance.

temp_lst_rec <- aggregate(dataRecov$yearAfterIndex[dataRecov$lst_tValue > 0.05],
                  by=list(dataRecov$yearDist[dataRecov$lst_tValue > 0.05],
                          dataRecov$region[dataRecov$lst_tValue > 0.05],
                          dataRecov$dist[dataRecov$lst_tValue > 0.05],dataRecov$x[dataRecov$lst_tValue > 0.05],
                          dataRecov$y[dataRecov$lst_tValue > 0.05]), FUN = min)
temp_et_rec <- aggregate(dataRecov$yearAfterIndex[dataRecov$et_tValue > 0.05],
                  by=list(dataRecov$yearDist[dataRecov$et_tValue > 0.05],
                          dataRecov$region[dataRecov$et_tValue > 0.05],
                          dataRecov$dist[dataRecov$et_tValue > 0.05],dataRecov$x[dataRecov$et_tValue > 0.05],
                          dataRecov$y[dataRecov$et_tValue > 0.05]), FUN = min)

names(temp_lst_rec) <- c('yearDist','region','dist','x','y','lst_rec')
names(temp_et_rec) <- c('yearDist','region','dist','x','y','et_rec')

temp_rec <- merge(temp_lst_rec,temp_et_rec,by=c('x','y','yearDist','region','dist'))

dataRecov2 <- temp_rec

temp <- aggregate(data$evi,by=list(data$x,data$y,data$yearDist,data$region,data$dist),FUN=min)
names(temp) <- c('x','y','yearDist','region','dist','evi')

dataRecov2 <- merge(dataRecov2,temp,by=c('x','y','yearDist','region','dist'))

#Write out to save some time.
write.table(dataRecov2,'G:/BFAST_v2/RadiativeBalance/recovery_data_actual_lstet.txt',sep='\t')

########Summarize recovery time by region (bar graphs with error bars).
#Calculate standard error. CHANGE TO INCLUDE YRS AFTER RECOVERY INDEX.
dataRecov2$lst_rec_se <- 0
dataRecov2$et_rec_se <- 0
for (i in 1:6){
  for (j in 2002:2012){
    for (k in 1:2){
      region = regionList[i]
      yearDist = j
      dist = c('fire','insect')[k]
      lst_se <- sd(dataRecov2$lst_rec[dataRecov2$region==region & dataRecov2$yearDist==yearDist & dataRecov2$dist==dist])/sqrt(length(dataRecov2$lst_rec[dataRecov2$region==region & dataRecov2$yearDist==yearDist & dataRecov2$dist==dist]))
      et_se <- sd(dataRecov2$et_rec[dataRecov2$region==region & dataRecov2$yearDist==yearDist & dataRecov2$dist==dist])/sqrt(length(dataRecov2$et_rec[dataRecov2$region==region & dataRecov2$yearDist==yearDist & dataRecov2$dist==dist]))
      dataRecov2$lst_rec_se[dataRecov2$region==region & dataRecov2$yearDist==yearDist & dataRecov2$dist==dist] <- lst_se
      dataRecov2$et_rec_se[dataRecov2$region==region & dataRecov2$yearDist==yearDist & dataRecov2$dist==dist] <- et_se
    }
  }
}

write.table(dataRecov2,'G:/BFAST_v2/RadiativeBalance/recovery_data_actual_lstet.txt',sep='\t')
#Read in table in the future.
dataRecov2 <- read.table('G:/BFAST_v2/RadiativeBalance/recovery_data_actual_lstet.txt',sep='\t',header=T)


#Plot.
lst_agg <- aggregate(dataRecov2$lst_rec,by=list(dataRecov2$region,dataRecov2$dist),FUN=mean)
lst_se_agg <- aggregate(dataRecov2$lst_rec_se,by=list(dataRecov2$region,dataRecov2$dist),FUN=mean)
et_agg <- aggregate(dataRecov2$et_rec,by=list(dataRecov2$region,dataRecov2$dist),FUN=mean)
et_se_agg <- aggregate(dataRecov2$et_rec_se,by=list(dataRecov2$region,dataRecov2$dist),FUN=mean)
names(lst_agg)<-c('region','dist','lst_rec')
names(lst_se_agg)<-c('region','dist','lst_rec_se')
names(et_agg)<-c('region','dist','et_rec')
names(et_se_agg)<-c('region','dist','et_rec_se')

temp <- merge(lst_agg,et_agg,by=c('region','dist'))
temp <- merge(temp,lst_se_agg,by=c('region','dist'))
temp <- merge(temp,et_se_agg,by=c('region','dist'))

#FIGURES!!!
ggplot(temp, aes(x=region, y=lst_rec, fill=dist)) + 
    geom_bar(stat="identity",colour = "black",position=position_dodge()) +
    geom_errorbar(aes(ymin=lst_rec-lst_rec_se, ymax=lst_rec+lst_rec_se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) + 
  theme(legend.title = element_blank(),legend.text = element_text(size=24)) + 
  scale_fill_discrete(breaks = c('fire','insect'), labels = c('Fire','Insects')) +
  xlab('Region') + 
  ylab('Years to LST Recovery')

ggplot(temp, aes(x=region, y=et_rec, fill=dist)) + 
    geom_bar(stat="identity",colour = "black",position=position_dodge()) +
    geom_errorbar(aes(ymin=et_rec-et_rec_se, ymax=et_rec+et_rec_se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) + 
  theme(legend.title = element_blank(),legend.text = element_text(size=24)) + 
  scale_fill_discrete(breaks = c('fire','insect'), labels = c('Fire','Insects')) + 
  xlab('Region') + 
  ylab('Years to ET Recovery')

##STATS!!!!!!
dist_lst_cov <- aov(lst_rec~dist + evi, data=dataRecov2)
Anova(dist_lst_cov, type='III')
dist_lst_posthocs <- glht(dist_lst_cov,linfct = mcp(dist = "Tukey"))
summary(dist_lst_posthocs)
confint(dist_lst_posthocs)

reg_lst_cov <- aov(lst_rec~region + evi, data=dataRecov2[dataRecov2$dist=='fire',])
Anova(reg_lst_cov, type='III')
reg_lst_posthocs <- glht(reg_lst_cov,linfct = mcp(region = "Tukey"))
summary(reg_lst_posthocs)
confint(reg_lst_posthocs)

dist_et_cov <- aov(et_rec~dist + evi, data=dataRecov2)
Anova(dist_et_cov, type='III')
dist_et_posthocs <- glht(dist_et_cov,linfct = mcp(dist = "Tukey"))
summary(dist_et_posthocs)
confint(dist_et_posthocs)

reg_et_cov <- aov(et_rec~region + evi, data=dataRecov2[dataRecov2$dist=='fire',])
Anova(reg_et_cov, type='III')
reg_et_posthocs <- glht(reg_et_cov,linfct = mcp(region = "Tukey"))
summary(reg_et_posthocs)
confint(reg_et_posthocs)

#######LST and ET differences over time (with error?)
recovTime <- aggregate(data$lst, by = list(data$yearDist,data$yearAfterDist,data$region,data$dist,data$yearAfterIndex),FUN = mean)
temp <- aggregate(data$et, by = list(data$yearDist,data$yearAfterDist,data$region,data$dist,data$yearAfterIndex),FUN = mean)
names(recovTime) <- c('yearDist','yearAfterDist','region','dist','yearAfterIndex','lst')
names(temp) <- c('yearDist','yearAfterDist','region','dist','yearAfterIndex','et')
recovTime <- merge(recovTime,temp,by=c('yearDist','yearAfterDist','region','dist','yearAfterIndex'))

recovTime$lst_diff <- NULL
recovTime$et_diff <- NULL
for (i in 1:nrow(recovTime)){
  year <- recovTime$yearAfterDist[i]
  region <- recovTime$region[i]
  recovTime$lst_diff[i] <- recovTime$lst[i] - mean(dataNoDist$lst[dataNoDist$region==region & 
                                        dataNoDist$yearDist==year])
  recovTime$et_diff[i] <- recovTime$et[i] - mean(dataNoDist$et[dataNoDist$region==region & 
                                        dataNoDist$yearDist==year])
}

lst <- aggregate(dataRecov2$lst,by=list(dataRecov2$yearAfterIndex,dataRecov2$dist,dataRecov2$region,dataRecov2$yearDist),FUN = mean)
names(lst) <- c('yearAfterIndex','dist','region','yearDist','lst')
et <- aggregate(dataRecov2$et,by=list(dataRecov2$yearAfterIndex,dataRecov2$dist,dataRecov2$region,dataRecov2$yearDist),FUN = mean)
names(et) <- c('yearAfterIndex','dist','region','yearDist','et')
ggplot(data=dataRecov2[et$dist=='fire' & et$yearDist=='2011',], aes(x=yearAfterIndex, y=et, group=region, colour=region)) +
    geom_line() +
    geom_point()

test <- aggregate(data$lst,by=list(data$region,data$dist,data$yearAfterIndex),FUN=mean)

############Maybe want this???
average_lst_impact <- aggregate(recovTime$lst_diff,by=list(recovTime$region,recovTime$dist,recovTime$yearAfterIndex),FUN=mean)
average_et_impact <- aggregate(recovTime$et_diff,by=list(recovTime$region,recovTime$dist,recovTime$yearAfterIndex),FUN=mean)

########Model LST and ET recovery time.
#Indicator plots for recovery.
par(mfrow=c(3,3),oma=c(0,0,0,0),            # Creates a 3x3 graphics window with no outer margin
  mar=c(0,0,3,0))                           #   (oma), and only a top margin (3) on each plot.
thresh <- round(quantile(dataRecov2$lst_rec,   # Computes indicator threshholds at each decile of
  probs=seq(.1,.9,.1)),2)                   #   of the coal ash percentages (10%, 20%, ..., 90%)
for (i in 1:9){                             # Loops through the 9 indicator plots.
  datadf <- data.frame(dataRecov2$x[dataRecov2$dist=='fire'],dataRecov2$y[dataRecov2$dist=='fire'], # Creates a data frame with x- and y-coordinates, and
    dataRecov2$lst_rec[dataRecov2$dist=='fire']>thresh[i])              #   indicators (1,0) if the coal ash % > threshhold
  plot(dataRecov2$x[dataRecov2$dist=='fire'],dataRecov2$y[dataRecov2$dist=='fire'],col=c("black","red")[as.factor(dataRecov2$lst_rec[dataRecov2$dist=='fire']>thresh[i])],main=paste("LST recovery % >",thresh[i],sep=" "))
}    


#Aggregate so only one version of a point.\
dataRecov2$lst_rec <- as.integer(dataRecov2$lst_rec)
dataRecov2$et_rec <- as.integer(dataRecov2$et_rec)

dataModel_lst <- aggregate(dataRecov2$evi, by=list(dataRecov2$x,dataRecov2$y,dataRecov2$yearDist,dataRecov2$dist,dataRecov2$lst_rec),FUN = max)
names(dataModel_lst) <- c('x','y','yearDist','dist','lst_rec','evi')
dataModel_lst$evi <- dataModel_lst$evi * -1
dataModel_lst$lst_recLog <- log(dataModel_lst$lst_rec)
dataModel_lst <- dataModel_lst[dataModel_lst$lst_recLog!=-Inf,]


###Create model.
indexes <- sample(1:nrow(dataModel_lst), size=round(0.7*nrow(dataModel_lst)))
training <- dataModel_lst[indexes,]
testing <- dataModel_lst[-indexes,]
lst_model <- glm(lst_recLog ~ x + y + evi + dist, data=training)

###Look at model fit.
summary(lst_model)
testing$model<-predict(lst_model,testing)
plot(lst_model$residuals)

plot(testing$lst_rec~exp(testing$model))
abline(a=0,b=1, col="blue")
validate<-lm(testing$lst_rec~exp(testing$model))
summary(validate) #0.8085 R^2
sqrt(mean((exp(testing$model)-testing$lst_rec)^2 , na.rm = TRUE )) #RMSE  of 0.3981494

qqplot(exp(testing$model),testing$lst_rec)
abline(a=0,b=1,col="blue")

#####Model ET.
dataModel_et <- aggregate(dataRecov2$evi, by=list(dataRecov2$x,dataRecov2$y,dataRecov2$yearDist,dataRecov2$dist,dataRecov2$et_rec),FUN = mean)
names(dataModel_et) <- c('x','y','yearDist','dist','et_rec','evi')
dataModel_et$evi <- dataModel_et$evi * -1
dataModel_et$et_recLog <- log(dataModel_et$et_rec)
dataModel_et <- dataModel_et[dataModel_et$et_recLog!=-Inf,]


###Create model.
indexes <- sample(1:nrow(dataModel_et), size=round(0.7*nrow(dataModel_et)))
training <- dataModel_et[indexes,]
testing <- dataModel_et[-indexes,]
et_model <- glm(et_recLog ~ x + y + evi + dist, data=training)

###Look at model fit.
summary(et_model)
testing$model<-predict(et_model,testing)
plot(et_model$residuals)

plot(testing$et_recLog~testing$model)
abline(a=0,b=1, col="blue")
validate<-lm(testing$et_recLog~testing$model)
summary(validate) #0.8432 R^2
sqrt(mean((testing$model-testing$et_recLog)^2 , na.rm = TRUE )) #RMSE  of 0.1672786

qqplot(testing$model,testing$et_recLog)
abline(a=0,b=1,col="blue")




###Predict over random points.
#Import points.
points <- read.table('G:/BFAST_v2/RadiativeBalance/randomPoints.txt',sep=',',header=T)
#Set variable names and levels of evi.
points <- points[,3:4]
names(points) <- c('x','y')

#Predict.
#low evi, fire.
points_low_fire <- points
points_low_fire$evi <- 0.08
points_low_fire$dist <- 'fire'

points_low_fire$et <- exp(predict(et_model,points_low_fire))
points_low_fire$lst <- exp(predict(lst_model,points_low_fire))

#mid evi, fire.
points_mid_fire <- points
points_mid_fire$evi <- 0.25
points_mid_fire$dist <- 'fire'

points_mid_fire$et <- exp(predict(et_model,points_mid_fire))
points_mid_fire$lst <- exp(predict(lst_model,points_mid_fire))

#high evi, fire.
points_high_fire <- points
points_high_fire$evi <- 0.60
points_high_fire$dist <- 'fire'

points_high_fire$et <- exp(predict(et_model,points_high_fire))
points_high_fire$lst <- exp(predict(lst_model,points_high_fire))

#low, insects.
points_low_insects <- points
points_low_insects$evi <- 0.08
points_low_insects$dist <- 'insect'

points_low_insects$et <- exp(predict(et_model,points_low_insects))
points_low_insects$lst <- exp(predict(lst_model,points_low_insects))

#mid evi, insects.
points_mid_insects <- points
points_mid_insects$evi <- 0.25
points_mid_insects$dist <- 'insect'

points_mid_insects$et <- exp(predict(et_model,points_mid_insects))
points_mid_insects$lst <- exp(predict(lst_model,points_mid_insects))

#high evi, insects.
points_high_insects <- points
points_high_insects$evi <- 0.60
points_high_insects$dist <- 'insect'

points_high_insects$et <- exp(predict(et_model,points_high_insects))
points_high_insects$lst <- exp(predict(lst_model,points_high_insects))

#Interpolate image.
library(geoR)
library(maps)
library(akima)
library(rgdal)
library(gstat)
library(lattice)
library(spatstat)
library(plyr)

#high, fire
map('state',c('washington','oregon','california','idaho','montana','nevada','utah','wyoming','colorado','arizona','new mexico'))
temp_interp<-interp(points_high_fire$x,points_high_fire$y,points_high_fire$et,duplicate="mean")
image(temp_interp,xlab="Longitude",ylab="Latitude",col=rev(heat.colors(24)),cex.main=1.6,cex.lab=1.4,add=T)
map('state',c('washington','oregon','california','idaho','montana','nevada','utah','wyoming','colorado','arizona','new mexico'),add=T)
image.legend(-107,49,zlim=range(points_high_fire$et),lwd=1,bty="n",cex=0.75,col=rev(heat.colors(24)))

#mid, fire
map('state',c('washington','oregon','california','idaho','montana','nevada','utah','wyoming','colorado','arizona','new mexico'))
temp_interp<-interp(points_mid_fire$x,points_mid_fire$y,points_mid_fire$lst,duplicate="mean")
image(temp_interp,xlab="Longitude",ylab="Latitude",col=rev(heat.colors(24)),cex.main=1.6,cex.lab=1.4,add=T)
map('state',c('washington','oregon','california','idaho','montana','nevada','utah','wyoming','colorado','arizona','new mexico'),add=T)
image.legend(-107,49,zlim=range(points_mid_fire$et),lwd=1,bty="n",cex=0.75,col=rev(heat.colors(24)))

#low,fire
map('state',c('washington','oregon','california','idaho','montana','nevada','utah','wyoming','colorado','arizona','new mexico'))
temp_interp<-interp(points_low_fire$x,points_low_fire$y,points_low_fire$et,duplicate="mean")
image(temp_interp,xlab="Longitude",ylab="Latitude",col=rev(heat.colors(24)),cex.main=1.6,cex.lab=1.4,add=T)
map('state',c('washington','oregon','california','idaho','montana','nevada','utah','wyoming','colorado','arizona','new mexico'),add=T)
image.legend(-107,49,zlim=range(points_low_fire$et),lwd=1,bty="n",cex=0.75,col=rev(heat.colors(24)))


#high,insect
map('state',c('washington','oregon','california','idaho','montana','nevada','utah','wyoming','colorado','arizona','new mexico'))
temp_interp<-interp(points_high_insects$x,points_high_insects$y,points_high_insects$et,duplicate="mean")
image(temp_interp,xlab="Longitude",ylab="Latitude",col=rev(heat.colors(24)),cex.main=1.6,cex.lab=1.4,add=T)
map('state',c('washington','oregon','california','idaho','montana','nevada','utah','wyoming','colorado','arizona','new mexico'),add=T)
image.legend(-107,49,zlim=range(points_high_insects$et),lwd=1,bty="n",cex=0.75,col=rev(heat.colors(24)))

#mid,insect
map('state',c('washington','oregon','california','idaho','montana','nevada','utah','wyoming','colorado','arizona','new mexico'))
temp_interp<-interp(points_mid_insects$x,points_mid_insects$y,points_mid_insects$et,duplicate="mean")
image(temp_interp,xlab="Longitude",ylab="Latitude",col=rev(heat.colors(24)),cex.main=1.6,cex.lab=1.4,add=T)
map('state',c('washington','oregon','california','idaho','montana','nevada','utah','wyoming','colorado','arizona','new mexico'),add=T)
image.legend(-107,49,zlim=range(points_mid_insects$et),lwd=1,bty="n",cex=0.75,col=rev(heat.colors(24)))


#low,insect
map('state',c('washington','oregon','california','idaho','montana','nevada','utah','wyoming','colorado','arizona','new mexico'))
temp_interp<-interp(points_low_insects$x,points_low_insects$y,points_low_insects$lst,duplicate="mean")
image(temp_interp,xlab="Longitude",ylab="Latitude",col=rev(heat.colors(24)),cex.main=1.6,cex.lab=1.4,add=T)
map('state',c('washington','oregon','california','idaho','montana','nevada','utah','wyoming','colorado','arizona','new mexico'),add=T)
image.legend(-107,49,zlim=range(points_low_insects$et),lwd=1,bty="n",cex=0.75,col=rev(heat.colors(24)))


#Plot.
