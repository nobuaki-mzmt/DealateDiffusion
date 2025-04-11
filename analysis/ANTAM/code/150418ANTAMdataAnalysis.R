## ANTAM trajectory data analysis
## 15_4_18

## omitting CF1 because she died at second day

## packages
library(car)
library(lme4)

Folder.path <- getwd()
f.namesplace <- list.files(Folder.path, pattern=".csv",full.names=T)

## analysis condition
FPS <- 5	## FPS using in analysis
msec <- 1/FPS * 1000	## one step / msec
steps <- 30*60*1000/msec + 1	## total steps for 30 min

## reading parameters
x = y = sec = time <-  matrix(0, ncol=length(f.namesplace), steps)
name = colony = sex = rep = day = totaldis = totaltime <- rep(0,length(f.namesplace))
angle = dis = move = speed = gostop = reorient = acceleration <- matrix(0, ncol=length(f.namesplace), steps)

## reading data
for(j in 1:length(f.namesplace)){
	d <- read.csv(f.namesplace[j],sep=",")[,1:3]
	namebegin <- regexpr("mouselog_",f.namesplace[j]) + 9
	name[j] <- substring(f.namesplace[j],namebegin, namebegin+7)
	colony[j] <- substring(f.namesplace[j],namebegin, namebegin)
	sex[j] <- substring(f.namesplace[j],namebegin+1, namebegin+1)
	rep[j] <- substring(f.namesplace[j],namebegin+2, namebegin+2)
	day[j] <- as.numeric(substring(f.namesplace[j], namebegin+7, namebegin+7))
	i <- 1500	## 何stepからのデータを使うか
	x[1,j] <- d[i,1]
	y[1,j] <- d[i,2]
	rownum <- 2
	i <- i + 1
	while(rownum < steps){
		sec[rownum, j] <- sec[rownum, j] + d[i, 3]
		if( abs(msec - sec[rownum, j]) < abs(msec - (sec[rownum, j]+d[i+1, 3]))){	# 次のsecを足したものの方が400より遠い
			x[rownum,j] <- d[i,1]
			y[rownum,j] <- d[i,2]
			time[rownum,j] <- sum(sec[,j])
			dis[rownum,j] <- ((x[rownum,j]-x[rownum-1,j])^2+(y[rownum,j]-y[rownum-1,j])^2)^0.5	# compute moving distance
			speed[rownum,j] <- dis[rownum,j] / sec[rownum,j]							# compute moving speed
			acceleration[rownum,j] <- speed[rownum,j] - speed[rownum-1,j]
			move[rownum,j] <- ((x[rownum,j]-x[1,j])^2+(y[rownum,j]-y[1,j])^2)^0.5		# compute distance from begining point
			# compute angle by vector
			Ax <- (x[rownum,j] - x[rownum-1,j])
			Bx <- (x[rownum+1,j] - x[rownum-1,j])
			Ay <- (y[rownum,j] - y[rownum-1,j])
			By <- (y[rownum+1,j] - y[rownum-1,j])
			hugo <- (Ax * By - Ay * Bx)/abs(Ax * By - Ay * Bx)
			cos <- (Ax * Bx + Ay * By) / ((Ax^2 + Ay^2)*(Bx^2 + By^2))^0.5
			angle[rownum,j] <- acos(cos) * hugo
			rownum <- rownum + 1
		}
		i <- i + 1	
	}
	totaldis[j] <- sum(dis[,j])
	totaltime[j] <- sum(sec[,j])
	print(name[j])
	#hist(r[,j]*(180/pi),180)
}




gostop[dis > 0] <- 1
ind <-8
plot(gostop[,ind]~time[,ind],type="l")
gorate <- apply(gostop,2,mean)
plot(gorate~day,col=sex2)
data.frame(sex=sex,day=day,gorate=gorate)
meanspeed <- apply(speed,2,mean)
plot(gorate~meanspeed,col=sex2)

reorient=gostop <- matrix(0, ncol=length(f.namesplace), length(steps)-1)
hist(angle)
reorient[angle > (pi/3)] <-1
plot(reorient[,8]~time[,8],type="l")
reorientrate <- apply(reorient,2,mean)
plot(reorientrate~day,col=sex2)
r<-lm(reorientrate ~ day*sex)
anova(r)

## ANTAM data analysis
## for Plot

## ♂は♀よりも方向転換角度が大きい
angle[sec>250]<-100
angle[sec<150]<-100
par(pin = c(4,4))
r <- hist(abs(angle[angle<100]),breaks=seq(0,pi,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), pch=20, xlim=c(0,3.2),
	xlab="angle", ylab = "frequency")
par(new=T)
r <- hist(abs(angle[,(1:length(sex))[sex=="M"]][angle[,(1:length(sex))[sex=="M"]] <250]),breaks=seq(0,pi,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), col="blue", pch=20, xlim=c(0,3.2),
	xlab="angle", ylab = "frequency")
par(new=T)
r <- hist(abs(angle[,(1:length(sex))[sex=="F"]][angle[,(1:length(sex))[sex=="F"]] <250]),breaks=seq(0,pi,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), col="red", pch=20, xlim=c(0,3.2),
	xlab="angle", ylab = "frequency")


## 日付と方向転換角度
par(pin = c(4,4))
r <- hist(abs(angle[,(1:length(day))[day==0]]),breaks=seq(0,pi,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), pch=20, xlim=c(0,3.2),
	xlab="angle", ylab = "frequency")
par(new=T)
r <- hist(abs(angle[,(1:length(day))[day==1]]),breaks=seq(0,pi,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), col="red", pch=20, xlim=c(0,3.2),
	xlab="angle", ylab = "frequency")
par(new=T)
r <- hist(abs(angle[,(1:length(day))[day==2]]),breaks=seq(0,pi,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), col="blue", pch=20, xlim=c(0,3.2),
	xlab="angle", ylab = "frequency")
par(new=T)
r <- hist(abs(angle[,(1:length(day))[day==3]]),breaks=seq(0,pi,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), col="green", pch=20, xlim=c(0,3.2),
	xlab="angle", ylab = "frequency")

angle2<-angle
angle2[sec>220] <- NaN
angle2 <- matrix(angle2,ncol=length(f.namesplace))
angle2[sec<180] <- NaN
angle2 <- matrix(angle2,ncol=length(f.namesplace))
hist(angle2)
rep2 <- as.numeric(rep)
rep2[(1:length(sex))[sex == "M"]] <- rep2[(1:length(sex))[sex == "M"]] +3
meanangle <- apply(abs(angle2),2,mean,na.rm=T)
# r <- lm(meanangle~day*sex) #+ colony/rep2)
r <- lmer(meanangle~day*sex + (1|colony/rep2))
Anova(r)#r 
plot(meanangle~day,col=sex2)

## ♂と♀の速度
par(pin = c(4,4))
r <- hist(abs(speed),breaks=seq(0,10,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), pch=20, xlim=c(0,10),
	xlab="speed", ylab = "frequency")
par(new=T)
r <- hist(abs(speed[,(1:length(sex))[sex=="M"]]),breaks=seq(0,10,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), col="blue", pch=20, xlim=c(0,10),
	xlab="speed", ylab = "frequency")
par(new=T)
r <- hist(abs(speed[,(1:length(sex))[sex=="F"]]),breaks=seq(0,10,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), col="red", pch=20, xlim=c(0,10),
	xlab="speed", ylab = "frequency")


## 日付と速度
par(pin = c(4,4))
r <- hist(abs(speed[,(1:length(day))[day==0]]),breaks=seq(0,6,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), pch=20, xlim=c(0,6),
	xlab="speed", ylab = "frequency")
par(new=T)
r <- hist(abs(speed[,(1:length(day))[day==1]]),breaks=seq(0,6,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), col="red", pch=20, xlim=c(0,6),
	xlab="speed", ylab = "frequency")
par(new=T)
r <- hist(abs(speed[,(1:length(day))[day==2]]),breaks=seq(0,6,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), col="blue", pch=20, xlim=c(0,6),
	xlab="speed", ylab = "frequency")
par(new=T)
r <- hist(abs(speed[,(1:length(day))[day==3]]),breaks=seq(0,6,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), col="green", pch=20, xlim=c(0,6),
	xlab="speed", ylab = "frequency")

rep2 <- as.numeric(rep)
rep2[(1:length(sex))[sex == "M"]] <- rep2[(1:length(sex))[sex == "M"]] +3
meanspeed <- apply(abs(speed),2,mean,na.rm=T)
#r <- lm(meanspeed~day*sex + colony/rep2)
r <- lmer(meanspeed~day*sex + (1|colony/rep2))
Anova(r)
plot(meanspeed~day,col=sex2)

## もっとも離れた場所
mostaway <- apply(move, 2, max)	#　各 trajectory の most away point
plot(mostaway ~ totaldis,col=as.numeric(day)+1,pch=20)
plot(mostaway ~ totaldis,col=sex2)
par(pin = c(4,4))
r <- hist(abs(move[,(1:length(day))[day==1]]),breaks=50,plot=T)
r <- hist(abs(move),breaks=50,plot=T)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-4,1), pch=20, xlim=c(0,500000),
	xlab="distance from begining point", ylab = "frequency")


## speed~angle 
angle2 <- angle#[,(1:length(sex))[sex=="F"]]
speed2 <- speed#[,(1:length(sex))[sex=="F"]]
dim(angle2) <- c(length(angle2),1)
dim(speed2) <- c(length(speed2),1)
speed3 <- round(speed2,2)
#plot(abs(angle2)~speed3)
angleperspeed <- tapply(abs(angle2),speed3,mean,na.rm=TRUE)
speedseq <- seq(0,max(speed3),length=length(angleperspeed))
plot(angleperspeed~speedseq,pch=1,xlim=c(0,6),ylim=c(0,1.5),col=1)

## max speed
sex2<-sex
sex2[sex2=="M"]<-1
sex2[sex2=="F"]<-2
plot(apply(speed,2,max)~apply(angle,2,max,na.rm=T),col=sex2)



ind<-45
library(rgl)
plot3d(x[1:steps,ind], y[1:steps,ind], time[1:steps,ind])


