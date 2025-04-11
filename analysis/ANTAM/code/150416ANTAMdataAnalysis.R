## ANTAM trajectory data analysis
## 15_4_16

Folder.path <- getwd()
f.namesplace <- list.files(Folder.path, pattern=".csv",full.names=T)

steps <- 3000:12001 * 5 
# total = more than 67000 steps (about 40 msec / 1 step)
# 1:12001 * 5 = about 200msec / 1 step

# reading parameters ()
name = colony = sex = rep = day = totaldis = totaltime <- rep(0,length(f.namesplace))
x = y = sec <-  matrix(0, ncol=length(f.namesplace), length(steps))
angle = dis = move = time = speed = gostop = reorient <- matrix(0, ncol=length(f.namesplace), length(steps)-1)

for(j in 1:length(f.namesplace)){
	name[j] <- substring(f.namesplace[j],66,73)
	colony[j] <- substring(f.namesplace[j],66,66)
	sex[j] <- substring(f.namesplace[j],67,67)
	rep[j] <- substring(f.namesplace[j],68,68)
	day[j] <- as.numeric(substring(f.namesplace[j],73,73))
	d <- read.csv(f.namesplace[j],sep=",")[,1:3]
	x[,j] <- d[steps,1]
	y[,j] <- d[steps,2]
	sec[,j] <- d[steps,3]
	for(i in 2:(length(steps)-1)){
		time[i,j] <- sec[i,j] + time[i-1,j]							# time in msec
		dis[i,j] <- ((x[i,j]-x[i-1,j])^2+(y[i,j]-y[i-1,j])^2)^0.5	# compute moving distance
		speed[i,j] <- dis[i,j] / sec[i,j]							# compute moving speed
		move[i,j] <- ((x[i,j]-x[1,j])^2+(y[i,j]-y[1,j])^2)^0.5		# compute distance from begining point
		# compute angle by vector
		Ax <- (x[i,j] - x[i-1,j])
		Bx <- (x[i+1,j] - x[i-1,j])
		Ay <- (y[i,j] - y[i-1,j])
		By <- (y[i+1,j] - y[i-1,j])
		hugo <- (Ax * By - Ay * Bx)/abs(Ax * By - Ay * Bx)
		cos <- (Ax * Bx + Ay * By) / ((Ax^2 + Ay^2)*(Bx^2 + By^2))^0.5
		angle[i,j] <- acos(cos) * hugo
		
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

## ‰‚ÍŠ‚æ‚è‚à•ûŒü“]Š·Šp“x‚ª‘å‚«‚¢
par(pin = c(4,4))
r <- hist(abs(angle),breaks=seq(0,pi,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), pch=20, xlim=c(0,3.2),
	xlab="angle", ylab = "frequency")
par(new=T)
r <- hist(abs(angle[,(1:length(sex))[sex=="M"]]),breaks=seq(0,pi,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), col="blue", pch=20, xlim=c(0,3.2),
	xlab="angle", ylab = "frequency")
par(new=T)
r <- hist(abs(angle[,(1:length(sex))[sex=="F"]]),breaks=seq(0,pi,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), col="red", pch=20, xlim=c(0,3.2),
	xlab="angle", ylab = "frequency")


## “ú•t‚Æ•ûŒü“]Š·Šp“x
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

## ‰‚ÆŠ‚Ì‘¬“x
par(pin = c(4,4))
r <- hist(abs(speed),breaks=seq(0,27,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), pch=20, xlim=c(0,27),
	xlab="speed", ylab = "frequency")
par(new=T)
r <- hist(abs(speed[,(1:length(sex))[sex=="M"]]),breaks=seq(0,27,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), col="blue", pch=20, xlim=c(0,27),
	xlab="speed", ylab = "frequency")
par(new=T)
r <- hist(abs(speed[,(1:length(sex))[sex=="F"]]),breaks=seq(0,27,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), col="red", pch=20, xlim=c(0,27),
	xlab="speed", ylab = "frequency")


## “ú•t‚Æ‘¬“x
par(pin = c(4,4))
r <- hist(abs(speed[,(1:length(day))[day==0]]),breaks=seq(0,27,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), pch=20, xlim=c(0,27),
	xlab="speed", ylab = "frequency")
par(new=T)
r <- hist(abs(speed[,(1:length(day))[day==1]]),breaks=seq(0,27,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), col="red", pch=20, xlim=c(0,27),
	xlab="speed", ylab = "frequency")
par(new=T)
r <- hist(abs(speed[,(1:length(day))[day==2]]),breaks=seq(0,27,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), col="blue", pch=20, xlim=c(0,27),
	xlab="speed", ylab = "frequency")
par(new=T)
r <- hist(abs(speed[,(1:length(day))[day==3]]),breaks=seq(0,27,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), col="green", pch=20, xlim=c(0,27),
	xlab="speed", ylab = "frequency")


## ‚à‚Á‚Æ‚à—£‚ê‚½êŠ
mostaway <- apply(move, 2, max)	#@Še trajectory ‚Ì most away point
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
plot(angleperspeed~speedseq,pch=1,xlim=c(0,26),ylim=c(0,1.5),col=1)

## max speed
sex2<-sex
sex2[sex2=="M"]<-1
sex2[sex2=="F"]<-2
plot(apply(speed,2,max)~apply(angle,2,max,na.rm=T),col=sex2)



ind<-45
library(rgl)
plot3d(x[1:(length(steps)-1),ind], y[1:(length(steps)-1),ind], time[,ind])


