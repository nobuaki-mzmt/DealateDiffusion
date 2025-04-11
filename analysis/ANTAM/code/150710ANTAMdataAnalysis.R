## ANTAM trajectory data analysis
## 15_7_06

## packages
library(car)
library(lme4)

Folder.path <- getwd()
f.namesplace <- list.files(Folder.path, pattern=".csv",full.names=T)

## analysis condition
FPS <- 2	## FPS using in analysis
secper <- 1/FPS	## 1コマあたりのsec
steps <- FPS * 60 *20 +1	## 20分のコマ数

## reading parameters
x = y = sec = time <-  matrix(0, ncol=length(f.namesplace), steps)
name = colony = sex = rep = day = totaldis = totaltime <- rep(0,length(f.namesplace))
angle = dis = move = speed = gostop = reorient = acceleration <- matrix(0, ncol=length(f.namesplace), steps)

	for(j in 1:length(f.namesplace)){
		
	d <- read.csv(f.namesplace[j],sep=",")[,c(1,2,4)]
	namebegin <- regexpr("mouselog_",f.namesplace[j]) + 9
	nameend <- regexpr(".csv",f.namesplace[j]) - 1

	name[j] <- substring(f.namesplace[j],namebegin, nameend)
	sex[j] <- substring(f.namesplace[j],nameend-2, nameend-2)
	rep[j] <- substring(f.namesplace[j],nameend-1, nameend-1)
	data <- (d[d[,3]>60*10,])	## tracking開始後10分以降のデータを使う
	data[,1] <- data[,1] - data[1,1]
	data[,2] <- data[,2] - data[1,2]
	data[,3] <- data[,3] - data[1,3]
	x[1,j] <- data[1,1]
	y[1,j] <- data[1,2]
	rownum <- 2
	i <- 2

	## distance: 直前の位置から、現在位置までの距離 （1列目は0）
	## speed: その際にかかった時間(sec)を割ったもの （1列目は0）
	## acceleration: その速度になるための速度の変化量
	## move: 初期地から離れた距離
	while(rownum < steps+1){
		if((abs(secper - (data[i,3]-time[rownum-1, j]))) < (abs(secper - (data[i+1,3]-time[rownum-1, j])))){	## もっとも1コマ分に近くなるとき
			x[rownum,j] <- data[i,1]
			y[rownum,j] <- data[i,2]
			time[rownum,j] <- data[i,3]
			sec[rownum,j] <- time[rownum,j] - time[rownum-1,j]
			dis[rownum,j] <- ((x[rownum,j]-x[rownum-1,j])^2+(y[rownum,j]-y[rownum-1,j])^2)^0.5	#
			speed[rownum,j] <- dis[rownum,j] / sec[rownum,j]						# compute moving speed
			acceleration[rownum,j] <- speed[rownum,j] - speed[rownum-1,j]
			move[rownum,j] <- ((x[rownum,j]-x[1,j])^2+(y[rownum,j]-y[1,j])^2)^0.5		# compute distance from begining point
			rownum <- rownum + 1
		}
		i <- i + 1	
	}
	## angle: 直前の位置と直後の位置から、現在位置での回転角を算出
	## 
	for(i in 2: (steps-1)){
		# compute angle by vector
		Ax <- (x[i,j] - x[i-1,j])
		Bx <- (x[i+1,j] - x[i,j])
		Ay <- (y[i,j] - y[i-1,j])
		By <- (y[i+1,j] - y[i,j])
		hugo <- (Ax * By - Ay * Bx)/abs(Ax * By - Ay * Bx)
		cos <- (Ax * Bx + Ay * By) / ((Ax^2 + Ay^2)*(Bx^2 + By^2))^0.5
		angle[i,j] <- acos(cos) * hugo
		if(x[i-1,j]==x[i+1,j] | y[i-1,j]==y[i+1,j] | 
			Ax*By == Ay*Bx){
			angle[i,j] <- 0
		}
	}
}





## ANTAM data analysis
## for Plot

## 方向転換角度
par(pin = c(4,4))
r <- hist(abs(angle),breaks=seq(0,pi,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), pch=20, xlim=c(0,3.2),
	xlab="angle (radian)", ylab = "frequency", las=1)
par(new=T)
r <- hist(abs(angle[,(1:length(sex))[sex=="M"]][angle[,(1:length(sex))[sex=="M"]] <250]),breaks=seq(0,pi,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), col="blue", pch=20, xlim=c(0,3.2),
	xlab="angle (radian)", ylab = "frequency", las=1)
par(new=T)
r <- hist(abs(angle[,(1:length(sex))[sex=="F"]][angle[,(1:length(sex))[sex=="F"]] <250]),breaks=seq(0,pi,length=50),plot=F)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-6,1), col="red", pch=20, xlim=c(0,3.2),
	xlab="angle (radian)", ylab = "frequency", las=1)
legend(legend=c("total","male","female"),"topright",pch=20,col=c(1,4,2))
maleangle <- as.vector(angle[,1:5])
femaleangle <- as.vector(angle[,6:10])
wilcox.test(maleangle,femaleangle)


## ♂と♀の速度
par(pin = c(4,4))
r <- hist(abs(speed/51.747),plot=F,breaks=seq(0,3500/51.747,length=50))
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	 pch=20, xlim=c(0,50),ylim = c(10^-5,1),ann=F,xaxt="n",yaxt="n")
#	xlab="speed", ylab = "frequency")
par(new=T)
r <- hist(abs(speed[,(1:length(sex))[sex=="M"]]/51.747),plot=F,breaks=seq(0,3500/51.747,length=50))
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-5,1), col="blue", pch=20, xlim=c(0,50),ann=F,xaxt="n",yaxt="n")
	#xlab="speed", ylab = "frequency")
par(new=T)
r <- hist(abs(speed[,(1:length(sex))[sex=="F"]]/51.747),plot=F,breaks=seq(0,3500/51.747,length=50))
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-5,1), col="red", pch=20, xlim=c(0,50),
	xlab="speed (mm/second)", ylab = "frequency",las=1)
legend(legend=c("total","male","female"),"topright",pch=20,col=c(1,4,2))

## スケーリング合わせ
## 5174.7 で 10cm
## 51.747/1 mm
sex
matplot(x/215.8,y/215.8,pch=".",xlim=c(-1500,1500),ylim=c(-1500,1500),col=sex2)
turn <- angle
turn[abs(turn) < 1] <- 0
turn[abs(turn) >= 1] <- 1
plot(x[,1]/215.8,y[,1]/215.8,pch=".",col=turn[,1]+1,cex=(turn[,1]*2+1))
plot(c(0,x[,1][turn[,1]>0])/215.8,c(0,y[,1][turn[,1]>0])/215.8,type="l")
abline(h=0,v=0)
plot(c(0,x[,2][turn[,2]>0])/215.8,c(0,y[,2][turn[,2]>0])/215.8,type="l")
plot(c(0,x[,3][turn[,3]>0])/215.8,c(0,y[,3][turn[,3]>0])/215.8,type="l")
plot(c(0,x[,4][turn[,4]>0])/215.8,c(0,y[,4][turn[,4]>0])/215.8,type="l")
plot(c(0,x[,5][turn[,5]>0])/215.8,c(0,y[,5][turn[,5]>0])/215.8,type="l")
plot(c(0,x[,6][turn[,6]>0])/215.8,c(0,y[,6][turn[,6]>0])/215.8,type="l")
plot(c(0,x[,7][turn[,7]>0])/215.8,c(0,y[,7][turn[,7]>0])/215.8,type="l")
plot(c(0,x[,8][turn[,8]>0])/215.8,c(0,y[,8][turn[,8]>0])/215.8,type="l")
plot(c(0,x[,9][turn[,9]>0])/215.8,c(0,y[,9][turn[,9]>0])/215.8,type="l")
plot(c(0,x[,10][turn[,10]>0])/215.8,c(0,y[,10][turn[,10]>0])/215.8,type="l")

theta <- seq(-pi, pi, length=100)
plot(100000*cos(theta)/51.747, 100000*sin(theta)/51.747, type="l", asp=1,ann=F, xaxt="n",yaxt="n"
	,xlim=c(-300000/51.747,300000/51.747) ,ylim=c(-300000/51.747,300000/51.747))
par(new=T)
plot(150000*cos(theta)/51.747, 150000*sin(theta)/51.747, type="l", asp=1,xlim=c(-300000/51.747,300000/51.747)
		,ylim=c(-300000/51.747,300000/51.747),xlab="x",ylab="y",ann=F, xaxt="n",yaxt="n")
par(new=T)
plot(200000*cos(theta)/51.747, 200000*sin(theta)/51.747, type="l", asp=1,xlim=c(-300000/51.747,300000/51.747)
		,ylim=c(-300000/51.747,300000/51.747),xlab="x (mm)",ylab="y (mm)",las=1)

for(i in 1:10){
	par(new=T)
	Sys.sleep(1)
	plot(x[,i]/51.747,y[,i]/51.747,type="l",xlim=c(-300000/51.747,300000/51.747),ylim=c(-300000/51.747,300000/51.747),
		asp=1,col=as.numeric(sex2[i]),ann=F, xaxt="n",yaxt="n",las=1)
}
abline(h=0,v=0)
legend(legend=c("male","female"),"topright",pch="-",col=c(1,2))

## もっとも離れた場所
sex2<-sex
sex2[sex2=="M"]<-1
sex2[sex2=="F"]<-2
mostaway <- apply(move, 2, max)	#　各 trajectory の most away point
totaldis <- apply(dis, 2, sum)
plot(mostaway ~ totaldis,col=as.numeric(day)+1,pch=20)
plot(mostaway ~ totaldis,col=sex2)
par(pin = c(4,4))
r <- hist(abs(move)/51.747,breaks=seq(0,371000/51.747,length=50),plot=T)
plot(r$breaks[-1], r$counts/sum(r$counts), log='y',
	ylim = c(10^-4,1), pch=20, xlim=c(0,371000/51.747),
	xlab="distance from begining point", ylab = "frequency",ann=F,xaxt="n",yaxt="n")
par(new=T)
rM <- hist(abs(move[,(1:length(sex))[sex=="M"]]/51.747),breaks=seq(0,371000/51.747,length=50),plot=F)
plot(rM$breaks[-1], rM$counts/sum(rM$counts), log='y',
	ylim = c(10^-4,1), pch=20, xlim=c(0,371000/51.747),col="blue",
	xlab="distance from begining point", ylab = "frequency",ann=F,xaxt="n",yaxt="n")
par(new=T)
rF <- hist(abs(move[,(1:length(sex))[sex=="F"]]/51.747),breaks=seq(0,371000/51.747,length=50),plot=F)
plot(rF$breaks[-1], rF$counts/sum(rF$counts), log='y',
	ylim = c(10^-4,1), pch=20, xlim=c(0,371000/51.747),col="red",
	xlab="distance from begining point (mm)", ylab = "frequency",las=1)
legend(legend=c("total","male","female"),"topright",pch=20,col=c(1,4,2))

plot(r$breaks[-1], log10(r$counts/sum(r$counts)),
	ylim = c(-4,0), pch=20, xlim=c(0,371000/51.747),
	xlab="distance from begining point", ylab = "frequency",ann=F,xaxt="n",yaxt="n")
par(new=T)
plot(rM$breaks[-1],log10(rM$counts/sum(rM$counts)),col=4,
	ylim = c(-4,0), pch=20, xlim=c(0,371000/51.747),ann=F,xaxt="n",yaxt="n")
resM<-lm(log10(rM$counts/sum(rM$counts))~rM$breaks[-1])
abline(resM$coefficients[1],resM$coefficients[2],col=4)
par(new=T)
plot(rF$breaks[2:22],log10(rF$counts/sum(rM$counts))[1:21],col=2,
	ylim = c(-4,0), pch=20, xlim=c(0,371000/51.747),
	xlab="distance from begining point", ylab = "frequency",las=1)
resM<-lm(log10(rF$counts/sum(rF$counts))[1:21]~rM$breaks[2:22])
abline(resM$coefficients[1],resM$coefficients[2],col=2)
legend(legend=c("total","male","female"),"topright",pch=20,col=c(1,4,2))


length(rM$breaks[-1])
dbegin <- c(rM$breaks[-1],rF$breaks[-1])
freq <- log(c(rM$counts/sum(rM$counts),rF$counts/sum(rF$counts)))
sexdis <- rep(c("M","F"),each=49)
freq2 <- freq[!is.infinite(freq)]
r<-lm(freq2~dbegin[1:length(freq2)]*sexdis[1:length(freq2)])
Anova(r)
summary(r)
r$coefficients
abline(r$coefficients[2],r$coefficients[1])
abline
x<-seq(0,10,length=1000)
y<-x
plot(log(y)~x,pch=".")

## speed~angle 
angle2 <- angle#[,(1:length(sex))[sex=="F"]]
speed2 <- speed#[,(1:length(sex))[sex=="F"]]
dim(angle2) <- c(length(angle2),1)
dim(speed2) <- c(length(speed2),1)
speed3 <- round(speed2,-1)
#plot(abs(angle2)~speed3)
angleperspeed <- tapply(abs(angle2),speed3,mean,na.rm=TRUE)
speedseq <- seq(0,max(speed3),length=length(angleperspeed))/51.747
plot(angleperspeed~speedseq,pch=1,col=1,las=1
	,ylab="turn angle (radian)",xlab="speed (mm/second)")



plot(apply(speed,2,max)~apply(angle,2,max,na.rm=T),col=sex2)
(apply(speed,2,max)~sex)
apply(speed,2,mean)
apply(angle,2,mean)
apply(move,2,mean)
apply(move,2,max)



