# R script for Oufiero et al. 2024
#Oufiero, C.E., L. Garikipati, E. McMillan, M.K. Sullivan and R. Turnbaugh. Modulation of prey capture kinematics in relation to prey distance helps predict success. Journal of Experimental Biology. 
#Analyses for dataset on S.lineola females (N=11) success/failure



data<-read.csv("Mantis 3D filming log.csv",nrows=240,stringsAsFactors=TRUE)
data$Tibia_AV_min_abs<-abs(data$Tibia_AV_min)
data$Tibia_LV_min_abs<-abs(data$Tibia_LV_min)
data<-data[,-c(7)]#remove empty column
data<-na.omit(data)#remove missing trials

data$suc<-as.character(data$success)
data$suc[data$suc=="yes"]<-1
data$suc[data$suc=="no"]<-0
data$suc<-as.numeric(data$suc)

data<-droplevels(data)
library(nlme)
library(lme4)
library(lmerTest)
library(ggplot2)
library(MuMIn)
attach(data)
library(sjPlot)
library(caret)
library(lavaan)
library(vegan)
library(plotrix)


#S.lineola 11 T20 values off, need to remove from dataset (body size too small, coxa way too long)
r<-which(data$Body_length<0.1)
data<-data[-r,]

#set up color vector

col<-hcl.colors(11,palette = "Set 2",alpha=0.5)

data$icol[data$Mantis==1] <- col[1]
data$icol[data$Mantis==2] <- col[2]
data$icol[data$Mantis==3] <- col[3]
data$icol[data$Mantis==4] <- col[4]
data$icol[data$Mantis==5] <- col[5]
data$icol[data$Mantis==6] <- col[6]
data$icol[data$Mantis==7] <- col[7]
data$icol[data$Mantis==8] <- col[8]
data$icol[data$Mantis==9] <- col[9]
data$icol[data$Mantis==11] <- col[10]
data$icol[data$Mantis==12] <- col[11]

data$shape[data$success=="yes"] <- 21
data$shape[data$success=="no"] <- 22


table(data$success)#count of success and failure

#split these into success and failure?
S<-(data[data$success=="yes",])
N<-(data[data$success=="no",])


#convert rcolors to rgb in plotting
scol=col2rgb("darkgoldenrod")
scol=scol/255
fcol<-col2rgb("darkgreen")
fcol<-fcol/255
#Figure 2

png(file="Prey Position histo.png", height=8, width=16,units="in",bg="white",res=300)
par(mfrow=c(1,2),mar=c(5.1,6,4.1,6.5))
hist(S$PP_dist_cm_st, col=rgb(scol[1],scol[2],scol[3],alpha=0.5),xlab="Prey Distance (cm)", main= "", xlim=c(0,6.5), ylim=c(0,40))
hist(N$PP_dist_cm_st, col=rgb(fcol[1],fcol[2],fcol[3],alpha=0.5),xlab="Prey Distance (cm)", main= "", add=TRUE)
ms<-mean(S$PP_dist_cm_st)
#abline=(v=ms)#,rgb(0.2,1,0.4,alpha=0.5),lwd=3)
abline(v=2.5, lwd=3)
abline(v=ms, lwd=3,col=rgb(scol[1],scol[2],scol[3]))
mf<-mean(N$PP_dist_cm_st)
abline(v=mf, lwd=3,col=rgb(fcol[1],fcol[2],fcol[3]))
segments(0,40,0.5,40)
text(0.5,40,"Mean from virtual studies",cex=0.7,pos=4)
segments(0,38,0.5,38,col=rgb(scol[1],scol[2],scol[3]))
text(0.5,38,"Mean from succesful",cex=0.7,pos=4,col=rgb(scol[1],scol[2],scol[3]))
segments(0,36,0.5,36,col=rgb(fcol[1],fcol[2],fcol[3]))
text(0.5,36,"Mean from unsuccesful",cex=0.7,pos=4,col=rgb(fcol[1],fcol[2],fcol[3]))

hist(S$PP_angle_st, col=rgb(scol[1],scol[2],scol[3],alpha=0.5),xlab="Prey Angle (degrees)", main= "", xlim=c(50,180), ylim=c(0,25))
hist(N$PP_angle_st, col=rgb(fcol[1],fcol[2],fcol[3],alpha=0.5),xlab="Prey Angle (degrees)", main= "", add=TRUE)
ms<-mean(S$PP_angle_st)
#abline=(v=ms)#,rgb(0.2,1,0.4,alpha=0.5),lwd=3)
#abline(v=2.5, lwd=3)
abline(v=ms, lwd=3,col=rgb(scol[1],scol[2],scol[3]))
mf<-mean(N$PP_angle_st)
abline(v=mf, lwd=3,col=rgb(fcol[1],fcol[2],fcol[3]))
#segments(0,40,0.5,40)
#text(0.5,40,"Mean from virtual studies",cex=0.7,pos=4)
segments(60,24,70,24,col=rgb(scol[1],scol[2],scol[3]))
text(70,24,"Mean from succesful",cex=0.7,pos=4,col=rgb(scol[1],scol[2],scol[3]))
segments(60,22,70,22,col=rgb(fcol[1],fcol[2],fcol[3]))
text(70,22,"Mean from unsuccesful",cex=0.7,pos=4,col=rgb(fcol[1],fcol[2],fcol[3]))
dev.off()


#ANOVA for success and failure in realtion to prey distance and angle
m<-aov(data$PP_dist_cm_st~data$success)
summary(m)

ma<-aov(data$PP_angle_st~data$success)
summary(ma)

#PCA of kinematic data
PCA.dat<-cbind(data[,c(17,20,27,28,31,32,35,87,45,46,65,66,67,23,24)])
PCA<-princomp(PCA.dat,cor=TRUE)
summary(PCA)

PCA.loadings<-PCA.loadings<-PCA$loadings
#write.csv(PCA.loadings, file="PCA.loadings.csv")
biplot(PCA)
bstick(PCA, tot.var=15)#note tot.var is total number of colums above 10 landmarks with x,y coordinate = 20 variables
screeplot(PCA, bstick = TRUE, type = "barplot",
         npcs = 15)
data$PC1<-PCA$scores[,1]         
data$PC2<-PCA$scores[,2]
data$PC3<-PCA$scores[,3]
data$PC4<-PCA$scores[,4]

# model predicting success with covariates
mod<-glm(success ~ PC1+PC2+PC3+PC4+mealy.size..g.+Body_length+Trial,data=data, family = binomial)
summary(mod)

# model predicting success with covariates
mod<-glm(success ~ PC1+PC2+PC3+PC4,data=data, family = binomial)
summary(mod)#no signifacnt covariates, leaving them out


#mixed model, singular with and without covariates, not reported
log.mod <- glmer(success ~ PC1+PC2+PC3+PC4+mealy.size..g.+Body_length+Trial+(1|Mantis),data=data, family = binomial,control = glmerControl(optimizer = "bobyqa"))
anova(log.mod)#singlular, don't need individual effect

Predicted_data <- data.frame(PC1 =seq(
  min(data$PC1), max(data$PC1 ),len=500),PC2 =seq(
  min(data$PC2), max(data$PC2 ),len=500),PC3 =seq(
  min(data$PC3), max(data$PC3 ),len=500),PC4 =seq(
  min(data$PC4), max(data$PC4 ),len=500))

Predicted_data$suc = predict(
  mod, Predicted_data, type="response")

#Figure 3
png(file="Predicted success PC1.png", height=8, width=8,units="in",res=300,bg="transparent")
par(mar=c(5.1,6,4.1,6.5))
plot(suc ~ PC1, data=data, pch=data$shape,bg=data$icol,xlab="PC1 22.4%", ylab="Prey Capture Success", yaxt='n',cex=1.5,cex.lab=1.5)
#mtext("No",side=2,at=0,line=2,cex=1.5,las=2)
#mtext("Yes",side=2,at=1,line=2, cex=1.5,las=2)
lines(suc ~ PC1, Predicted_data, lwd=2, col="black")
mtext("No",side=2,at=0,line=1,cex=1.5,las=2)
mtext("Yes",side=2,at=1,line=1, cex=1.5,las=2)
axis(side = 4, at = pretty((data$suc)))
mtext("Predicted Probablity Success", side = 4, line = 3)
mtext("Farther, expansive, fast", side = 1, at = 3, line = 2, cex=0.7)
mtext("prey capture strikes", side = 1, at = 3, line = 2.5, cex=0.7)
mtext("Closer, reduced, slow", side = 1, at = -7, line = 2, cex=0.7)
mtext("prey capture strikes", side = 1, at = -7, line = 2.5, cex=0.7)
dev.off()


#color this by predicted success
data$pred.suc<-predict(mod, type="response")
pred.suc.rnd<-round(100*(data$pred.suc-min(data$pred.suc))/(max(data$pred.suc)-min(data$pred.suc)), digits=3)
 pallete<-heat.colors(100)
data$pred.col<-(pallete[pred.suc.rnd+1])


s.col<-sort(data$pred.col, decreasing=TRUE)

#Figure 4
png(file="PC1 PC2 pred suc.png", height=8, width=8,units="in",res=300,bg="transparent")
par(mar=c(5.1,6,4.1,6.5))
plot(data$PC1,data$PC2,pch=data$shape,bg=data$pred.col, xlab="PC1 22.85%", ylab="PC2 18.70%")

gradient.rect(-8,7.5,-5, 8, col=smoothColors(s.col),nslices=100, border=TRUE, gradient="x") 
text(-8,7, "Increased", pos=3, cex=0.7)
text(-5, 7, "Decreased", pos=3,cex=0.7)
text(-6.5, 6.5, "Probability Success", pos=3,cex=0.8)
mtext("Farther, expansive, fast", side = 1, at = 3, line = 2, cex=0.7)
mtext("prey capture strikes", side = 1, at = 3, line = 2.5, cex=0.7)
mtext("Closer, reduced, slow", side = 1, at = -7, line = 2, cex=0.7)
mtext("prey capture strikes", side = 1, at = -7, line = 2.5, cex=0.7)

mtext("Increased prey angle,", side = 2, at = 7, line = 3, cex=0.7)
mtext("foreleg lateral displacement", side = 2, at = 7, line = 2.5, cex=0.7)
mtext("body movement", side = 2, at = 7, line = 2, cex=0.7)

mtext("Faster foreleg", side = 2, at = -2, line = 2.5, cex=0.7)
mtext("angular velocities", side = 2, at = -2, line = 2, cex=0.7)
dev.off()


#####just prey distance predicting success
mod<-glm(suc ~ PP_dist_cm_st,data=data, family = binomial)
summary(mod)

Predicted_data <- data.frame(PP_dist_cm_st =seq(
  min(data$PP_dist_cm_st), max(data$PP_dist_cm_st ),len=500))

Predicted_data$suc = predict(
  mod, Predicted_data, type="response")


#FIgure 5
png(file="Predicted success prey dist.png", height=8, width=8,units="in",res=300,bg="transparent")
par(mar=c(5.1,6,4.1,6.5))
plot(suc ~ PP_dist_cm_at, data=data, pch=data$shape,bg=data$icol,xlab="Prey Distance (cm)", ylab="Prey Capture Success", yaxt='n',cex=1.5,cex.lab=1.5)
lines(suc ~ PP_dist_cm_st, Predicted_data, lwd=2, col="black")
mtext("No",side=2,at=0,line=1,cex=1.5,las=2)
mtext("Yes",side=2,at=1,line=1, cex=1.5,las=2)
axis(side = 4, at = pretty((data$suc)))
mtext("Predicted Probablity Success", side = 4, line = 3)

segments(2.5,-1,2.5,0.5186913,lty="dashed",col="darkorange3",lwd=2)
segments(2.5,0.5186913,6.4,0.5186913,lty="dashed",col="darkorange3",lwd=2)
mpd<-mean(data$PP_dist_cm_st)
segments(mpd,-1,mpd,0.4231015,lty="dotted",col="darkorchid4",lwd=2)
segments(mpd,0.4231015,6.4,0.4231015,lty="dotted",col="darkorchid4",lwd=2)

text(2.4,0.2,"Virtual Preferred Distance",cex=0.7,col="darkorange3",srt=90)
text(mpd-.1,0.2,"Preferred Distance This Study",cex=0.7,col="darkorchid4",srt=90)

dev.off()


