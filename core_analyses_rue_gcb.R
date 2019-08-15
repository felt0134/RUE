#analysis of the relationship between rain use efficiency and growing season precipitation####
#testing models of pue 
pue<-read.csv(file.choose(),header=TRUE) #the pue-precip relationship
pue.2<-pue[-6,] #remove 1989 because seond year of a drought 

## basic summary stats ##

#look at average experimental soil moisture
lm.vwc<-lm(VWC~GSP,data=exp)
summary(lm.vwc)

#mean ANPP
anpp.rue<-aggregate(ANPP~Variability,mean,data=pue.2)
# 24% reduction in ANPP 
(512.12-391.27)/512.12

## mean GSP ##
mm.rue<-aggregate(GSP~Variability,mean,data=pue.2)
# 6.5 % difference in average growing season precipitation
(471.5333 - 442.7355)/442.7355

#subset to create and test model fits for the two rainfall patterns
exp<-subset(pue.2,Variability=="Reduced")
obs<-subset(pue.2,Variability=="Ambient")

#look at and compare anpp-gsp linear models
anpp.exp<-lm(ANPP~GSP,data=exp)
summary(anpp.exp)

anpp.obs<-lm(ANPP~GSP,data=obs)
summary(anpp.obs)

#experimental PUE-GSP models
lm.exp<-lm(PUE~GSP,data=exp)
plot(PUE~GSP,data=exp)
poly.exp<-lm(PUE~GSP+I(GSP^2),exp)

#compare models with analysis of variance
anova(poly.exp,lm.exp)
#nonlinear a better fit.

#observational
lm.obs<-lm(PUE~GSP,data=obs)
poly.obs<-lm(PUE~GSP + I(GSP^2),obs)
AIC(lm.obs,poly.obs) 
29.29 - 27.5 #AIC difference
anova(poly.obs,lm.obs)
#linear a better fit

#driest year prediction
driest<-data.frame(GSP=182)
predict(poly.exp,driest) #2.08
predict(lm.obs,driest) #1.2
(2.08-1.2)/2.08

#comparing years
library(nlme)

#looking at 2015 and 2016 RUE-GSP individually (data aavailable upon request)
#data 
exp.both.years<-read.csv(file.choose())
head(exp.both.years)

#2015
fifteen<-subset(exp.both.years,Variability.2=="Even.2015")
plot(PUE~GSP,data=fifteen)
lm.exp.2015<-lm(PUE~GSP,data=fifteen)
plot(PUE~GSP,data=fifteen)
poly.exp.2015<-lm(PUE~GSP+I(GSP^2),fifteen)
AIC(lm.exp.2015,poly.exp.2015) 
anova(lm.exp.2015,poly.exp.2015) #nonlinear far better fit

#2016
sixteen<-subset(exp.both.years,Variability.2=="Even.2016")
lm.exp.2016<-lm(PUE~GSP,data=sixteen)
plot(PUE~GSP,data=sixteen)
poly.exp.2016<-lm(PUE~GSP+I(GSP^2),sixteen)
AIC(lm.exp.2016,poly.exp.2016) 
anova(lm.exp.2016,poly.exp.2016) #nonlinear far better

#analysis of variance for even versus variable rainfall patterns
library(car)
altered_var<-subset(pue.2,Variability=="Reduced")
natural_var<-subset(pue.2,Variability=="Ambient")

#assumptions samples are drawn from a population with an approximately normal distribution
shapiro.test(altered_var$ANPP) #normal
qqPlot(altered_var$ANPP) #approximately normal, some skew at the first tail
shapiro.test(sqrt(natural_var$ANPP)) #normal
qqPlot(natural_var$ANPP) #approximatel ynormal, some skew towards end tail

#homogenity of variances
leveneTest(ANPP~Variability,data=pue.2) #homegenity of variance assumption not refuted

#analysis of variance
anova_model <- aov(ANPP~Variability,data=pue.2)
summary(anova_model)

######model selection for dry years#######
library(nlme)
library(MASS)

head(pue.summary)
pue.summary<-read.csv(file.choose(),header=TRUE) #summary stats for dry years
pue.drought<-subset(pue.summary, gsp.mm.may.august< 500)
pue.drought.2<-pue.drought[-5,] #remove 1989 from analysis
head(pue.drought.2)

#assessing multi-colinearity
#initial model
vif.test.1<-lm(RUE ~ gsp.mm.may.august  + max.event  + cv.event.size + mean.event  + X..rain.after.30 + CDD,data=pue.drought.2)
vif(vif.test.1) #colinearity among predictor variables, remove CV
vif.test.2<-lm(RUE ~ gsp.mm.may.august  + max.event + mean.event  + X..rain.after.30 + CDD,data=pue.drought.2)
vif(vif.test.2) #still high colinearity, remove next highest: mean event size

#After removing two most colinear variables (mean and cv of event size)
vif.test.3<-lm(RUE ~ gsp.mm.may.august + max.event + X..rain.after.30 + CDD,data=pue.drought.2)
vif(vif.test.3) #multicolinearity reduced below 5 for all variables, use this for model selection

#model selection after reducing muticolinearity including total precip
stepAIC(lm(RUE ~ gsp.mm.may.august + max.event + X..rain.after.30 + CDD,data=pue.drought.2,na.action = na.exclude))

#multicolinear model
stepAIC(lm(RUE ~ gsp.mm.may.august + max.event  + cv.event.size + mean.event+ X..rain.after.30 + CDD,data=pue.summary.2,na.action = na.exclude))

#testing for colinearity in the selected model
vif.test.selected<-lm(RUE ~ gsp.mm.may.august + max.event,data=pue.drought.2) 
vif(vif.test.selected)
#low colinearity
summary(vif.test.selected)

#partial correllation for GSP controlling for max.event
mm1 = lm(gsp.mm.may.august ~ max.event,data=pue.drought.2)
res1 = mm1$residuals
mm2 = lm(RUE~max.event,data=pue.drought.2)
res2 = mm2$residuals
cor(res1,res2)

#partial correlation of max.event controlling for GSP
mm3 = lm(max.event ~ gsp.mm.may.august,data=pue.drought.2)
res3 = mm3$residuals
mm4 = lm(RUE~gsp.mm.may.august,data=pue.drought.2)
res4 = mm4$residuals
cor(res3,res4)


#####correlations#######
library(Hmisc)

#load dataset
precip.exp.attributes<-read.csv(file.choose(),header=TRUE) #comparing CV, mean event size, and GSP among datasets

#correlation coefficients: experimental event size and cv of event size
precip.attributes.2<- precip.exp.attributes[-c(15),] #remove 1989
head(precip.exp.attributes.2)
exp.pue<-subset(precip.attributes.2,dataset=="Experimental")

#percent differnce of CV of event size in ambient versus experimental regimes
mean(exp.pue$CV) #6.26
mean(full.dataset.2$CV) #90.3
(90.3 - 6.26)/90.3 #93% difference in event size variability

#pearson correlation coefficients
exp.pue.2 <- exp.pue[, c(3,4)] 
exp.corr.2<-rcorr(as.matrix(exp.pue.2),type="pearson") #make correlation matrix
exp.corr.2$r #get coefficients
exp.corr.2$P #get p-value

#correlation coefficients: experimental event size and GSP
head(precip.attributes.2)
exp.pue.3 <- exp.pue[, c(4,5)] 
exp.corr.3<-rcorr(as.matrix(exp.pue.3),type="pearson") #make correlation matrix
exp.corr.3$r #get coefficients
exp.corr.3$P

#correlations with RUE during dry years
head(pue.drought.2)
pue.drought.3 <- pue.drought.2[, c(2,3,4,5,6,9,10)] #isolate variables of interest for the matric
# print the first 6 rows
res<-rcorr(as.matrix(pue.drought.3),type="pearson") #make correlation matrix
table1.corr<-res$r #get coefficients
res$P #get p values

#make into table
write.csv(table1.corr,file="corr.coeff.csv")
plot(mean.event~cv.event.size,data=pue.summary.2)
