###### ECON 706 - Problem Set 1 - Luduvice, Pietrosanti and Silverstein #######

# Libraries

library(vars)
library(tseries)
library(forecast)

# Data

start <- read.csv("HOUST.csv", header = TRUE)
startts <- ts(start["VALUE"], start = c(1968, 1), end = c(2014, 12), freq = 12)
starttsf <- window(startts, start = c(1968,1), end = c(2014,8))

comp <- read.csv("COMPUTSA.csv", header = TRUE)
compts <- ts(comp["VALUE"], start = c(1968, 1), end = c(2014, 12), freq = 12)
comptsf <- window(compts, start = c(1968,1), end = c(2014,8))

summary(startts)
summary(compts)

###### Univariate Analysis ######

#### Autocorrelations and Unit Root Detection ####

acfcompf=acf(comptsf, lag.max=100, main="ACF Completions");
pacfcompf=pacf(comptsf, lag.max=100, main="PACF Completions")

acfstartf=acf(starttsf, lag.max=100, main="ACF Starts")
pacfstartf=pacf(starttsf, lag.max=100, main="PACF Starts")

# Taking first difference of both series

difcomptsf=diff(comptsf) 
difstarttsf=diff(starttsf)

acfdifcompf=acf(difcomptsf, lag.max=100, main="ACF Difference in Completions")
acfdifstartf=acf(difstarttsf, lag.max=100, main="ACF Difference in Starts")

pacfdifcompf=pacf(difcomptsf, lag.max=100, main="PACF Difference in Completions")
pacfdifstartf=pacf(difstarttsf, lag.max=100, main="PACF Difference in Starts")

# Performing Augmented Dickey-Fuller unit root test on original series

test1=ur.df(comptsf, type=c("none"), selectlags = c("AIC"))
summary(test1)
test2=ur.df(comptsf, type=c("none"), selectlags = c("BIC")) 
summary(test2)
test3=ur.df(comptsf, type=c("drift"), selectlags = c("AIC"))
summary(test3)
test4=ur.df(comptsf, type=c("drift"), selectlags = c("BIC"))
summary(test4)
test5=ur.df(starttsf, type=c("none"), selectlags = c("AIC"))
summary(test5)
test6=ur.df(starttsf, type=c("none"), selectlags = c("BIC"))
summary(test6)
test7=ur.df(starttsf, type=c("drift"), selectlags = c("AIC"))
summary(test7)
test8=ur.df(starttsf, type=c("drift"), selectlags = c("BIC"))
summary(test8)

# Performing Augmented Dickey-Fuller unit root test on diff series

test9=ur.df(difstarttsf, type=c("none"), selectlags = c("AIC")) 
summary(test9)
test9b=ur.df(difstarttsf, type=c("none"), selectlags = c("BIC"))
summary(test10)
test10=ur.df(difcomptsf, type=c("none"), selectlags = c("AIC"))
summary(test10)
test10b=ur.df(difcomptsf, type=c("none"), selectlags = c("BIC"))
summary(test10b)
test11=ur.df(difstarttsf, type=c("drift"), selectlags = c("AIC"))
summary(test11)
test11b=ur.df(difstarttsf, type=c("drift"), selectlags = c("BIC"))
summary(test11b)
test12=ur.df(difcomptsf, type=c("drift"), selectlags = c("AIC"))
summary(test12)
test12b=ur.df(difcomptsf, type=c("drift"), selectlags = c("BIC"))
summary(test12b)

#### Model Selection and Residual Analysis ####

# Selecting the model with the SIC(BIC) criterion

p=4;
q=4;
fitcompBIC=matrix(,p,q)
for (i in 1:p){
  for (j in 1:q){
    fitcomp=arima(comptsf, order=c(i-1,1,j-1))
    fitcompBIC[i,j]=AIC(fitcomp,k=log(length(comptsf)))
  } 
}

fitstartBIC=matrix(,p,q)
for (i in 1:p){
  for (j in 1:q){
    fitstart=arima(starttsf, order=c(i-1,1,j-1))
    fitstartBIC[i,j]=AIC(fitstart,k=log(length(starttsf)))
  } 
}


# Estimation of ARIMA for both series

bfitcomp12 = arima(comptsf, order=c(1,1,2))
bfitstart01=arima(starttsf, order=c(0,1,1))
bfitstart12 = arima(starttsf, order=c(1,1,2))

# Extracting residuals and evaluating ACF

rescomp12 = resid(bfitcomp12)
resstart01 = resid(bfitstart01)
resstart12 = resid(bfitstart12)

acfrescomp12 = acf(rescomp12)
acfresstart01 = acf(resstart01)
acfresstart12 = acf(resstart12)

# Testing serial correlation of residuals with the Box-Pierce test

p=4;
q=4;
fitcompLBOX=matrix(,p,q)
for (i in 1:p){
  for (j in 1:q){
    fitcomp=arima(comptsf, order=c(i-1,1,j-1))
    rescomp=residuals(fitcomp)
    btestcomp=Box.test(rescomp, lag = 23, type = c("Box-Pierce"), fitdf = i+j)
    fitcompLBOX[i,j]=btestcomp$p.value
  } 
}

fitstartLBOX=matrix(,p,q)
for (i in 1:p){
  for (j in 1:q){
    fitstart=arima(starttsf, order=c(i-1,1,j-1))
    resstart=residuals(fitstart)
    bteststart=Box.test(resstart, lag = 23, type = c("Box-Pierce"), fitdf = i+j)
    fitstartLBOX[i,j]=bteststart$p.value
  } 
}

###### Multivariate Analysis ######

colltsf=cbind(starttsf,comptsf)
colltsf1=cbind(comptsf,starttsf)

#### Cross-correlations and Identification ####

# Cross-correlations and Partial Auto-correlations

ccf(starttsf[,1], comptsf[,1], lag.max=100, main="Cross-correlation Start with Completion")
ccf(comptsf[,1], starttsf[,1], lag.max=100, main="Cross-correlation Completion with Start")

# Lag selection

VARselect(colltsf, lag.max=50, type="const")
VARselect(colltsf1, lag.max=50, type="const")

# Estimation

var3 = VAR(colltsf, p=3, type = "const")
summary(var3)
var14 = VAR(colltsf, p=14, type ="const")
summary(var14)

var3inv = VAR(colltsf1, p=3, type="const")
summary(var3inv)
var14inv = VAR(colltsf1, p=14, type="const")
summary(var14inv)


# Residual analysis

serial.test(var3, type="PT.asymptotic", lags.pt=23)
serial.test(var14, type="PT.asymptotic", lags.pt=23)


# Cointegration test

coint1 = ca.jo(colltsf, type="trace", K=14)
summary(coint1)

#### Granger causality and IRFs ####

bfitcol = VAR(colltsf, p=14, type ="const") # Just renaming
bfitcolinv = VAR(colltsf1, p=14, type="const")

caus1 = causality(bfitcol)
caus1
caus2 = causality(bfitcol, cause= "comptsf")
caus2


# IRFs

irf1 = irf(bfitcol, n.ahead=50)
plot(irf1, main)
irfinv= irf(bfitcolinv, n.ahead=50)
plot(irfinv)



###### Forecasts ######


### Univariate ###

bfitcomp12 = arima(comptsf, order=c(1,1,2))
bfitstart12 = arima(starttsf, order=c(1,1,2))

forcomp12=predict(bfitcomp12, 4);
forstart12=predict(bfitstart12, 4);


plot(startts, xlab="t", ylab="House Start Index",main="Overlay Forecasts & Actuals,  ARIMA(1,1,2)")
lines(forstart12$pred, col="blue")
lines(forstart12$pred+2*forstart12$se, col="red")
lines(forstart12$pred-2*forstart12$se, col="red")


plot(compts, xlab="t", ylab="House Completion Index",main="Overlay Forecasts & Actuals,  ARIMA(1,1,2)")
lines(forcomp12$pred, col="blue")
lines(forcomp12$pred+2*forcomp12$se, col="red")
lines(forcomp12$pred-2*forcomp12$se, col="red")


### Multivariate ###


varfor=predict(bfitcol, n.ahead=4);

varprstart=ts(varfor$fcst$start[,1], start=c(2014, 9), end=c(2014, 12), frequency=12);
varCIstart=ts(varfor$fcst$start[,4], start=c(2014, 9), end=c(2014, 12), frequency=12);
varprcomp=ts(varfor$fcst$comp[,1], start=c(2014, 9), end=c(2014, 12), frequency=12);
varCIcomp=ts(varfor$fcst$comp[,4], start=c(2014, 9), end=c(2014, 12), frequency=12);

par(mfrow=c(1,1)) 

plot(startts, xlab="t", ylab="House Start Index",main="Overlay Forecasts & Actuals, VAR(14)")
lines(varprstart, col="blue")
lines(varprstart+varCIstart, col="red")
lines(varprstart-varCIstart, col="red")

plot(compts, xlab="t", ylab="House Completion Index",main="Overlay Forecasts & Actuals,  VAR(14)")
lines(varprcomp, col="blue")
lines(varprcomp+varCIcomp, col="red")
lines(varprcomp-varCIcomp, col="red")

startts_small = window(startts, start=2013);
compts_small = window(compts, start=2013);

plot(startts_small, ylim=c(700,1400), xlab="t", ylab="Housing Starts (Thousands)",main="Overlay Forecasts & Actuals,  ARIMA(1,1,2)")
lines(forstart12$pred, col="blue")
lines(forstart12$pred+2*forstart12$se, col="red")
lines(forstart12$pred-2*forstart12$se, col="red")

plot(compts_small, ylim=c(700,1400), xlab="t", ylab="Housing Completions (Thousands)",main="Overlay Forecasts & Actuals,  ARIMA(1,1,2)")
lines(forcomp12$pred, col="blue")
lines(forcomp12$pred+2*forcomp12$se, col="red")
lines(forcomp12$pred-2*forcomp12$se, col="red")

plot(startts_small, ylim=c(700,1400), xlab="t", ylab="Housing Starts (Thousands)",main="Overlay Forecasts & Actuals, VAR(14)")
lines(varprstart, col="blue")
lines(varprstart+varCIstart, col="red")
lines(varprstart-varCIstart, col="red")

plot(compts_small, ylim=c(700,1400), xlab="t", ylab="Housing Completions (Thousands)", main="Overlay Forecasts & Actuals, VAR(14)")
lines(varprcomp, col="blue")
lines(varprcomp+varCIcomp, col="red")
lines(varprcomp-varCIcomp, col="red")
