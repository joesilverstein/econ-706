###### ECON 706 - Problem Set 2 - Luduvice, Pietrosanti and Silverstein #######

# Libraries

library(vars)
library(tseries)
library(forecast)
library(TSA)

source("lagWindow.R") # causes the parser to read the specified file
source("var.spec.R")
source("spectralWindow.R")
source("mySpecPgram.R")

# Data

start = read.csv("HOUST.csv", header = TRUE)
startTS = ts(start["VALUE"], start = c(1970, 1), end = c(2015, 1), freq = 12)

comp = read.csv("COMPUTSA.csv", header = TRUE)
compTS = ts(comp["VALUE"], start = c(1970, 1), end = c(2015, 1), freq = 12)

compNSA = read.csv("COMPUTNSA.csv", header = TRUE)
compNSATS = ts(compNSA["VALUE"], start = c(1970, 1), end = c(2015, 1), freq = 12)

###### Univariate Analysis ######

# Taking first difference of both series

diffStartTS = diff(startTS)
diffCompTS = diff(compTS) 

# Raw Periodograms:
spectrum(diffStartTS, log="no")
spectrum(diffCompTS, log="no")

# Estimate univariate spectral density

# Autoregressive spectral density estimation:
# No need to set bandwidth because the estimation is parametric
# Use ARIMA(1,1,2) as in the last homework
spec.ar(diffStartTS, log="no")
spec.ar(diffCompTS, log="no")
# It makes sense that AIC chooses AR(12) and AR(13) because in the previous homework it choose 12 lags for VAR (which only includes AR lags)

# Spectral Window (Smoothed Periodogram) Spectral Estimation

# First plot the raw spectrum for comparison
spectrum(diffStartTS) # spectral representation of raw time series data

# Try several different bandwidths and number of passes for the Daniell kernel and compare
# with the raw periodogram and the AR spectrum periodogram to
# determine which Daniell kernel to use

# First try different bandwidths at each pass
for (i in 3:10) {
  for (j in 3:10) {
    T=540
    bandwidthStart = (2*i+2*j+1)/T
    par(mfrow=c(1,2))
    spectrum(diffStartTS, log="no")  
    k = kernel("daniell",c(i,j))
    print(i)
    print(j)
    mySpecPgram(diffStartTS, k, taper=0, log="no", bandwidthDisplay=bandwidthStart)
    readline("Press <return> to continue")
  }
}
# Conclusion: Larger values of the time window m tend to oversmooth. 
# Bias-Variance tradeoff.

# Now try different numbers of passes, where each pass has the same bandwidth
spectralWindow(diffStartTS,seq(3,10,1),1)
spectralWindow(diffStartTS,seq(3,10,1),2)
spectralWindow(diffStartTS,seq(3,10,1),3)

### By inspection, we decide on using c(8,8) for starts.
k = kernel("daniell",c(8,8))
T=540
par(mfrow=c(1,1))
bandwidthStart = (4*8+1)/T # this is the correct bandwidth - not what is displayed
mySpecPgram(diffStartTS, k, taper=0, log="no", plot=TRUE, bandwidthDisplay=bandwidthStart)

# Now do the same for completions:

# First try different bandwidths at each pass
for (i in 1:12) {
  for (j in 1:8) {
    par(mfrow=c(2,2))
    spectrum(diffCompTS, log="no")
    spec.ar(diffCompTS, log="no")
    k = kernel("daniell",c(i,j))
    print(i)
    print(j)
    spec.pgram(diffCompTS, k, taper=0, log="no")
    out = readline("Press <return> to continue")
  }
}

# Now try different numbers of passes, where each pass has the same bandwidth
spectralWindow(diffCompTS,seq(3,10,1),1)
spectralWindow(diffCompTS,seq(3,10,1),2)
spectralWindow(diffCompTS,seq(3,10,1),3)

### By inspection, we decide on using c(7,7) for completions.
k = kernel("daniell",c(7,7))
spec.pgram(diffCompTS, k, taper=0, log="no")
bandwidthComp = (4*7+1)/T # this is the correct bandwidth - not what is displayed
mySpecPgram(diffStartTS, k, taper=0, log="no", plot=TRUE, bandwidthDisplay=bandwidthComp)

# Now perform lag window spectral estimation:

# Starts:

lagWindow(diffStartTS, seq(20,60,5), "", layout="group")
lagWindow(diffStartTS, seq(22,28,1), "", layout="group")

# By inspection, we decide to use M=25. The final graph is:
lagWindow(diffStartTS, 25, "Series: diffStartTS\nLag Window", layout="alone")

# Completions:

lagWindow(diffCompTS, seq(20,60,5), "", layout="group")
lagWindow(diffCompTS, seq(27,33,1), "", layout="group")

# By inspection, we decide to use M=30.
lagWindow(diffCompTS, 30, "Series: diffCompTS\nLag Window", layout="alone")

###### Multivariate Analysis ######

# VAR Spectral Density Estimation:

diffStartCompTS = cbind(diffStartTS, diffCompTS)

# estimation of VAR model using AIC to choose number of lags
X.p = ar(diffStartCompTS,aic=TRUE)

# make sure T is set to a number before running this, or else it will evaluate as TRUE
T = length(diffStartCompTS)
fr = rep(0,T/2)
for (j in 1:T/2) {
  fr[j] = (6/pi)*(2*pi*j/T) 
}

# Parametric spectral analysis
X.p.sp = var.spec(fr,X.p)

# Spectrum of Starts:
# Note: this throws a warning message because the values it's plotting are saved as complex numbers, but
# since they all have 0 imaginary component the coercion to real numbers is not problematic.
plot(fr,X.p.sp[,1,1],type="l",xlab=" ",ylab=" ")

# Spectrum of Completions:
# Note: this throws a warning message because the values it's plotting are saved as complex numbers, but
# since they all have 0 imaginary component the coercion to real numbers is not problematic.
plot(fr,X.p.sp[,2,2],type="l",xlab="Frequency",ylab="spectrum")

# Cross-spectrum
plot(fr,X.p.sp[,1,2],type="l",main="Cross Spectrum")

# Coherence
plot(fr,Mod(X.p.sp[,1,2])^2/( X.p.sp[,1,1]*X.p.sp[,2,2]), 
     type="l",main="Coherence",xlab="Frequency",ylab="spectrum")

# Phase
plot(fr,Arg(X.p.sp[,1,2]),type="l",ylim=c(-pi,pi),main="Phase",xlab="Frequency",ylab="spectrum")

##### NSA Analysis #####

# We have chosen the completions series to do the NSA analysis because it
# is better 'well-behaved' than the starts series.

### Autocorrelations and Unit Root Detection ###

acfcompNSATS=acf(compNSATS, lag.max=100, main="ACF Completions")
pacfcompNSATS=pacf(compNSATS, lag.max=100, main="PACF Completions")

# Taking first difference of the series

diffCompNSATS=diff(compNSATS) 

# Performing Augmented Dickey-Fuller unit root test on original series

test1=ur.df(compNSATS, type=c("none"), selectlags = c("AIC"))
summary(test1)
test2=ur.df(compNSATS, type=c("none"), selectlags = c("BIC")) 
summary(test2)
test3=ur.df(compNSATS, type=c("drift"), selectlags = c("AIC"))
summary(test3)
test4=ur.df(compNSATS, type=c("drift"), selectlags = c("BIC"))
summary(test4)

# Performing Augmented Dickey-Fuller unit root test on diff series

test5=ur.df(diffcompNSATS, type=c("none"), selectlags = c("AIC")) 
summary(test5)
test6=ur.df(diffcompNSATS, type=c("none"), selectlags = c("BIC"))
summary(test6)
test7=ur.df(diffcompNSATS, type=c("drift"), selectlags = c("AIC"))
summary(test7)
test8=ur.df(diffcompNSATS, type=c("drift"), selectlags = c("BIC"))
summary(test8)

# As expected, we still observe that the series is non-stationary. We take one difference
# and reapply the tests observing a strong rejection indicating that the non-stationary
# component was successfully ruled out with this filter.

### Spectral density Estimation ###

# Raw periodogram
spectrum(diffCompNSATS, log="no")

# Autoregressive spectral density estimation:
spec.ar(diffCompNSATS, log="no")

# Spectral-window (Smoothed Periodogram) Spectral Estimation

spectralWindow(diffCompNSATS,seq(3,10,1),1)
spectralWindow(diffCompNSATS,seq(3,10,1),2)
spectralWindow(diffCompNSATS,seq(3,10,1),3)

# By inspection, we decide on using c(1,1) for NSA completions.
k = kernel("daniell",c(1,1))
par(mfrow=c(1,1))
spec.pgram(diffCompNSATS, k, taper=0, log="no", plot=TRUE)

# Lag-window spectral estimation:

lagWindow(diffCompNSATS, seq(20,110,5), "", layout="group") #should change in the function
lagWindow(diffCompNSATS, seq(70,110,5), "", layout="group") #the kernel for c(1,1)

# By inspection, we see that from 80 lags and on the difference between the estimates
# is not visually significant so we chosse 90 for the sake of parcimony
lagWindow(diffCompNSATS, 90, "Series: diffcompNSATS\nLag Window", layout="alone")