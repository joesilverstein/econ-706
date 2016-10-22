lagWindow = function(series, sequence, title, layout){
  
for (M in sequence) {
  
  # M = 50; # start with 12 lags because number of months in a year
  gammaHat = acf(series, lag.max=M, type="covariance", plot=FALSE) # first element of acf vector is the variance! 
  T = length(series)
  lambda = rep(0,T) # since lambda is symmetric, only need to define from 0 to T
  # Note: first element of lambda is lambda(0)=1
  
  # construct triangular kernel:
  for (t in 1:T) {
    if (t<=M) {
      lambda[t] = 1 - abs(t)/M 
    }
    else {
      lambda[t] = 0
    }
  }
  
  # Create grid of frequencies omega. The indices of the grid are the the indices j
  # mentioned on slide 78 at the top of the slide.
  gridOmega = rep(0,T/2)
  for (j in 1:T/2) {
    gridOmega[j] = (6/pi)*(2*pi*j/T) #rescaling for coherence with the other plots
  }
  f = rep(0,T/2) # density function (paired with gridOmega)
  
  for (j in 1:T/2) { # for each value of omega in the grid
    sum = 0
    # Calculate the sum. Since lambda is 0 if tau>M, we only need to add up the sum for indices <=M here
    for (t in 1:M) {
      sum = sum + lambda[t]*gammaHat$acf[t+1]*cos(2*pi*j*t/T) # the i goes away because e^(i*pi)=-1
      #print(lambda[t])
      #print(gammaHat$acf[t])
    }
    sum = 2*sum # by symmetry
    sum = sum + var(series) # need to add the term for tau=0, which is the variance gammaHat(0)
    f[j] = (2*pi)^(-1)*sum
  }
  
  if (layout=="group") {
    par(mfrow=c(2,2))
    spectrum(series, log="no")
    spec.ar(series, log="no")
    k = kernel("daniell",c(7,7))
    print(M)
    spec.pgram(series, k, taper=0, log="no")
    plot(gridOmega, f, type="l", main=title, xlab="frequency", ylab="spectrum")
    
    out = readline("Press <return> to continue or enter <q> to quit: ")
    if (out=="q") { break }
  }
  else {
    par(mfrow=c(1,1))
    print(M)
    plot(gridOmega, f, type="l", main=title, xlab="frequency", ylab="spectrum")
  }
}
  
}