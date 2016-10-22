# Nonparametric density estimation using standard normal kernel
# There are problems with this, so don't use it for now.
gauss = function(x) 1/sqrt(2*pi) * exp(-(x^2)/2) # standard normal kernel
diffStartSpecGrid = seq(from = min(diffStartSpec$freq), to = max(diffStartSpec$freq), by = 0.01)
difstart = as.numeric(diffStartSpec)
n = length(diffStartSpec$freq)
h = 0.1 # bandwidth selection
bumps = sapply(diffStartSpec$freq, function(a) gauss((diffStartSpecGrid - a)/h)/(n * h))
plot(diffStartSpecGrid, rowSums(bumps), ylab = "Density",type = "l", xlab = "Frequency", lwd = 2)
# It looks like the standard normal kernel smooths too much, and in the wrong way
