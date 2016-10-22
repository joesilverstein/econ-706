spectralWindow = function(series,sequence,numPasses) 
  # sequence: vector of values of m to loop through
  # numPasses: number of times to run the filter
  {
  
for (i in sequence) {
  par(mfrow=c(2,2))
  spectrum(series, log="no")
  spec.ar(series, log="no")
  passesSeq = i*rep(1,numPasses)
  k = kernel("daniell",passesSeq)
  spec.pgram(series, k, taper=0, log="no")
  print(passesSeq)
  
  out = readline("Press <return> to continue or enter <q> to quit: ")
  if (out=="q") { break }
}

}