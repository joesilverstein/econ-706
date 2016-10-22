var.spec = function(fr,X.p)
  # fr . vector of frequencies
  # X.p . AR(p) model estimated by R base function ar
{ 
  
p = X.p$order; nf=length(fr)
im=complex(real=0,imaginary=1); # i
sigma = X.p$var.pred
k = length(sigma[1,]) 
Id = diag(1,nrow=k,ncol=k)
sp = array(dim=c(nf,k,k))
for (w in 1:nf) # for each of the frequencies in the grid
{ A = Id
  for (l in 1:p) { A = A-X.p$ar[l,,]*exp(-im*(1/2)*fr[w]*l) }
  A = solve(A)
  sp[w,,] = A%*%sigma%*%t(Conj(A))/(2*pi)
}
return(sp)

} 