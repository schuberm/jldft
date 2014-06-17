## Code to solve Poisson's equation

## Compute distances dr to center point in cell
dr= sqrt(sum((ones(prod(S),1)*sum(R,2)'/2-r).^2,2)); ## <=== CODE INSERTION # 1

## Compute two normalized Gaussians (widths 0.50 and 0.75)
sigma1=0.75
g1=exp(-dr.^2/(2*sigma1^2))/sqrt(2*pi*sigma1^2)^3

sigma2=0.50
g2=exp(-dr.^2/(2*sigma2^2))/sqrt(2*pi*sigma2^2)^3


## Define charge density as the difference
n=g2-g1

## Check norms and integral (should be near 1 and 0, respectively)
@sprintf("Normalization check on g1: %20.16f\n",sum(g1)*det(R)/prod(S))
@sprintf("Normalization check on g2: %20.16f\n",sum(g2)*det(R)/prod(S))
@sprintf("Total charge check: %20.16f\n",sum(n)*det(R)/prod(S))

## Visualize slices through center of cell
for dir=1:3
  text=sprintf("n#d=#d slice of n",dir,S(dir)/2)
  title(text)
  fprintf("#s (Hit <enter>... )\n",text)
  mesh(slice(n,S,S(dir)/2,dir)); pause
end

#### Solve Poisson's equation
##
phi=cI(Linv(-4*pi*O(cJ(n))))
####Due to rounding; tiny imaginary parts creep into the solution.  Eliminate
####by taking the real part.
phi=real(phi)
##
#### Visualize slices through center of cell
for dir=1:3
  text=sprintf("n#d=#d slice of phi",dir,S(dir)/2)
  title(text)
  fprintf("#s (Hit <enter>... )\n",text)
  mesh(slice(phi,S,S(dir)/2,dir)); pause
end
##
#### Check total Coulomb energy
Unum=0.5*real(cJ(phi)'*O(cJ(n)))
Uanal=((1/sigma1+1/sigma2)/2-sqrt(2)/sqrt(sigma1^2+sigma2^2))/sqrt(pi)
@sprintf("Numeric, analytic Coulomb energy: %20.16f,%20.16f\n",Unum,Uanal)
