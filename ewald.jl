## Code to solve Poisson's equation

## Compute distances dr to center point in cell
dr= sqrt(sum((ones(prod(S),1)*sum(R,2)'/2-r).^2,2)); ## <=== CODE INSERTION # 1

## Compute one normalized Gaussian (widths 0.25)
sigma1=0.25
g1=exp(-dr.^2/(2*sigma1^2))/(sqrt(2*pi*sigma1^2)^3*Z)

## Use structure factor to create all atoms
n=cI(cJ(g1).*Sf); 
#n=cJ(g1)#.*Sf); 
n=real(n)

## Check norms and integral (should be near 1 and 0, respectively)
@sprintf("Normalization check on g1: %20.16f\n",sum(g1)*det(R)/prod(S)*Z)
@sprintf("Total charge check: %20.16f\n",sum(n)*det(R)/prod(S))

## Visualize slices through center of cell
#for dir=1:3
#  text=sprintf("n#d=#d slice of n",dir,S(dir)/2)
#  title(text)
#  fprintf("#s (Hit <enter>... )\n",text)
#  mesh(slice(n,S,S(dir)/2,dir)); pause
#end

#### Solve Poisson's equation
##
phi=cI(Linv(-4*pi*O(cJ(n))))
####Due to rounding; tiny imaginary parts creep into the solution.  Eliminate
####by taking the real part.
phi=real(phi)
##
#### Visualize slices through center of cell
#for dir=1:3
#  text=sprintf("n#d=#d slice of phi",dir,S(dir)/2)
#  title(text)
# fprintf("#s (Hit <enter>... )\n",text)
#  mesh(slice(phi,S,S(dir)/2,dir)); pause
#end
##
#### Check total Coulomb energy
Unum=0.5*real(cJ(phi)'*O(cJ(n)))
Uself=Z^2/(2*sqrt(pi))*(1/sigma1)*size(X,1)
#@sprintf("Ewald energy: %20.16f\n",Unum-Uself)
