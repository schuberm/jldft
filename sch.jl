## Code to solve Shrodinger's equation
global gbl_f
gbl_f=1
## Compute square of distances dr to center point in cell
dr2= sum((ones(prod(S),1)*sum(R,2)'/2-r).^2,2)

dr=r-ones(prod(S),1)*diag(R)'/2
dr=sqrt(sum(dr.^2,2))

## Compute harmonic oscillator potential of frequency omega
omega=4
V=0.5*omega*dr2
Vc=0.5*4*dr.^2

if(V!=Vc)
	@sprintf("V is different" )
end
##compute dual representation of potential
#cJdag(O(cJ(V)))
gbl_Vdual=cJdag(O(cJ(V)))

## Finite difference test
Ns=4; ## Number of states
#randn("seed",0.2004)
#randn(1)
W=(randn(prod(S),Ns)+im*randn(prod(S),Ns))
#more off; ## View output as it is computed
#fdtest(W)

W=W*sqrtm(inv(W'*O(W)))

sd(W,250)

# ## Allow for more digits in printouts
# format long
# ## Converge using steepest descents
# W=sd(W,400)
# ## Extract and display final results
# [Psi, epsilon]=getPsi(W)
# for st=1:Ns
# printf("=== State # #d, Energy = #f ===\",st,epsilon(st))
# dat=abs(cI(Psi(:,st))).^2
# for k=1:3
# sl=slice(dat,S,S(k)/2,k)
# name=sprintf("psi#d_m#d.ppm",st,k)
# ppm(name,sl*0.3,sl,sl); system(["display " name "&"])
# end
# end
