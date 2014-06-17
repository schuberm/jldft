# VWN parameterization of the exchange correlation energy
function out=excVWN(n)
  ## Constants
  X1 = 0.75*(3.0/(2.0*pi))^(2.0/3.0)
  A  =  0.0310907
  x0 = -0.10498
  b  = 3.72744
  c  = 12.9352
  Q  = sqrt(4*c-b*b)
  X0 = x0*x0+b*x0+c

  rs=(4*pi/3*n).^(-1/3); ## Added internal conversion to rs
  
  x=sqrt(rs); X=x.*x+b*x+c

  out=-X1./rs ...
    + A*( ...
	+log(x.*x./X)+2*b/Q*atan(Q./(2*x+b)) ...
	-(b*x0)/X0*( ...
		    log((x-x0).*(x-x0)./X)+2*(2*x0+b)/Q*atan(Q./(2*x+b)) ...
		    ) ...
	)
endfunction
