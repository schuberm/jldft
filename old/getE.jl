# 
#
# Usage: out=getE(in)
#
# in: W; expansion coefficients for Ns unconstrained wavefunctions
# out: output 3d data set
#

function out=getE(in)
  global gbl_Vdual; ## Must declare all globals with such statements to access them
  Y=in*inv(sqrtm(in'*O(in)))
  n=real(diagouter(cI(Y),cI(Y)))
  out=real(-0.5*trace(Y"*L(Y))+gbl_Vdual"*n)
endfunction
