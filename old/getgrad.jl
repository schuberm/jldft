# Gradient of E
#
# Usage: out=getgrad(in)
#
# in: W; expansion coefficients for Ns unconstrained wavefunctions
# out: output 3d data set
#

function out=getgrad(in)
  global gbl_Vdual; ## Must declare all globals with such statements to access them
  U=in'*O(in); 
  Uinv=inv(U); 
  #out=(H(in)-O(in*inv(in"*O(in)))*(in'*H(in)))*inv(in"*O(in))
  out=(H(in)-O(in)*Uinv*(in'*H(in)))*Uinv
endfunction
