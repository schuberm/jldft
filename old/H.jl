# 
#
# Usage: out=H(in)
#
# in: W; expansion coefficients for Ns unconstrained wavefunctions
# out: output 3d data set
#

function out=H(in)
  global gbl_Vdual; ## Must declare all globals with such statements to access them
  #Vtild=cJdag(O(cJ(gbl_Vdual)))
  out=-0.5*L(in)+cIdag(Diagprod(gbl_Vdual,cI(in)))
  #out=-0.5*L(in);#+cIdag(Diagprod(Vtild,cI(in)))
endfunction
