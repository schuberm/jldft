# Laplacian operator (acting on 3d data sets)
#
# Usage: out=L(in)
#
# in: input 3d data set
# out: output 3d data set
#

function L(in)
  global gbl_R; ## Must declare all globals with such statements to access them
  global gbl_G2;
  ## Operator definition (multiplication by volume)
  #out= -det(gbl_R)*diag(gbl_G2)*in; ## <=== YOUR CODE HERE
  out= -det(gbl_R)*gbl_G2*ones(1,size(in,2)).*in; ## <=== YOUR CODE HERE
  return out
end