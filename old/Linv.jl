# Linv operator (acting on 3d data sets)
#
# Usage: out=Linv(in)
#
# in: input 3d data set
# out: output 3d data set
#

function Linv(in)
  global gbl_R; ## Must declare all globals with such statements to access them
  global gbl_G2
  ## Operator definition (multiplication by volume)
  #out= (-det(gbl_R)*diagm(sparse(complex(gbl_G2))))\in;
  out= \((-det(gbl_R)*diagm(sparse(complex(gbl_G2[:])))),in); ## <=== YOUR CODE HERE
  return out
end