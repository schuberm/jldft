# Dual inverse transform (acting on 3d data sets)
#
# Usage: out=cJdag(in)
#
# in: input 3d data set
# out: output 3d data set
#

function cJdag(in)
  global gbl_S; ## Must declare all globals with such statements to access them
  out=zeros(size(in))
  for col=1:size(in,2)
    out[:,col]= fft3(in[:,col],gbl_S,1)/prod(gbl_S); ## <=== YOUR CODE HERE
  end
  return out
end
