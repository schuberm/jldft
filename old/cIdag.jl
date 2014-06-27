# Dual Forward transform (acting on 3d data sets)
#
# Usage: out=cIdag(in)
#
# in: input 3d data set
# out: output 3d data set
#

function cIdag(in)
  global gbl_S; ## Must declare all globals with such statements to access them
  out=zeros(size(in))
  for col=1:size(in,2)
    out[:,col]= fft3(in[:,col],gbl_S,-1); ## <=== YOUR CODE HERE
  end
  return out
end
