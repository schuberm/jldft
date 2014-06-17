# Diagonal elements of the outer product of two matrices
#
# Usage: out=diagouter(A,B)
#

function diagouter(A,B)
  out=sum(A.*conj(B),2)
  return out
end
