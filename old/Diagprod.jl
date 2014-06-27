# out=diag(a)*B
#

function Diagprod(a,B)
  out=(a*ones(1,size(B,2))).*B
  return out
end
