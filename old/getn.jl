function n=getn(Y,f)
  n=real(diagouter(cI(Y)*f,cI(Y))); ## Charge density
#  n=0
#  for col=1:size(Y,2)
#    tmp=cI(Y(:,col))
##    n=n+f*tmp.*conj(tmp)
#    n=n+f(col)*tmp.*conj(tmp)
#  end
end
