function out=Q(in,U)
  [V,mu]=eig(U); mu=diag(mu)

  denom=sqrt(mu)*ones(1,length(mu)); denom=denom+denom'
  out=V*( (V"*in*V)./denom )*V"
end
