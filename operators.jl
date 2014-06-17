#Define globals
global gbl_S;
global gbl_Vdual;
global gbl_R;
global gbl_G2;

# Forward transform (acting on 3d data sets)
#
# Usage: out=cI(in)
#
# in: input 3d data set
# out: output 3d data set
#

function cI(in)
  #global gbl_S; ## Must declare all globals with such statements to access them
  ## Operator definition (multiplication by volume)
  ##out=zeros(size(in))
  out=Array(Complex{Float64},size(in,1),size(in,2))
  for col=1:size(in,2)
    out[:,col]= fft3(in[:,col],gbl_S,1); ## <=== YOUR CODE HERE
  end
  return out
end

# Dual Forward transform (acting on 3d data sets)
#
# Usage: out=cIdag(in)
#
# in: input 3d data set
# out: output 3d data set
#

function cIdag(in)
  #global gbl_S; ## Must declare all globals with such statements to access them
  #out=zeros(size(in))
  out=Array(Complex{Float64},size(in,1),size(in,2))
  for col=1:size(in,2)
    out[:,col]= fft3(in[:,col],gbl_S,-1); ## <=== YOUR CODE HERE
  end
  return out
end

# Inverse transform (acting on 3d data sets)
#
# Usage: out=cJ(in)
#
# in: input 3d data set
# out: output 3d data set
#

function cJ(in)
  #global gbl_S; ## Must declare all globals with such statements to access them
  out=Array(Complex{Float64},size(in,1),size(in,2))
  #out=zeros(size(in,1),size(in,2))
  for col=1:size(in,2)
    out[:,col]= fft3(in[:,col],gbl_S,-1)/(gbl_S[1]*gbl_S[2]*gbl_S[3]); ## <=== YOUR CODE HERE
  end
  return out
end

# Dual inverse transform (acting on 3d data sets)
#
# Usage: out=cJdag(in)
#
# in: input 3d data set
# out: output 3d data set
#

function cJdag(in)
  #global gbl_S; ## Must declare all globals with such statements to access them
  #out=zeros(size(in))
  out=Array(Complex{Float64},size(in,1),size(in,2))
  for col=1:size(in,2)
    out[:,col]= fft3(in[:,col],gbl_S,1)/prod(gbl_S); ## <=== YOUR CODE HERE
  end
  return out
end

# Diagonal elements of the outer product of two matrices
#
# Usage: out=diagouter(A,B)
#

function diagouter(A,B)
  out=sum(A.*conj(B),2)
  return out
end

# out=diag(a)*B
#

function Diagprod(a,B)
  out=(a*ones(1,size(B,2))).*B
  return out
end

## out=fft3(dat,N,s) - computes 3d fft (dimensions in N) of sign s
## out[l,m,n)=sum_{a,b,c} exp(2 pi s i *(a*l/N(1)+b*m/N(2)+c*n/N(3))*in(a,b,c].
## Notes: 1) fortran/matlab ordering assumed, ordering in mem is out[1,1,1]
##           out[2,1,1), ..., out(N(1),1,1), out(1,2,1], ...
##        2) a=fft3(b,N,1) => b=fft3(a,N,-1)/prod(N); ie., fft(dat,N,1) and
##           fft(dat,N,-1) are inverses except for normalization by
##           N(1)*N(2)*N(3)
function fft3(dat,N,s)
#  tic
  if s==1
    out=reshape(ifftn(reshape(dat,N[1],N[2],N[3]))*prod(N),size(dat))
  else
    out=reshape(fftn(reshape(dat,N[1],N[2],N[3])),size(dat))
  end

#  Ntot=prod(size(dat))
#  MFLOPS=5*Ntot*log2(Ntot)/1e6/toc
  return out
end

# 
#
# Usage: out=getE(in)
#
# in: W; expansion coefficients for Ns unconstrained wavefunctions
# out: output 3d data set
#

function getE(in)
  #global gbl_Vdual; ## Must declare all globals with such statements to access them
  Y=in*inv(sqrtm(in'*O(in)))
  n=real(diagouter(cI(Y),cI(Y)))
  out=real(-0.5*trace(Y'*L(Y))+gbl_Vdual'*n)
  return out
end

# Gradient of E
#
# Usage: out=getgrad(in)
#
# in: W; expansion coefficients for Ns unconstrained wavefunctions
# out: output 3d data set
#

function getgrad(in)
  #global gbl_Vdual; ## Must declare all globals with such statements to access them
  U=in'*O(in); 
  Uinv=inv(U); 
  #out=(H(in)-O(in*inv(in"*O(in)))*(in'*H(in)))*inv(in"*O(in))
  out=(H(in)-O(in)*Uinv*(in'*H(in)))*Uinv
  return out
end

function getn(Y,f)
  n=real(diagouter(cI(Y)*f,cI(Y))); ## Charge density
#  n=0
#  for col=1:size(Y,2)
#    tmp=cI(Y(:,col))
##    n=n+f*tmp.*conj(tmp)
#    n=n+f(col)*tmp.*conj(tmp)
#  end
  return n
end

function getPsi(W)
  ## Orthonormalize
  Y=W*inv(sqrtm(W'*O(W)))

  ## Wompute subspace Hamiltonian and diagonalize
  mu=Y'*H(Y)

  (U, epsilon)=eig(mu); epsilon=real(diag(epsilon))
  Psi=Y*U
  (Psi, epsilon)
end

# 
#
# Usage: out=H(in)
#
# in: W; expansion coefficients for Ns unconstrained wavefunctions
# out: output 3d data set
#

function H(in)
  #global gbl_Vdual; ## Must declare all globals with such statements to access them
  #Vtild=cJdag(O(cJ(gbl_Vdual)))
  out=-0.5*L(in)+cIdag(Diagprod(gbl_Vdual,cI(in)))
  #out=-0.5*L(in);#+cIdag(Diagprod(Vtild,cI(in)))
  return out
end

# Laplacian operator (acting on 3d data sets)
#
# Usage: out=L(in)
#
# in: input 3d data set
# out: output 3d data set
#

function L(in)
  #global gbl_R; ## Must declare all globals with such statements to access them
  #global gbl_G2;
  ## Operator definition (multiplication by volume)
  #out= -det(gbl_R)*diag(gbl_G2)*in; ## <=== YOUR CODE HERE
  out= -det(gbl_R)*gbl_G2*ones(1,size(in,2)).*in; ## <=== YOUR CODE HERE
  return out
end

# Linv operator (acting on 3d data sets)
#
# Usage: out=Linv(in)
#
# in: input 3d data set
# out: output 3d data set
#

function Linv(in)
  #global gbl_R; ## Must declare all globals with such statements to access them
  #global gbl_G2
  ## Operator definition (multiplication by volume)
  #out= (-det(gbl_R)*diagm(sparse(complex(gbl_G2))))\in;
  out= \((-det(gbl_R)*diagm(sparse(complex(gbl_G2[:])))),in); ## <=== YOUR CODE HERE
  return out
end

# Overlap operator (acting on 3d data sets)
#
# Usage: out=O(in)
#
# in: input 3d data set
# out: output 3d data set
#
# Uses GLOBAL variable(s) ---
# gbl_R: Lattice vectors

function O(in)
  #global gbl_R; ## Must declare all globals with such statements to access them
  out= det(gbl_R)*in
  return out
end

function Q(in,U)
  (V,mu)=eig(U); mu=diag(mu)

  denom=sqrt(mu)*ones(1,length(mu)); denom=denom+denom'
  out=V*( (V*in*V)./denom )*V
  return out
end

# Steepest descent
#
# Usage: out=sd(in,Nit)
#
# in: W; expansion coefficients for Ns unconstrained wavefunctions
# out: output 3d data set

function sd(in,Nit)
  alpha=3E-5
  out=in
  for iter=1:1:Nit
    out=out-alpha*getgrad(out)
    @sprintf("Energy: %20.16f\n",getE(out))
    #old=out
  end
  return out
end
