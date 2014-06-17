function getPsi(W)
  ## Orthonormalize
  Y=W*inv(sqrtm(W'*O(W)))

  ## Wompute subspace Hamiltonian and diagonalize
  mu=Y'*H(Y)

  [U, epsilon]=eig(mu); epsilon=real(diag(epsilon))
  Psi=Y*U
endfunction
  [Psi, epsilon]
end
