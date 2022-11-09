def GramSchmidt(A):
  m = A.shape[0]
  n = A.shape[1]
  
  Q = np.zeros((m,n))
  R = np.zeros((n,n))
  
  for j in range(n):
    v = A[:,j]
    if j > 0:
      for i in range(j):
        R[i,j] = np.dot(Q[:,i], A[:,j])
        v = v - R[i,j] * Q[:,i]
    R[j,j] = np.linalg.norm(v, 2)
    Q[:,j] = v / R[j,j]
  
  return Q, R
