from __future__ import print_function
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import lsqr

row = [0, 0, 1, 1, 2, 2]
col = [0, 1, 0, 1, 0, 1]
data = [1, 10, 3, 4, 2, 6]

A = coo_matrix((data, (row, col)), shape=(3, 2)).tocsc()
b = [1, 2, 3]
print(lsqr(A, b, show=True)[0])
