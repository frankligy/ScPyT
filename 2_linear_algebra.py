import numpy as np
import scipy.linalg as sl

# vector space part
a = np.random.randn(2,3)
a.T

b = np.random.randn(3,2)
a @ b

# matrix_rank function
x = np.random.rand(2,3)
np.linalg.matrix_rank(x)

# determinant
x = np.random.rand(4,4)
np.linalg.det(x)

# inner product, outer product, hadamard product of two 1D array
a = np.array([4,5,6])
b = np.array([7,8,9])
np.inner(a,b)
np.outer(a,b)
a * b

# how to project a to b
a = np.array([4,5,6])
b = np.array([7,8,9])
proj_b_a = np.inner(a,b) / np.inner(b,b) * b

# LU decomposition
a = np.random.randn(3,4)
p,l,u = sl.lu(a)

# QR decomposition
a = np.random.randn(3,4)
q,r = np.linalg.qr(a)

# eigendecompostion
ei = np.random.randn(4,4)
w,v = np.linalg.eig(ei)

# characteristic polynomial
a = np.random.randn(5,5)
np.linalg.det(np.trace(a)*np.identity(5)-a)

# SVD
svd = np.random.rand(4,5)
u,s,vh = np.linalg.svd(svd)

# norm
a = np.array([4,5,6])
np.linalg.norm(a,ord=3)

x = np.random.rand(2,3)
np.linalg.norm(x,ord='fro')

# einsum
x = np.random.rand(2,3)
np.einsum('ij -> ji',x)   # transpose
np.einsum('ij ->',x)    # sum
np.einsum('ij -> i',x)   # column sum
np.einsum('ij -> j',x)  # row sum

x = np.random.rand(2,3)
y = np.random.rand(5,3)
np.einsum('ij,kj -> ik',x,y)  # matrix multiplication

a = np.array([4,5,6])
b = np.array([7,8,9])
np.einsum('i,i ->',a,b)   # inner product
np.einsum('i,j ->ij',a,b)  # outer product
np.einsum('i,i ->i',a,b)   # hadamard product

y = np.random.rand(5,3)
np.einsum('ij -> j',y)   # diagonal
np.einsum('ij ->',y)   # trace

svd = np.random.rand(4,5)
u,s,vh = np.linalg.svd(svd)

ei = np.random.randn(4,4)
w,v = np.linalg.eig(ei)






