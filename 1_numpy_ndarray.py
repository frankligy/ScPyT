import numpy as np

# understanding numpy array
a = np.array([[1,2,3],[4,5,6]])
a.strides
a.dtype


# slicing and indexing
b = np.array([[1,2,3,4,5,6,7,8,9,10],
             [4,5,6,7,8,9,20,11,12,13],
             [1,2,3,4,5,6,7,8,9,9]])

# what is a slice object
test = slice(1,10,2)

# basic indexing only return a "view", it will change the original array content
b0 = b[1:3,4:7]
b0[0,1] = 99
b1 = b[:,(4,7)]
b1[0,1] = 99


# understand dtype
d_type = np.dtype('<U8')
d_type.byteorder
d_type.itemsize
d_type.name

# structural array and record array
sa = np.array([('John',[88,95,100]),('Mary',[77,88,68])],
              dtype=[('student','<U8'),('grades','<i4',(3,))])
sa['student']
ra = sa.view(np.recarray)
ra.student