#! /usr/bin/env python

import setpath
from PyTrilinos import Epetra

# Problem size
n = 4

# SerialDenseVector, default constructor
v1 = Epetra.SerialDenseVector()
assert(v1.CV() == Epetra.Copy)
assert(v1.Length() == 0)
v1.Size(2*n)
assert(v1.Length() == 2*n)
v1.Resize(n)
assert(v1.Length() == n)
for i in range(n):
    v1[i] = i*i
print "v1 =", v1
for i in range(n):
    assert(v1[i] == i*i)
    assert(v1(i) == i*i)

# SerialDenseVector, sized constructor
v2 = Epetra.SerialDenseVector(n)
assert(v2.CV() == Epetra.Copy)
assert(v2.Length() == n)
for i in range(n):
    v2[i] = i*(n-i)
print "v2 =", v2
for i in range(n):
    assert(v2[i] == i*(n-i))
    assert(v2(i) == i*(n-i))

# SerialDenseVector, copy constructor
v3 = Epetra.SerialDenseVector(v2)
assert(v3.CV() == Epetra.Copy)
assert(v3.Length() == v2.Length())
print "v3 =", v3
for i in range(v3.Length()):
    assert(v3[i] == v2[i])
    assert(v3(i) == v2(i))

# SerialDenseMatrix, default constructor
m1 = Epetra.SerialDenseMatrix()
assert(m1.CV() == Epetra.Copy)
assert(m1.M()  == 0)
assert(m1.N()  == 0)
m1.Shape(2,1)
assert(m1.M()  == 2)
assert(m1.N()  == 1)
m1.Reshape(4,2)
assert(m1.M()  == 4)
assert(m1.N()  == 2)

# SerialDenseMatrix, sized constructor
m2 = Epetra.SerialDenseMatrix(n,n)
assert(m2.LDA() == m2.M())
for i in range(n):
    if (i>0): m2[i,i-1] = 1
    m2[i,i] = -2
    if (i<n-1): m2[i,i+1] = 1
print "m2    =", m2

# SerialDenseMatrix, copy constructor
scalar = 2
m3 = Epetra.SerialDenseMatrix(m2)
m3.Scale(scalar)
for i in range(m3.M()):
    for j in range(m3.N()):
        assert(m3[i,j] == scalar*m2[i,j])
assert(m3.NormOne() == scalar*m2.NormOne())
assert(m3.NormInf() == scalar*m2.NormInf())
# m4 = Epetra.SerialDenseMatrix(m2)
# m4 += m2    Why does this seg fault?
# assert(m4.NormOne() == m2.NormOne())
# assert(m4.NormInf() == m2.NormInf())

# SerialDenseSolver
m2inv = Epetra.SerialDenseMatrix(m2)
for i in range(n):
    if (i>0): assert(m2[i,i-1] == 1)
    assert(m2[i,i] == -2)
    if (i<n-1): assert(m2[i,i+1] == 1)
    for j in range(n):
        assert(m2inv[i,j] == m2[i,j])
m2sys = Epetra.SerialDenseSolver()
m2sys.SetMatrix(m2inv)
m2sys.Invert()
a = m2sys.FactoredMatrix()
for i in range(n):
    for j in range(n):
        assert(a[i,j] == m2inv[i,j])
print "m2inv =", m2inv

# Matrix multiply
identity = Epetra.SerialDenseMatrix(n,n)
identity.Multiply("N","N",1,m2,m2inv,0)
eps = 1.0e-15
for i in range(identity.M()):
    for j in range(identity.N()):
        if (i==j):
            assert(abs(identity[i,j]-1) < eps)
        else:
            assert(abs(identity[i,j]  ) < eps)

# Test exception handling for extended __setitem__
# and __getitem__ methods
try:
    m2[0,1,2] = 3.0
    raise RuntimeError, "Too many indexes not caught"
except IndexError:
    pass

try:
    m2[3.14]  = 3.14
    raise RuntimeError, "Float index not caught"
except IndexError:
    pass

try:
    v = m2[0,1,2]
    raise RuntimeError, "Too many indexes not caught"
except IndexError:
    pass

try:
    v = m2[3.14]
    raise RuntimeError, "Float index not caught"
except IndexError:
    pass
