#! /usr/bin/env python

import setpath
from PyTrilinos import Epetra
from Numeric import *

length  = 12

print "\n--------------"
shape = (2,2,3)
print "Create a", shape, "array:"
pyArray = arange(length, typecode=Float64)
pyArray.shape = shape
print "pyArray length         =", length
print "pyArray.shape          =", pyArray.shape
print "pyArray.typecode()     =", pyArray.typecode()
print "pyArray.iscontiguous() =", pyArray.iscontiguous()
print "pyArray =", pyArray
print "pyArray[0,0,0] = ", pyArray[0,0,0]
print "pyArray[0,0,1] = ", pyArray[0,0,1]
print "pyArray[0,0,2] = ", pyArray[0,0,2]
print "pyArray[0,1,0] = ", pyArray[0,1,0]
print "pyArray[0,1,1] = ", pyArray[0,1,1]
print "pyArray[0,1,2] = ", pyArray[0,1,2]
print "pyArray[1,0,0] = ", pyArray[1,0,0]
print "pyArray[1,0,1] = ", pyArray[1,0,1]
print "pyArray[1,0,2] = ", pyArray[1,0,2]
print "pyArray[1,1,0] = ", pyArray[1,1,0]
print "pyArray[1,1,1] = ", pyArray[1,1,1]
print "pyArray[1,1,2] = ", pyArray[1,1,2]

print "\n--------------"
print "Reorder axes of pyArray to create pyArray_2"
print "and then use it to create numPyarray_2."
print "This creates a non-contiguous array WITHOUT "
print "stride-one acces on the last dimenssion:"
pyArray_2 = swapaxes(pyArray,0,2)
print "pyArray_2 =", pyArray_2

print "\n--------------"
print "Reorder axes of pyArray to create pyArray_3"
print "and then use it to create numPyarray_3:"
print "This creates a non-contiguous array WITH "
print "stride-one acces on the last dimenssion:"
pyArray_3 = swapaxes(pyArray,0,1)
print "pyArray_3 =", pyArray_3

print "\n-- Epetra Load Test ----------------------"
print "\n----------------------------------------"
print   "-- Load Tests --------------------------"
print   "----------------------------------------"

length  = 12
comm = Epetra.SerialComm()
pbm = Epetra.BlockMap(length,1,0,comm)

# 1-D Contiguous:
print "\n-- 1-D Contiguous ----------------------"
pyArray1DCon         = arange(length, typecode=Float64)
pyArray1DCon.shape   = (length,)
epetraVector1DCon    = Epetra.Vector(pbm, pyArray1DCon)
print "epetraVector1DCon   = \n", epetraVector1DCon

# 2-D Contiguous:
print "\n-- 2-D Contiguous ----------------------"
pyArray2DCon         = arange(length, typecode=Float64)
pyArray2DCon.shape   = (2,6)
epetraVector2DCon    = Epetra.Vector(pbm, pyArray2DCon)
print "epetraVector2DCon   = \n", epetraVector2DCon

# 2-D Non-Contiguous:
print "\n-- 2-D Non-Contiguous ------------------"
pyArray2DNonCon      = swapaxes(pyArray2DCon,0,1)
epetraVector2DNonCon = Epetra.Vector(pbm, pyArray2DNonCon)
print "epetraVector2DNonCon= \n", epetraVector2DNonCon

# 3-D Contiguous:
print "\n-- 3-D Contiguous ----------------------"
pyArray3DCon         = arange(length, typecode=Float64)
pyArray3DCon.shape   = (2,2,3)
epetraVector3DCon    = Epetra.Vector(pbm, pyArray3DCon)
print "epetraVector3DCon   = \n", epetraVector3DCon

# 3-D Non-Contiguous:
print "\n-- 3-D Non-Contiguous ------------------"
pyArray3DNonCon      = swapaxes(pyArray3DCon,0,2)
epetraVector3DNonCon = Epetra.Vector(pbm, pyArray3DNonCon)
print "epetraVector3DNonCon= \n", epetraVector3DNonCon

# 3-D Non-Contiguous Stride 1:
print "\n-- 3-D Non-Contiguous Stride 1 ---------"
pyArray3DNonCon     = swapaxes(pyArray3DCon,0,1)
epetraVector3DNonCon  = Epetra.Vector(pbm, pyArray3DNonCon)
print "epetraVector3DNonCon= \n", epetraVector3DNonCon

length              = 24
pbm = Epetra.BlockMap(length,1,0,comm)

# 4-D Contiguous:
print "\n-- 4-D Contiguous ----------------------"
pyArray4DCon        = arange(length, typecode=Float64)
pyArray4DCon.shape  = (2,2,2,3)
epetraVector4DCon    = Epetra.Vector(pbm, pyArray4DCon)
print "epetraVector4DCon   = \n", epetraVector4DCon

# 4-D Non-Contiguous:
print "\n-- 4-D Non-Contiguous ------------------"
pyArray4DNonCon     = swapaxes(pyArray4DCon,0,3)
epetraVector4DNonCon = Epetra.Vector(pbm, pyArray4DNonCon)
print "epetraVector4DNonCon= \n", epetraVector4DNonCon

# 4-D Non-Contiguous Stride 1:
print "\n-- 4-D Non-Contiguous Stride 1 ---------"
pyArray4DNonCon     = swapaxes(pyArray4DCon,0,1)
epetraVector4DNonCon = Epetra.Vector(pbm, pyArray4DNonCon)
print "epetraVector4DNonCon= \n", epetraVector4DNonCon

print "\n----------------------------------------"
print   "-- Unload Tests ------------------------"
print   "----------------------------------------"

length  = 12
pbm = Epetra.BlockMap(length,1,0,comm)

# 1-D Contiguous:
print "\n-- 1-D Contiguous ----------------------"
pyArray1DCon         = arange(length, typecode=Float64)
pyArray1DCon.shape   = (length,)
epetraVector1DCon    = Epetra.Vector(pbm, pyArray1DCon)
epetraVector1DCon.Scale(3.0)
print "epetraVector1DCon   = \n", epetraVector1DCon
pyArray1DCon[:] = epetraVector1DCon

# 2-D Contiguous:
print "\n-- 2-D Contiguous ----------------------"
pyArray2DCon         = arange(length, typecode=Float64)
pyArray2DCon.shape   = (2,6)
epetraVector2DCon    = Epetra.Vector(pbm, pyArray2DCon)
epetraVector2DCon.Scale(5.0)
print "epetraVector2DCon   = \n", epetraVector2DCon
pyArray2DCon[:] = epetraVector2DCon

# 2-D Non-Contiguous:
print "\n-- 2-D Non-Contiguous ------------------"
pyArray2DCon         = arange(length, typecode=Float64)
pyArray2DCon.shape   = (2,6)
pyArray2DNonCon      = swapaxes(pyArray2DCon,0,1)
epetraVector2DNonCon = Epetra.Vector(pbm, pyArray2DNonCon)
epetraVector2DNonCon.Scale(7.0)
print "epetraVector2DNonCon= \n", epetraVector2DNonCon
pyArray2DNonCon[:] = epetraVector2DNonCon

# 3-D Contiguous:
print "\n-- 3-D Contiguous ----------------------"
pyArray3DCon         = arange(length, typecode=Float64)
pyArray3DCon.shape   = (2,2,3)
epetraVector3DCon    = Epetra.Vector(pbm, pyArray3DCon)
epetraVector3DCon.Scale(11.0)
print "epetraVector3DCon   = \n", epetraVector3DCon
pyArray3DCon[:] = epetraVector3DCon

# 3-D Non-Contiguous:
print "\n-- 3-D Non-Contiguous ------------------"
pyArray3DCon         = arange(length, typecode=Float64)
pyArray3DCon.shape   = (2,2,3)
pyArray3DNonCon      = swapaxes(pyArray3DCon,0,2)
epetraVector3DNonCon = Epetra.Vector(pbm, pyArray3DNonCon)
epetraVector3DNonCon.Scale(13.0)
print "epetraVector3DNonCon= \n", epetraVector3DNonCon
pyArray3DNonCon[:] = epetraVector3DNonCon

# 3-D Non-Contiguous Stride 1:
print "\n-- 3-D Non-Contiguous Stride 1 ---------"
pyArray3DCon         = arange(length, typecode=Float64)
pyArray3DCon.shape   = (2,2,3)
pyArray3DNonCon     = swapaxes(pyArray3DCon,0,1)
epetraVector3DNonCon  = Epetra.Vector(pbm, pyArray3DNonCon)
epetraVector3DNonCon.Scale(17.0)
print "epetraVector3DNonCon= \n", epetraVector3DNonCon
pyArray3DNonCon[:] = epetraVector3DNonCon

length              = 24
pbm = Epetra.BlockMap(length,1,0,comm)

# 4-D Contiguous:
print "\n-- 4-D Contiguous ----------------------"
pyArray4DCon        = arange(length, typecode=Float64)
pyArray4DCon.shape  = (2,2,2,3)
epetraVector4DCon    = Epetra.Vector(pbm, pyArray4DCon)
epetraVector4DCon.Scale(19.0)
print "epetraVector4DCon   = \n", epetraVector4DCon
pyArray4DCon[:] = epetraVector4DCon

# 4-D Non-Contiguous:
print "\n-- 4-D Non-Contiguous ------------------"
pyArray4DCon        = arange(length, typecode=Float64)
pyArray4DCon.shape  = (2,2,2,3)
pyArray4DNonCon     = swapaxes(pyArray4DCon,0,3)
epetraVector4DNonCon = Epetra.Vector(pbm, pyArray4DNonCon)
epetraVector4DNonCon.Scale(23.0)
print "epetraVector4DNonCon= \n", epetraVector4DNonCon
pyArray4DNonCon[:] = epetraVector4DNonCon

# 4-D Non-Contiguous Stride 1:
print "\n-- 4-D Non-Contiguous Stride 1 ---------"
pyArray4DCon        = arange(length, typecode=Float64)
pyArray4DCon.shape  = (2,2,2,3)
pyArray4DNonCon     = swapaxes(pyArray4DCon,0,1)
epetraVector4DNonCon = Epetra.Vector(pbm, pyArray4DNonCon)
epetraVector4DNonCon.Scale(27.0)
print "epetraVector4DNonCon= \n", epetraVector4DNonCon
pyArray4DNonCon[:] = epetraVector4DNonCon
