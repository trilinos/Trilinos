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
print "Create a NumPyArray from above pyArray:"
numPyArray = Epetra.NumPyArray(pyArray)
print "numPyArray =\n", numPyArray

print "\n--------------"
print "Reorder axes of pyArray to create pyArray_2"
print "and then use it to create numPyarray_2."
print "This creates a non-contiguous array WITHOUT "
print "stride-one acces on the last dimenssion:"
pyArray_2 = swapaxes(pyArray,0,2)
numPyArray_2 = Epetra.NumPyArray(pyArray_2)
print "pyArray_2 =", pyArray_2
print "numPyArray_2 =\n", numPyArray_2

print "\n--------------"
print "Reorder axes of pyArray to create pyArray_3"
print "and then use it to create numPyarray_3:"
print "This creates a non-contiguous array WITH "
print "stride-one acces on the last dimenssion:"
pyArray_3 = swapaxes(pyArray,0,1)
numPyArray_3 = Epetra.NumPyArray(pyArray_3)
print "pyArray_3 =", pyArray_3
print "numPyArray_3 =\n", numPyArray_3

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
numPyArray1DCon      = Epetra.NumPyArray(pyArray1DCon)
print "numPyArray1DCon     = \n", numPyArray1DCon
epetraVector1DCon    = Epetra.Vector(pbm, pyArray1DCon)
print "epetraVector1DCon   = \n", epetraVector1DCon

# 2-D Contiguous:
print "\n-- 2-D Contiguous ----------------------"
pyArray2DCon         = arange(length, typecode=Float64)
pyArray2DCon.shape   = (2,6)
numPyArray2DCon      = Epetra.NumPyArray(pyArray2DCon)
print "numPyArray2DCon     = \n", numPyArray2DCon
epetraVector2DCon    = Epetra.Vector(pbm, pyArray2DCon)
print "epetraVector2DCon   = \n", epetraVector2DCon

# 2-D Non-Contiguous:
print "\n-- 2-D Non-Contiguous ------------------"
pyArray2DNonCon      = swapaxes(pyArray2DCon,0,1)
numPyArray2DNonCon   = Epetra.NumPyArray(pyArray2DNonCon)
print "numPyArray2DNonCon  = \n", numPyArray2DNonCon
epetraVector2DNonCon = Epetra.Vector(pbm, pyArray2DNonCon)
print "epetraVector2DNonCon= \n", epetraVector2DNonCon

# 3-D Contiguous:
print "\n-- 3-D Contiguous ----------------------"
pyArray3DCon         = arange(length, typecode=Float64)
pyArray3DCon.shape   = (2,2,3)
numPyArray3DCon      = Epetra.NumPyArray(pyArray3DCon)
print "numPyArray3DCon     = \n", numPyArray3DCon
epetraVector3DCon    = Epetra.Vector(pbm, pyArray3DCon)
print "epetraVector3DCon   = \n", epetraVector3DCon

# 3-D Non-Contiguous:
print "\n-- 3-D Non-Contiguous ------------------"
pyArray3DNonCon      = swapaxes(pyArray3DCon,0,2)
numPyArray3DNonCon   = Epetra.NumPyArray(pyArray3DNonCon)
print "numPyArray3DNonCon  = \n", numPyArray3DNonCon
epetraVector3DNonCon = Epetra.Vector(pbm, pyArray3DNonCon)
print "epetraVector3DNonCon= \n", epetraVector3DNonCon

# 3-D Non-Contiguous Stride 1:
print "\n-- 3-D Non-Contiguous Stride 1 ---------"
pyArray3DNonCon     = swapaxes(pyArray3DCon,0,1)
numPyArray3DNonCon  = Epetra.NumPyArray(pyArray3DNonCon)
print "numPyArray3DNonCon Stride 1= \n", numPyArray3DNonCon
epetraVector3DNonCon  = Epetra.Vector(pbm, pyArray3DNonCon)
print "epetraVector3DNonCon= \n", epetraVector3DNonCon

length              = 24
pbm = Epetra.BlockMap(length,1,0,comm)

# 4-D Contiguous:
print "\n-- 4-D Contiguous ----------------------"
pyArray4DCon        = arange(length, typecode=Float64)
pyArray4DCon.shape  = (2,2,2,3)
numPyArray4DCon     = Epetra.NumPyArray(pyArray4DCon)
print "numPyArray4DCon     = \n", numPyArray4DCon
epetraVector4DCon    = Epetra.Vector(pbm, pyArray4DCon)
print "epetraVector4DCon   = \n", epetraVector4DCon

# 4-D Non-Contiguous:
print "\n-- 4-D Non-Contiguous ------------------"
pyArray4DNonCon     = swapaxes(pyArray4DCon,0,3)
numPyArray4DNonCon  = Epetra.NumPyArray(pyArray4DNonCon)
print "numPyArray4DNonCon  = \n", numPyArray4DNonCon
epetraVector4DNonCon = Epetra.Vector(pbm, pyArray4DNonCon)
print "epetraVector4DNonCon= \n", epetraVector4DNonCon

# 4-D Non-Contiguous Stride 1:
print "\n-- 4-D Non-Contiguous Stride 1 ---------"
pyArray4DNonCon     = swapaxes(pyArray4DCon,0,1)
numPyArray4DNonCon  = Epetra.NumPyArray(pyArray4DNonCon)
print "numPyArray4DNonCon Stride 1 = \n", numPyArray4DNonCon
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
epetraVector1DCon.unloadViaCopy(pyArray1DCon)
numPyArray1DCon      = Epetra.NumPyArray(pyArray1DCon)
print "numPyArray1DCon     = \n", numPyArray1DCon

# 2-D Contiguous:
print "\n-- 2-D Contiguous ----------------------"
pyArray2DCon         = arange(length, typecode=Float64)
pyArray2DCon.shape   = (2,6)
epetraVector2DCon    = Epetra.Vector(pbm, pyArray2DCon)
epetraVector2DCon.Scale(5.0)
print "epetraVector2DCon   = \n", epetraVector2DCon
epetraVector2DCon.unloadViaCopy(pyArray2DCon)
numPyArray2DCon      = Epetra.NumPyArray(pyArray2DCon)
print "numPyArray2DCon     = \n", numPyArray2DCon

# 2-D Non-Contiguous:
print "\n-- 2-D Non-Contiguous ------------------"
pyArray2DCon         = arange(length, typecode=Float64)
pyArray2DCon.shape   = (2,6)
pyArray2DNonCon      = swapaxes(pyArray2DCon,0,1)
epetraVector2DNonCon = Epetra.Vector(pbm, pyArray2DNonCon)
epetraVector2DNonCon.Scale(7.0)
print "epetraVector2DNonCon= \n", epetraVector2DNonCon
epetraVector2DNonCon.unloadViaCopy(pyArray2DNonCon)
numPyArray2DNonCon   = Epetra.NumPyArray(pyArray2DNonCon)
print "numPyArray2DNonCon  = \n", numPyArray2DNonCon

# 3-D Contiguous:
print "\n-- 3-D Contiguous ----------------------"
pyArray3DCon         = arange(length, typecode=Float64)
pyArray3DCon.shape   = (2,2,3)
epetraVector3DCon    = Epetra.Vector(pbm, pyArray3DCon)
epetraVector3DCon.Scale(11.0)
print "epetraVector3DCon   = \n", epetraVector3DCon
epetraVector3DCon.unloadViaCopy(pyArray3DCon)
numPyArray3DCon      = Epetra.NumPyArray(pyArray3DCon)
print "numPyArray3DCon     = \n", numPyArray3DCon

# 3-D Non-Contiguous:
print "\n-- 3-D Non-Contiguous ------------------"
pyArray3DCon         = arange(length, typecode=Float64)
pyArray3DCon.shape   = (2,2,3)
pyArray3DNonCon      = swapaxes(pyArray3DCon,0,2)
epetraVector3DNonCon = Epetra.Vector(pbm, pyArray3DNonCon)
epetraVector3DNonCon.Scale(13.0)
print "epetraVector3DNonCon= \n", epetraVector3DNonCon
epetraVector3DNonCon.unloadViaCopy(pyArray3DNonCon)
numPyArray3DNonCon   = Epetra.NumPyArray(pyArray3DNonCon)
print "numPyArray3DNonCon  = \n", numPyArray3DNonCon

# 3-D Non-Contiguous Stride 1:
print "\n-- 3-D Non-Contiguous Stride 1 ---------"
pyArray3DCon         = arange(length, typecode=Float64)
pyArray3DCon.shape   = (2,2,3)
pyArray3DNonCon     = swapaxes(pyArray3DCon,0,1)
epetraVector3DNonCon  = Epetra.Vector(pbm, pyArray3DNonCon)
epetraVector3DNonCon.Scale(17.0)
print "epetraVector3DNonCon= \n", epetraVector3DNonCon
epetraVector3DNonCon.unloadViaCopy(pyArray3DNonCon)
numPyArray3DNonCon  = Epetra.NumPyArray(pyArray3DNonCon)
print "numPyArray3DNonCon Stride 1= \n", numPyArray3DNonCon

length              = 24
pbm = Epetra.BlockMap(length,1,0,comm)

# 4-D Contiguous:
print "\n-- 4-D Contiguous ----------------------"
pyArray4DCon        = arange(length, typecode=Float64)
pyArray4DCon.shape  = (2,2,2,3)
epetraVector4DCon    = Epetra.Vector(pbm, pyArray4DCon)
epetraVector4DCon.Scale(19.0)
print "epetraVector4DCon   = \n", epetraVector4DCon
epetraVector4DCon.unloadViaCopy(pyArray4DCon)
numPyArray4DCon     = Epetra.NumPyArray(pyArray4DCon)
print "numPyArray4DCon     = \n", numPyArray4DCon

# 4-D Non-Contiguous:
print "\n-- 4-D Non-Contiguous ------------------"
pyArray4DCon        = arange(length, typecode=Float64)
pyArray4DCon.shape  = (2,2,2,3)
pyArray4DNonCon     = swapaxes(pyArray4DCon,0,3)
epetraVector4DNonCon = Epetra.Vector(pbm, pyArray4DNonCon)
epetraVector4DNonCon.Scale(23.0)
print "epetraVector4DNonCon= \n", epetraVector4DNonCon
epetraVector4DNonCon.unloadViaCopy(pyArray4DNonCon)
numPyArray4DNonCon  = Epetra.NumPyArray(pyArray4DNonCon)
print "numPyArray4DNonCon  = \n", numPyArray4DNonCon

# 4-D Non-Contiguous Stride 1:
print "\n-- 4-D Non-Contiguous Stride 1 ---------"
pyArray4DCon        = arange(length, typecode=Float64)
pyArray4DCon.shape  = (2,2,2,3)
pyArray4DNonCon     = swapaxes(pyArray4DCon,0,1)
epetraVector4DNonCon = Epetra.Vector(pbm, pyArray4DNonCon)
epetraVector4DNonCon.Scale(27.0)
print "epetraVector4DNonCon= \n", epetraVector4DNonCon
epetraVector4DNonCon.unloadViaCopy(pyArray4DNonCon)
numPyArray4DNonCon  = Epetra.NumPyArray(pyArray4DNonCon)
print "numPyArray4DNonCon Stride 1 = \n", numPyArray4DNonCon
