#! /usr/bin/env python

import setpath
from PyTrilinos import Epetra

# Create a communicator object
comm = Epetra.SerialComm()
assert(comm.MyPID()   == 0)
assert(comm.NumProc() == 1)
print "comm = {", comm, "}"

# Create a map
pm  = Epetra.Map(4,0,comm)
# assert(pm.Comm()                == comm            )
assert(pm.NumGlobalElements()   == 4               )
assert(pm.NumMyElements()       == 4               )
assert(pm.ElementSize()         == 1               )
assert(pm.IndexBase()           == 0               )
assert(pm.NumGlobalPoints()     == 4               )
assert(pm.NumMyPoints()         == 4               )
assert(pm.MinMyElementSize()    == pm.ElementSize())
assert(pm.MaxMyElementSize()    == pm.ElementSize())
assert(pm.MinElementSize()      == pm.ElementSize())
assert(pm.MaxElementSize()      == pm.ElementSize())
assert(pm.ConstantElementSize() == True            )
assert(pm.SameAs(pm)            == True            )
assert(pm.DistributedGlobal()   == False           )
assert(pm.MinAllGID()           == 0               )
assert(pm.MaxAllGID()           == 3               )
assert(pm.MinMyGID()            == 0               )
assert(pm.MaxMyGID()            == 3               )
assert(pm.MinLID()              == 0               )
assert(pm.MaxLID()              == 3               )
for i in range(pm.NumMyElements()):
    assert(pm.LID(i)   == pm.GID(i))
    assert(pm.MyGID(i) == True     )
    assert(pm.MyLID(i) == True     )
print "pm = {", pm, "}"

# Create a block map
pbm = Epetra.BlockMap(2,2,0,comm)
# assert(pbm.Comm()                == comm             )
assert(pbm.NumGlobalElements()   == 2                )
assert(pbm.NumMyElements()       == 2                )
assert(pbm.ElementSize()         == 2                )
assert(pbm.IndexBase()           == 0                )
assert(pbm.NumGlobalPoints()     == 4                )
assert(pbm.NumMyPoints()         == 4                )
assert(pbm.MinMyElementSize()    == pbm.ElementSize())
assert(pbm.MaxMyElementSize()    == pbm.ElementSize())
assert(pbm.MinElementSize()      == pbm.ElementSize())
assert(pbm.MaxElementSize()      == pbm.ElementSize())
assert(pbm.ConstantElementSize() == True             )
assert(pbm.SameAs(pm)            == False            )
assert(pbm.DistributedGlobal()   == False            )
assert(pbm.MinAllGID()           == 0                )
assert(pbm.MaxAllGID()           == 1                )
assert(pbm.MinMyGID()            == 0                )
assert(pbm.MaxMyGID()            == 1                )
assert(pbm.MinLID()              == 0                )
assert(pbm.MaxLID()              == 1                )
for i in range(pbm.NumMyElements()):
    assert(pbm.LID(i)   == pbm.GID(i))
    assert(pbm.MyGID(i) == True      )
    assert(pbm.MyLID(i) == True      )
print "pbm = {", pbm, "}"

# Create a single vector from a non-blocked map
pv1 = Epetra.Vector(pm)
pv1.PutScalar(1.0)
print "pv1 = {\n", pv1, "}"

# Create a single vector from a blocked map
pv2 = Epetra.Vector(pbm)
pv2.PutScalar(2.0)
print "pv2 = {\n", pv2, "}"

# Create a multi-vector from a non-blocked map
pv3 = Epetra.MultiVector(pm,2)
pv3.PutScalar(3.0)
print "pv3 = {\n", pv3, "}"

# Create a multi-vector from a blocked map
pv4 = Epetra.MultiVector(pbm,2)
pv4.PutScalar(4.0)
print "pv4 = {\n", pv4, "}"

#-------------------------------------------------------------
# NumPy Array to Epetra.Vector testing
#-------------------------------------------------------------

from Numeric import *

# Create a NumPy array
length = 10
# Create a block map for the Epetra Vector
pbm2 = Epetra.BlockMap(length,1,0,comm)

nparray =  arange(length,typecode=Float64)
nparray *= 11.0
print "Created nparray =", nparray

# Create an empty Vector using just the BlockMap
pvector_1 = Epetra.Vector(pbm2)
print "Using nparray, created pvector_1 = {\n", pvector_1, "}"

# Create a Vector using the NumPy array
nparray /= 11.0
print "Scaled nparray by 1/11.0"
pvector_2 = Epetra.Vector(pbm2, nparray)
print "Using nparray, created pvector_2 = {\n", pvector_2, "}"

## # Create a Vector using the NumPy array
## pvector_3 = Epetra.Vector(pbm2, Epetra.PyObjectHolder(nparray))
## print "pvector_3 = {\n", pvector_3, "}"

# Create a Vector using the NumPy array
nparray *= 11.0
print "Scaled nparray by 11.0"
## pvector_2 = Epetra.Vector(pbm2, Epetra.PyObjectHolder(nparray))
pvector_2 = Epetra.Vector(pbm2, nparray)
print "Using nparray, created pvector_2 = {\n", pvector_2, "}"

# Load Vector using loadViaCopy
nparray /= 11.0
print "Scaled nparray by 1/11.0"
## pvector_2.loadViaCopy(Epetra.PyObjectHolder(nparray))
pvector_2.loadViaCopy(nparray)
print "Using nparray, loaded pvector_2 = {\n", pvector_2, "}"

# Unload Vector using loadViaCopy
pvector_2.Scale(11.0) # Mulitply Epetra_Vector by 11.0
print "Scaled pvector_2 in place by 11.0"
print "pvector_2 = {", pvector_2, "}"
## pvector_2.unloadViaCopy(Epetra.PyObjectHolder(nparray))
pvector_2.unloadViaCopy(nparray)
print "Unloaded nparray from modified pvector_2 =", nparray
