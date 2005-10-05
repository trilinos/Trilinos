#! /usr/bin/env python

try:
  import setpath
  import Epetra, Galeri
except ImportError:
  from PyTrilinos import Epetra, Galeri
  print "Using system-installed Epetra, Galeri"

Comm = Epetra.PyComm()

# ===========================================================================
if Comm.MyPID() == 0:
  print "Create `Linear' Map..."
n = 4 * Comm.NumProc()
List = { "n": n }
Map = Galeri.CreateMap("Linear", Comm, List)

# ===========================================================================
if Comm.MyPID() == 0:
  print "Create `Cartesian2D' Map..."
nx = 4
ny = 4 * Comm.NumProc()
mx = 1
my = Comm.NumProc()
List = { "nx": nx, "ny": ny,
         "mx": mx, "my": my }
Map = Galeri.CreateMap("Cartesian2D", Comm, List)

# ===========================================================================
if Comm.MyPID() == 0:
  print "Create `Cartesian3D' Map..."
nx = 4
ny = 4 
nz = 4 * Comm.NumProc()
mx = 1
my = 1
mz = Comm.NumProc()
List = { "nx": nx, "ny": ny, "nz": nz,
         "mx": mx, "my": my, "mz": mz }
Map = Galeri.CreateMap("Cartesian3D", Comm, List)

# ===========================================================================
if Comm.MyPID() == 0:
  print "Create `Random' Map..."
n = 64 * Comm.NumProc()
List = { "n": n }
Map = Galeri.CreateMap("Random", Comm, List)

# ===========================================================================
if Comm.MyPID() == 0:
  print "Create `Interlaced' Map..."
n = 64 * Comm.NumProc()
List = { "n": n }
Map = Galeri.CreateMap("Interlaced", Comm, List)


