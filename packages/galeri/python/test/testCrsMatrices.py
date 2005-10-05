#! /usr/bin/env python

try:
  import setpath
  import Epetra, Galeri
except ImportError:
  from PyTrilinos import Epetra, Galeri
  print "Using system-installed Epetra, Galeri"

Comm = Epetra.PyComm()

nx = 16 * Comm.NumProc()
ny = 16 
nz = 16
mx = Comm.NumProc()
my = 1
mz = 1

# ============================================================================ 

List = {"n": nx }

Map = Galeri.CreateMap("Interlaced", Comm, List)

if Comm.MyPID() == 0:
  print "Creating `Diag'..."
Matrix = Galeri.CreateCrsMatrix("Diag", Map, List)

if Comm.MyPID() == 0:
  print "Creating `Tridiag'..."
Matrix = Galeri.CreateCrsMatrix("Tridiag", Map, List)

if Comm.MyPID() == 0:
  print "Creating `Laplace1D'..."
Matrix = Galeri.CreateCrsMatrix("Laplace1D", Map, List)

if Comm.MyPID() == 0:
  print "Creating `Laplace1DNeumann'..."
Matrix = Galeri.CreateCrsMatrix("Laplace1DNeumann", Map, List)

if Comm.MyPID() == 0:
  print "Creating `Cauchy'..."
Matrix = Galeri.CreateCrsMatrix("Cauchy", Map, List)

if Comm.MyPID() == 0:
  print "Creating `Fielder'..."
Matrix = Galeri.CreateCrsMatrix("Fielder", Map, List)

if Comm.MyPID() == 0:
  print "Creating `Hanowa'..."
Matrix = Galeri.CreateCrsMatrix("Hanowa", Map, List)

if Comm.MyPID() == 0:
  print "Creating `Hilbert'..."
Matrix = Galeri.CreateCrsMatrix("Hilbert", Map, List)

if Comm.MyPID() == 0:
  print "Creating `JordanBlock'..."
Matrix = Galeri.CreateCrsMatrix("JordanBlock", Map, List)

if Comm.MyPID() == 0:
  print "Creating `KMS'..."
Matrix = Galeri.CreateCrsMatrix("KMS", Map, List)

if Comm.MyPID() == 0:
  print "Creating `Lehmer'..."
Matrix = Galeri.CreateCrsMatrix("Lehmer", Map, List)

if Comm.MyPID() == 0:
  print "Creating `Ones'..."
Matrix = Galeri.CreateCrsMatrix("Ones", Map, List)

if Comm.MyPID() == 0:
  print "Creating `Pei'..."
Matrix = Galeri.CreateCrsMatrix("Pei", Map, List)

if Comm.MyPID() == 0:
  print "Creating `Ris'..."
Matrix = Galeri.CreateCrsMatrix("Ris", Map, List)

# ============================================================================ 

List = {"n": nx * ny,
        "nx": nx, "ny": ny,
        "mx": mx, "my": my }
Map = Galeri.CreateMap("Random", Comm, List)

if Comm.MyPID() == 0:
  print "Creating `Cross2D'..."
Matrix = Galeri.CreateCrsMatrix("Cross2D", Map, List)

if Comm.MyPID() == 0:
  print "Creating `BigCross2D'..."
Matrix = Galeri.CreateCrsMatrix("BigCross2D", Map, List)

if Comm.MyPID() == 0:
  print "Creating `Star2D'..."
Matrix = Galeri.CreateCrsMatrix("Star2D", Map, List)

if Comm.MyPID() == 0:
  print "Creating `BigStar2D'..."
Matrix = Galeri.CreateCrsMatrix("BigStar2D", Map, List)

if Comm.MyPID() == 0:
  print "Creating `Laplace2D'..."
Matrix = Galeri.CreateCrsMatrix("Laplace2D", Map, List)

if Comm.MyPID() == 0:
  print "Creating `Stretched2D'..."
Matrix = Galeri.CreateCrsMatrix("Stretched2D", Map, List)

if Comm.MyPID() == 0:
  print "Creating `UniFlow2D'..."
Matrix = Galeri.CreateCrsMatrix("UniFlow2D", Map, List)

if Comm.MyPID() == 0:
  print "Creating `Recirc2D'..."
Matrix = Galeri.CreateCrsMatrix("Recirc2D", Map, List)

if Comm.MyPID() == 0:
  print "Creating `Biharmonic2D'..."
Matrix = Galeri.CreateCrsMatrix("Biharmonic2D", Map, List)

if Comm.MyPID() == 0:
  print "Creating `Laplace2DFourthOrder'..."
Matrix = Galeri.CreateCrsMatrix("Laplace2DFourthOrder", Map, List)

# ============================================================================ 

List = {"n": nx * ny * nz,
        "nx": nx, "ny": ny, "nz": nz,
        "mx": mx, "my": my, "mz": mz }
Map = Galeri.CreateMap("Random", Comm, List)

if Comm.MyPID() == 0:
  print "Creating `Laplace3D'..."
Matrix = Galeri.CreateCrsMatrix("Laplace3D", Map, List)
