#! /bin/sh
# -------------------------------------------------------------------------- #
# Simple script to create shared libraries on LINUX/GCC and LAM/MPI.
#
# This is the procedure:
# - Configure trilinos with MPI support, for example using `--enable-mpi'
# - Compile and install using `make install'. This compiles the python 
#   modules as well, but static libraries are used. 
# - Go to the packages installation directory (for example,
#   Trilinos/G4_MPI/, as specified by INSTALL_LIB)
# - Fix the value of INSTALL_LIB below.
# - Execute this script. This creates shared libraries, copies them in
#   the installation subdirectory, then re-compiles the python modules
#   using these shared libraries.
# - You still need to extend you DYLD_LIBRARY_PATH and your PYTHONPATH
#   accondingly to the install directory of Trilinos to use these
#   modules.
#
# Notes:
# - NOX/LOCA are not supported.
# - These packages are supposed to be enabled: 
#   *) teuchos
#   *) epetra
#   *) epetraext,
#   *) triutils
#   *) amesos
#   *) aztecoo
#   *) ifpack
#   *) ml 
#   If other packages have been
#   enabled, they will not be included in the Python interface. In particular,
#   NOX/LOCA are not supported by this script.
# - You MUST have compiled LAM/MPI with support for shared libraries. If
#   you don't, the compilation and installation will run smoothly, but
#   python will generate an exception while using some modules.
#
# Marzio Sala, SNL 9214
# Last updated on 17-Jul-05.
# -------------------------------------------------------------------------- #

INSTALL_LIB=${HOME}/Trilinos/G4_MPI

# ------- #
# Teuchos #
# ------- #
#
cd packages/teuchos/src
g++ -dynamiclib -o libteuchos.dylib *.o \
  -framework vecLib -lmpi -llam -single_module -lf2c
cp libteuchos.dylib $INSTALL_LIB/lib
cd ../../..

# ------ #
# Epetra #
# ------ #
#
cd packages/epetra/src
g++ -dynamiclib -o libepetra.dylib *.o \
  -framework vecLib -lmpi -llam -single_module -lf2c
cp libepetra.dylib $INSTALL_LIB/lib
cd ../python/src
#make clean
#make install
cd ../../../..

# --------- #
# EpetraExt #
# --------- #
#
cd packages/epetraext/src
g++ -dynamiclib -o libepetraext.dylib *.o \
  ~/Trilinos/G4_MPI/packages/epetra/src/libepetra.dylib \
  -framework vecLib -lmpi -llam -single_module -lf2c
cp libepetraext.dylib $INSTALL_LIB/lib
cd ../python/src
make clean
make install
cd ../../../..

# -------- #
# Triutils #
# -------- #
#
cd packages/triutils/src
g++ -dynamiclib -o libtriutils.dylib *.o \
  ~/Trilinos/G4_MPI/packages/epetra/src/libepetra.dylib \
  -framework vecLib -lmpi -llam -single_module -lf2c
cp libtriutils.dylib $INSTALL_LIB/lib
cd ../python/src
make clean
make install
cd ../../../..

# ------- #
# AztecOO #
# ------- #
#
cd packages/aztecoo/src
g++ -dynamiclib -o libaztecoo.dylib *.o \
  ~/Trilinos/G4_MPI/packages/epetra/src/libepetra.dylib \
  ~/Trilinos/G4_MPI/packages/teuchos/src/libteuchos.dylib \
  -framework vecLib -lmpi -llam -single_module -lf2c
cp libaztecoo.dylib $INSTALL_LIB/lib
cd ../python/src
make clean
make install
cd ../../../..

# ------ #
# Amesos #
# ------ #
#
cd packages/amesos/src
g++ -dynamiclib -o libamesos.dylib *.o \
  ~/Trilinos/G4_MPI/packages/epetraext/src/libepetraext.dylib \
  ~/Trilinos/G4_MPI/packages/epetra/src/libepetra.dylib \
  ~/Trilinos/G4_MPI/packages/teuchos/src/libteuchos.dylib \
  -framework vecLib \
  -lmpi -llam -llammpi++ -single_module -lf2c
cp libamesos.dylib $INSTALL_LIB/lib
cd ../python/src
make clean
make install
cd ../../../..

# ------ #
# IFPACK #
# ------ #
#
cd packages/ifpack/src
g++ -dynamiclib -o libifpack.dylib *.o \
  ~/Trilinos/G4_MPI/packages/epetraext/src/libepetraext.dylib \
  ~/Trilinos/G4_MPI/packages/epetra/src/libepetra.dylib \
  ~/Trilinos/G4_MPI/packages/aztecoo/src/libaztecoo.dylib \
  ~/Trilinos/G4_MPI/packages/amesos/src/libamesos.dylib \
  ~/Trilinos/G4_MPI/packages/teuchos/src/libteuchos.dylib \
  -framework vecLib -lmpi -llam -single_module -lf2c
cp libifpack.dylib $INSTALL_LIB/lib
cd ../python/src
make clean
make install
cd ../../../..

# -- #
# ML #
# -- #
#
cd packages/ml/src
g++ -dynamiclib -o libml.dylib *.o \
  ~/Trilinos/G4_MPI/packages/epetra/src/libepetra.dylib \
  ~/Trilinos/G4_MPI/packages/epetraext/src/libepetraext.dylib \
  ~/Trilinos/G4_MPI/packages/triutils/src/libtriutils.dylib \
  ~/Trilinos/G4_MPI/packages/aztecoo/src/libaztecoo.dylib \
  ~/Trilinos/G4_MPI/packages/amesos/src/libamesos.dylib \
  ~/Trilinos/G4_MPI/packages/ifpack/src/libifpack.dylib \
  ~/Trilinos/G4_MPI/packages/teuchos/src/libteuchos.dylib \
  -framework vecLib -lmpi -llam -single_module -lf2c
cp libml.dylib $INSTALL_LIB/lib
cd ../python/src
make clean
make install
cd ../../../..
