#!/ bin/sh
# -------------------------------------------------------------------------- #
# Simple script to create shared libraries on LINUX/GCC and LAM/MPI.
#
# This is the procedure:
# - configure trilions with --enable-mpi
# - compile. This compiles the python modules as well, but static libraries
#   are used.
# - go to the packages installation directory (for example,
#   Trilinos/LINUX_MPI/). 
# - execute this script. This creates shared libraries, copies them in
#   the installation subdirectory, then re-compiles the python modules
#   using these shared libraries.
# - You still need to extend you LD_LIBRARY_PATH and your PYTHONPATH
#   accondingly to the install directory of Trilinos to use these
#   modules.
#
# Notes:
# - NOX/LOCA are not supported.
# - You MUST have compiled LAM/MPI with support for shared libraries. If
#   you don't, the compilation and installation will run smoothly, but
#   python will generate an exception while using some modules.
# - I wasn't able to let python run with MPICH, because MPICH adds several
#   parameters to the input line that python doesn't like.
#
# MS 
# Last updated on 11-Jun-05.
# -------------------------------------------------------------------------- #

LIBS="teuchos epetra epetraext triutils amesos ifpack aztecoo ml"

for i in ${LIBS}
do
  if test -d packages/$i ; then
    echo "Processing $i..."
    cd packages/$i/src
    gcc -shared -Wl,-soname,lib$i.so -o lib$i.so *.o
    cp lib$i.so ../../../lib
    if test -d ../python/src ; then
      cd ../python/src
      make clean
      make install
      cd -
    fi
    cd ../../..
  fi
done
