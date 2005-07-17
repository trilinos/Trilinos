#! /bin/sh
# -------------------------------------------------------------------------- #
# Simple script to create shared libraries on LINUX/GCC and LAM/MPI.
# This is the procedure:
#
# - Configure Trilinos with MPI support, for example using `--enable-mpi'
#
# - Compile and install using `make install'. This compiles the python 
#   modules as well, but static libraries are used. 
#
# - Go to the packages installation directory (for example,
#   Trilinos/LINUX_MPI/)
#
# - execute this script. This creates shared libraries for the following
#   packages:
#   *) teuchos
#   *) epetra
#   *) epetraext,
#   *) triutils
#   *) amesos
#   *) aztecoo
#   *) ifpack
#   *) ml 
#   and copies them in the installation subdirectory, then re-compiles 
#   the python modules using these shared libraries. If one of the above
#   packages has not been enabled, it is simply ignored by this script.
#   Note however that PyTrilinos without Epetra and Teuchos has very
#   limited capabilities.
#
# - You still need to extend you LD_LIBRARY_PATH and your PYTHONPATH
#   accondingly to the install directory of Trilinos to use these
#   modules. For example, in BASH, you might need to write:
#   % export PYTHONPATH=$PYTHONPATH:~/Trilinos/LINUX_MPI/lib/python2.3/site-packages/PyTrilinos
#   % export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/Trilinos/LINUX_MPI/lib
#
# - Test your installation by typing at the shell prompt:
#   % python -c "from PyTrilinos import Epetra"
#   If no error shows up, then the Epetra module has been successfully
#   configured. Then repeat the same by replacing Epetra with 
#   the enabled packages reported above.
#
#
# Notes:
#
# - NOX/LOCA are not supported by this script.
#
# - You MUST have compiled LAM/MPI with support for shared libraries. If
#   you don't, the compilation and installation will run smoothly, but
#   python will generate an exception while using some modules.
#
# - I wasn't able to let python run with MPICH, because MPICH adds several
#   parameters to the input line that python doesn't like.
#
# Marzio Sala, SNL 9214
# Last updated on 17-Jul-05.
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
