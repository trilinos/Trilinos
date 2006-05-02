#! /bin/sh

# simple script, make a shared version of Epetra, then recompile
# all the packages defined in PACKAGES, clean and re-make the python
# extensions. The script is delicate, it does not check if something
# fails!
#
# I can use most of PyTrilinos after running this. Some problem persists,
# for example, the Amesos interface of ML does not appear to work properly.

BASEDIR=$HOME/Trilinos/LINUX_MPI
cd packages/epetra/src
gcc -fPIC -shared -Wl,-soname,lib$i.so -o lib$i.so \
$BASEDIR/packages/PyTrilinos/src/libpytrilinos.a *.o
cp libepetra.so $BASEDIR/lib/
cd ../../..

PACKAGES="epetraext aztecoo galeri ifpack triutils amesos ml"
for i in $PACKAGES
do
  if test -d packages/$i ; then
    echo "Processing $i..."
    cd packages/$i/src
    /bin/rm lib$i.a
    make
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
