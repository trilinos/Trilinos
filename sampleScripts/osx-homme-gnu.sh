#!/bin/sh                                                                                                                                                                                          
                                                                                                                                                                                                   
if [ ! $TRILINOS_SRC_DIR ]; then
  echo Must specify location of Trilinos source in variable TRILINOS_SRC_DIR
  exit 1
fi

if   [ ! $TRILINOS_INSTALL_DIR -a $Trilinos_DIR ]; then
  echo Setting TRILINOS_INSTALL_DIR to $Trilinos_DIR
  TRILINOS_INSTALL_DIR="$Trilinos_DIR"
elif [ ! $TRILINOS_INSTALL_DIR -a ! $Trilinos_DIR]; then
  echo Must set TRILINOS_INSTALL_DIR.
  exit 1
fi

EXTRA_ARGS=$@

cmake \
  -D CMAKE_INSTALL_PREFIX:PATH="$TRILINOS_INSTALL_DIR" \
  -D CMAKE_BUILD_TYPE=RELEASE \
  -D Trilinos_ENABLE_Epetra:BOOL=ON \
  -D Trilinos_ENABLE_EpetraExt=ON \
  -D Trilinos_ENABLE_NOX:BOOL=ON \
  -D Trilinos_ENABLE_Stratimikos=ON \
  -D Trilinos_ENABLE_Belos=ON \
  -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON \
  -D Trilinos_ENABLE_INSTALL_CMAKE_CONFIG_FILES:BOOL=ON \
  -D TPL_ENABLE_MPI:BOOL=ON \
  -D MPI_BASE_DIR=/opt/openmpi/current \
  $* \
  ${TRILINOS_SRC_DIR}

# \
#   -D BLAS_LIBRARY_NAMES="sci_gnu" \
#   -D LAPACK_LIBRARY_NAMES="sci_gnu" \
#   -D BLAS_LIBRARY_DIRS=$CRAY_LIBSCI_PREFIX_DIR/lib \
#   -D LAPACK_LIBRARY_DIRS=$CRAY_LIBSCI_PREFIX_DIR/lib \
