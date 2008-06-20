dnl @synopsis TAC_ARG_CHECK_MPI
dnl
dnl Check to make sure any definitions set in TAC_ARG_CONFIG_MPI
dnl are valid, set the MPI flags.  Test MPI compile using C++ compiler.
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_CHECK_MPI],
[

if test "X${HAVE_PKG_MPI}" = "Xyes"; then

  if test -n "${MPI_DIR}" && test -z "${MPI_INC}"; then
    MPI_INC="${MPI_DIR}/include"
  fi

  if test -n "${MPI_INC}"; then
    CPPFLAGS="${CPPFLAGS} -I${MPI_INC}"
  fi

  AC_LANG_CPLUSPLUS 
  AC_MSG_CHECKING(for mpi.h)
  AC_TRY_CPP([#include "mpi.h"],
    [AC_MSG_RESULT(yes)], 
    [
     AC_MSG_RESULT(no)  
     echo "-----"
     echo "Cannot link simple MPI program."
     echo "Try --with-mpi-compilers to specify MPI compilers."
     echo "Or try --with-mpi-libs, --with-mpi-incdir, --with-mpi-libdir"
     echo "to specify all the specific MPI compile options."
     echo "-----"
     AC_MSG_ERROR(MPI cannot link)
    ])

  if test -n "${MPI_DIR}" && test -z "${MPI_LIBDIR}"; then
    MPI_LIBDIR="${MPI_DIR}/lib"
  fi

  if test -n "${MPI_LIBDIR}"; then
    LDFLAGS="${LDFLAGS} -L${MPI_LIBDIR}"
  fi

  if test -z "${MPI_LIBS}" && test -n "${MPI_LIBDIR}"; then
    MPI_LIBS="-lmpi"
  fi

  if test -n "${MPI_LIBS}"; then
    LIBS="${MPI_LIBS} ${LIBS}"
  fi

#   AC_LANG_CPLUSPLUS 
#   AC_MSG_CHECKING(whether MPI will link using C++ compiler)
#   AC_TRY_LINK([#include <mpi.h>],
#   [int c; char** v; MPI_Init(&c,&v);],
#   [AC_MSG_RESULT(yes)], 
#   [AC_MSG_RESULT(no)  
#    echo "-----"
#    echo "Cannot link simple MPI program."
#    echo "Try --with-mpi-cxx to specify MPI C++ compile script."
#    echo "Or try --with-mpi-libs, --with-mpi-incdir, --with-mpi-libdir"
#    echo "to specify all the specific MPI compile options."
#    echo "-----"
#    AC_MSG_ERROR(MPI cannot link)]
#   )

fi
])
