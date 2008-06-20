dnl @synopsis TAC_ARG_WITH_BLASLIB
dnl
dnl Test for --with-blaslib="name".
dnl 
dnl Prepends the specified name to the list of files to check for BLAS
dnl routines.  
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_BLASLIB],
[
AC_ARG_WITH(blaslib,
AC_HELP_STRING([--with-blaslib], 
[name of library containing BLAS: will search lib directories for
-lname]),
[
USE_BLASLIB=yes
NEWBLASLIB=${withval}
]
)

   BLASLIBS="cxml blas complib.sgimath"

if test "X${USE_BLASLIB}" = "Xyes"; then

   BLASLIBS="${NEWBLASLIB} ${BLASLIBS}"

fi
])

