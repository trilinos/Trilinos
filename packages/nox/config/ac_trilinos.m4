AC_DEFUN([AC_TRILINOS],
[
AC_ARG_ENABLE(trilinos,
[  --enable-trilinos-arch  Specify arch directory for TRILINOS],
[
case $enableval in
  no)
    USE_TRILINOS_ARCH=no
  ;;
  *)
    USE_TRILINOS=yes
    USE_TRILINOS_ARCH=yes
    TRILINOS_ARCH_DIR=$enableval
  ;;
esac
],
[USE_TRILINOS_ARCH=no]
)

AC_ARG_ENABLE(trilinos,
[  --disable-trilinos      Disable TRILINOS],
[
case $enableval in
  yes)
    USE_TRILINOS=yes
  ;;
  no)
    USE_TRILINOS=no
  ;;
  *)
    AC_MSG_ERROR([Invalid value for --disable-trilinos ($enableval)])
  ;;
esac
],
[USE_TRILINOS=yes]
)

if test ${USE_TRILINOS} = yes && test ${USE_TRILINOS_ARCH} = no; then

   AC_MSG_CHECKING(whether TRILINOS_HOME is defined)
   if test $TRILINOS_HOME; then
      AC_MSG_RESULT([yes (${TRILINOS_HOME})])
   else
      AC_MSG_RESULT([no (using ${HOME}/Trilinos)])
      TRILINOS_HOME=${HOME}/Trilinos
   fi

   AC_MSG_CHECKING(whether TRILINOS_HOME is valid)
   if test -d ${TRILINOS_HOME}; then
      AC_MSG_RESULT(yes)
   else
      AC_MSG_RESULT([no (disabling Trilinos)])
      USE_TRILINOS=no
   fi
fi

if test ${USE_TRILINOS} = yes && test ${USE_TRILINOS_ARCH} = no; then

   AC_MSG_CHECKING(whether TRILINOS_ARCH is defined)
   if test ${TRILINOS_ARCH}; then
      AC_MSG_RESULT([yes (${TRILINOS_ARCH})])
   else
      case $target in
	rs6000-ibm-aix3.2.* | rs6000-ibm-aix4.* | powerpc-ibm-aix4.*)
	TRILINOS_ARCH="IBMSP"
	;;	
	*)
	TRILINOS_ARCH="LINUX"
  	;;
      esac
      AC_MSG_RESULT([no (using ${TRILINOS_ARCH})])
   fi

   AC_MSG_CHECKING(whether TRILINOS_COMM is defined)
   if test ${USE_MPI} = yes; then
      TRILINOS_COMM="MPI"
   else
      TRILINOS_COMM="SERIAL"
   fi
   AC_MSG_RESULT([overriding with ${TRILINOS_COMM}])

   AC_MSG_CHECKING(whether TRILINOS_ID is defined)
   if test ${TRILINOS_ID}; then
      AC_MSG_RESULT([yes (${TRILINOS_ID})])
   else
      AC_MSG_RESULT(no)
   fi

   AC_MSG_CHECKING(whether TRILINOS_TARGET is defined)
   if test ${TRILINOS_TARGET}; then
      AC_MSG_RESULT([yes (${TRILINOS_TARGET})])
   else
      TRILINOS_TARGET="${TRILINOS_ARCH}.${TRILINOS_COMM}${TRILINOS_ID}"
      AC_MSG_RESULT([no (using ${TRILINOS_TARGET})])
   fi

   TRILINOS_CXXFLAGS="-DEPETRA_${TRILINOS_COMM} -I${TRILINOS_HOME}/packages/epetra/src"
   TRILINOS_LDADD="-L${TRILINOS_HOME}/lib/${TRILINOS_TARGET}"

fi


if test ${USE_TRILINOS} = yes && test ${USE_TRILINOS_ARCH} = yes; then

   AC_MSG_CHECKING(whether TRILINOS_ARCH_DIR is valid)
   if test -d ${TRILINOS_ARCH_DIR}; then
      AC_MSG_RESULT(yes)
   else
      AC_MSG_RESULT([no (disabling Trilinos)])
      USE_TRILINOS=no
   fi

fi

if test ${USE_TRILINOS} = yes && test ${USE_TRILINOS_ARCH} = yes; then

   if test ${USE_MPI} = yes; then
      TRILINOS_COMM="MPI"
   else
      TRILINOS_COMM="SERIAL"
   fi

   TRILINOS_CXXFLAGS="-DEPETRA_${TRILINOS_COMM} -I${TRILINOS_ARCH_DIR}/include/epetra/"
   TRILINOS_LDADD="-L${TRILINOS_ARCH_DIR}/lib/"

fi


])

AC_SUBST(TRILINOS_CXXFLAGS)
AC_SUBST(TRILINOS_LDADD)

