AC_DEFUN([AC_TRILINOS],
[
AC_ARG_ENABLE(trilinos,
[  --disable-trilinos           Do not use TRILINOS],
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

if test ${USE_TRILINOS} = yes; then

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

if test ${USE_TRILINOS} = yes; then
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
   if test ${TRILINOS_COMM}; then
      AC_MSG_RESULT([yes (${TRILINOS_COMM})])
   else
      TRILINOS_COMM="SERIAL"
      AC_MSG_RESULT([no (using ${TRILINOS_COMM})])
   fi

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
fi
])

AC_SUBST(TRILINOS_HOME)
AC_SUBST(TRILINOS_ARCH)
AC_SUBST(TRILINOS_COMM)
AC_SUBST(TRILINOS_ID)
AC_SUBST(TRILINOS_TARGET)
