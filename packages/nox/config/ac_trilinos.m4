AC_DEFUN([AC_TRILINOS],
[
if test ${USE_TRILINOS} = yes; then

   AC_MSG_CHECKING(whether TRILINOS_HOME is defined)
   if test $TRILINOS_HOME; then
      AC_MSG_RESULT([yes (${TRILINOS_HOME})])
   else
      AC_MSG_RESULT([no (using default of ${HOME}/Trilinos)])
      TRILINOS_HOME=${HOME}/Trilinos
   fi

   AC_MSG_CHECKING(whether TRILINOS_HOME is valid)
   if test -d ${TRILINOS_HOME}; then
      AC_MSG_RESULT(yes)
   else
      AC_MSG_RESULT(no)
      USE_TRILINOS=no
   fi

   AC_MSG_CHECKING(whether TRILINOS_ARCH is defined)
   if test ${TRILINOS_ARCH}; then
      AC_MSG_RESULT([yes (${TRILINOS_ARCH})])
   else
      TRILINOS_ARCH="LINUX"
      AC_MSG_RESULT([no (using default of ${TRILINOS_ARCH})])
   fi

   AC_MSG_CHECKING(whether TRILINOS_COMM is defined)
   if test ${TRILINOS_COMM}; then
      AC_MSG_RESULT([yes (${TRILINOS_COMM})])
   else
      TRILINOS_COMM="SERIAL"
      AC_MSG_RESULT([no (using default of ${TRILINOS_COMM})])
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
      AC_MSG_RESULT([no (using default of ${TRILINOS_TARGET})])
   fi
   

fi
])