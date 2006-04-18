AC_DEFUN([TAC_ENABLE_PYTHON],
[
AZ_PYTHON_DEFAULT( )  # Default to no python support
AZ_PYTHON_ENABLE( )   # Check for --enable-python, --disable-python
AZ_PYTHON_WITH( )     # Check for --with-python, --without-python

AC_MSG_CHECKING(whether we should build the python wrappers)
if test -n "$PYTHON"; then
  BUILD_PYTHON=yes
  AC_MSG_RESULT(yes)
  AM_CONDITIONAL(BUILD_PYTHON, true)

  # Ensure that we have python version 2.2 or greater (for distutils)
  AZ_PYTHON_VERSION_ENSURE( [2.2] )

  # Python compiler and linker flags
  AZ_PYTHON_CSPEC
  AZ_PYTHON_LSPEC

  # Check that Python.h is available
  save_CPPFLAGS=$CPPFLAGS
  CPPFLAGS="$save_CPPFLAGS $PYTHON_CSPEC"
  AC_LANG([C++])
  AC_CHECK_HEADER(
  [Python.h],
  break,
  AC_MSG_ERROR([You must have Python.h in order to build the Python support!!]))
  CPPFLAGS="$save_CPPFLAGS"

  # If user specifies prefix, use it for the PYTHON_PREFIX
  if test "$prefix" != "$ac_default_prefix"; then
    PYTHON_PREFIX=$prefix
  fi
  if test "$exec_prefix" != "$ac_default_prefix"; then
    PYTHON_EXECPREFIX=$exec_prefix
  fi

else
  AC_MSG_RESULT(no)
  BUILD_PYTHON=no
  AM_CONDITIONAL(BUILD_PYTHON, false)
fi

# ------------------------------------------------------------------------
# If the python wrappers are to be built, then SWIG (Simple Wrapper
# Interface Generator) is required
# ------------------------------------------------------------------------

if test -n "$PYTHON"; then

  # Check for --with-swig[=path]
  AC_MSG_CHECKING(for --with-swig)
  AC_ARG_WITH(swig,
              [AC_HELP_STRING([--with-swig@<:@=SWIG@:>@],
                              [enable swig and set swig binary])],
              [AC_MSG_RESULT(yes)
               WITH_SWIG=yes
               if test X${withval} != Xyes ; then
                 SWIG=$withval
               fi],
              [AC_MSG_RESULT(no)
               AC_CHECK_PROG(WITH_SWIG,swig,yes,no)])

  # Report error if no swig found
  if test ${WITH_SWIG} = no; then
     AC_MSG_ERROR(
     [Python wrappers require swig (Simple Wrapper Interface Generator).
      See http://www.swig.org])
  fi

  # SWIG configuration
  AC_PROG_SWIG(1.3.23)
  SWIG_ENABLE_CXX
  SWIG_MULTI_MODULE_SUPPORT
  SWIG_PYTHON
fi
AM_CONDITIONAL(HAVE_SWIG,test X${WITH_SWIG} = Xyes)
AC_SUBST(GNU_HAVE_SWIG, ${WITH_SWIG})

])
