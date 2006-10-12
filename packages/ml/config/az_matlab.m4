dnl @synopsis AZ_MATLAB_EXEC
dnl @synopsis AZ_MATLAB_ROOT
dnl
dnl @summary New and revised Matlab support.
dnl
dnl This file provides autoconf support for those applications that
dnl want to use matlab. This set of macros help define $MATLAB (the full 
dnl executable path for Matlab), $MEX (the full executable path for mex),
dnl and $MATLAB_ROOT (location of the matlab install directory)
dnl
dnl $MATLAB_USE and $MATLAB_ROOT_USE are automake variable that defines whether Matlab 
dnl support should be included or not in your application.
dnl
dnl The following is an example of how to set up for matlab usage
dnl within your application in your configure.in:
dnl
dnl   AZ_MATLAB_DEFAULT( )
dnl   AZ_MATLAB_EXEC( )               # Optional
dnl   AZ_MATLAB_ROOT( )               # Optional
dnl   # if $MATLAB is not defined, then the following do nothing.
dnl
dnl The AZ_MATLAB_DEFAULT sets the $MATLAB_USE and $MATLAB_ROOT_USE to false. Thereby, 
dnl excluding it if it was optional.
dnl
dnl The AZ_MATLAB_ROOT looks for the optional configure parameters of
dnl --with-matlab-root/--without-matlab and establishes the $MATLAB_ROOT and
dnl $MATLAB_ROOT USE variables accordingly.
dnl
dnl The AZ_MATLAB_EXEC looks for the optional configure parameters of
dnl --with-matlab-exec/--without-matlab and establishes the $MATLAB, $MEX and
dnl $MATLAB_USE variables accordingly.
dnl
dnl
dnl Once that these macros have be run, we can use MATLAB_USE within
dnl the makefile.am file to conditionally add the Matlab support such
dnl as:
dnl
dnl Makefile.am example showing optional inclusion of directories:
dnl
dnl  if MATLAB_USE
dnl  plugins = plugins
dnl  src = src
dnl  else
dnl  plugins =
dnl  src =
dnl  endif
dnl
dnl  SUBDIRS = . $(plugins) $(src)
dnl
dnl Makefile.am example showing optional shared library build:
dnl
dnl  if MATLAB_USE
dnl  lib_LTLIBRARIES        = libElemList.la
dnl  libElemList_la_SOURCES = libElemList.c
dnl  libElemList_la_CFLAGS  = @MATLAB_CSPEC@
dnl  libElemList_la_LDFLAGS = @MATLAB_LSPEC@
dnl  endif
dnl
dnl Makefile.am example showing optional program build:
dnl
dnl  if MATLAB_USE
dnl  bin_PROGRAMS    = runFunc
dnl  runFunc_SOURCES = runFunc.c
dnl  runFunc_CFLAGS  = @MATLAB_CSPEC@
dnl  runFunc_LDFLAGS = @MATLAB_LSPEC@
dnl  endif
dnl
dnl The above compiles the modules only if MATLAB_USE was specified as
dnl true. Also, the else portion of the if was optional.
dnl
dnl @category InstalledPackages
dnl @author Chris Siefert <csiefer@sandia.gov>
dnl @version 2006-06-26
dnl Based on (read: modied heavily from) the python version by:
dnl @author Robert White <kranki@mac.com>
dnl @author Dustin Mitchell <dustin@cs.uchicago.edu>
dnl @license GPLWithACException

# AZ_MATLAB_DEFAULT( )
# -----------------
# Sets the default to not include Matlab support.

AC_DEFUN([AZ_MATLAB_DEFAULT],
[
    az_matlab_use=false
    AM_CONDITIONAL(MATLAB_USE, test x"$az_matlab_use" = x"true")
    AM_CONDITIONAL(MATLAB_ROOT_USE, test x"$az_matlab_use" = x"true")
])


# AZ_MATLAB_EXEC( [path] )
# -----------------------------------------------------------------
# Handles the various --with-matlab-exec commands.
# Input:
#   $1 is the optional search path for the matlab executable if needed
# Ouput:
#   MATLAB_USE (AM_CONDITIONAL) is true if matlab executable found
#   and --with-matlab was requested; otherwise false.
#   $MATLAB_BIN contains the full executable path to matlab if MATLAB_USE
#   is true.
#
# Example:
#   AZ_MATLAB_EXEC( )
#   or
#   AZ_MATLAB_EXEC("/usr/bin")

AC_DEFUN([AZ_MATLAB_EXEC],
[
    AC_ARG_VAR([MATLAB_EXE],[Matlab Executable Path])

    # unless MATLAB_EXE was supplied to us (as a precious variable),
    # see if --with-matlab[=MatlabExecutablePath], --with-matlab,
    # --without-matlab or --with-matlab=no was given.
    if test -z "$MATLAB_EXE"
    then
        AC_MSG_CHECKING(for --with-matlab-exec)
        AC_ARG_WITH(
            matlab-exec,
            AC_HELP_STRING([--with-matlab-exec@<:@=MATLAB_EXE@:>@],
                [absolute path of Matlab executable]
            ),
            [
                if test "$withval" = "yes"
                then
                    # "yes" was specified, but we don't have a path
                    # for the executable.
                    # So, let's searth the PATH Environment Variable.
                    AC_MSG_RESULT(yes)
                    AC_PATH_PROG([MEX],mex,[],$1)
                    AC_PATH_PROG([MATLAB_BIN],matlab,[],$1)
                    if test -z "$MATLAB_EXE"
                    then
                        AC_MSG_ERROR(no path to matlab found)
                    fi
                    az_matlab_use=true
                    AM_CONDITIONAL(MATLAB_USE, test x"$az_matlab_use" = x"true")
                elif test "$withval" = "no"
                then
                    AC_MSG_RESULT(no)
                    az_matlab_use=false
                    AM_CONDITIONAL(MATLAB_USE, test x"$az_matlab_use" = x"true")
                else
                    # $withval must be the executable path then.
                    AC_SUBST([MATLAB_EXE], ["${withval}"])
                    AC_MSG_RESULT($withval)
                    AC_PATH_PROG([MEX],mex,[],$MATLAB_EXE)
                    AC_PATH_PROG([MATLAB_BIN],matlab,[],$MATLAB_EXE)
                    az_matlab_use=true
                    AM_CONDITIONAL(MATLAB_USE, test x"$az_matlab_use" = x"true")
                fi
            ],
            [
                # --with-matlab-exec was not specified.
                AC_MSG_RESULT(no)
                az_matlab_use=false
                AM_CONDITIONAL(MATLAB_USE, test x"$az_matlab_use" = x"true")
            ]
        )
    fi

])






# AZ_MATLAB_ROOT( [path] )
# -----------------------------------------------------------------
# Handles the various --with-matlab-root commands.
# Input:
#   $1 is the optional search path for the matlab root
# Ouput:
#   MATLAB_ROOT_USE (AM_CONDITIONAL) is true if matlab root is present and
#   and --with-matlab-root was requested; otherwise false.
#   $MATLAB_ROOT contains the full executable path to matlab if MATLAB_ROOT_USE
#   is true.
#
# Example:
#   AZ_MATLAB_ROOT( )
#   or
#   AZ_MATLAB_ROOT("/usr/bin")

AC_DEFUN([AZ_MATLAB_ROOT],
[
    AC_ARG_VAR([MATLAB_ROOT],[Matlab Root Directory])

    # unless MATLAB_ROOT was supplied to us (as a precious variable),
    # see if --with-matlab-root[=MatlabExecutablePath], --with-matlab-root,
    # --without-matlab or --with-matlab=no was given.
    if test -z "$MATLAB_ROOT"
    then
        AC_MSG_CHECKING(for --with-matlab-root)
        AC_ARG_WITH(
            matlab-root,
            AC_HELP_STRING([--with-matlab-root@<:@=MATLAB_ROOT@:>@],
                [absolute path name of Matlab root directory]
            ),
            [
                if test "$withval" = "yes"
                then
                    # "yes" was specified, not an actual path.  Well, we'd better hope that $1
                    # has the root path for matlab.  We'll use mexext to check (not to mention 
                    # we need to find it anyway
                    AC_MSG_RESULT(yes)
                    AC_PATH_PROG([MEXEXT],mexext,[],$1/bin)
                    if test -z "$MEXEXT"
                    then
                        AC_MSG_ERROR(matlab root not found)
                    else 
                      MV=`$MEXEXT`
                      AC_SUBST(MEXEXT_VALUE,["$MV"])
                    fi
                    az_matlab_root_use=true
                    AM_CONDITIONAL(MATLAB_ROOT_USE, test x"$az_matlab_root_use" = x"true")
                elif test "$withval" = "no"
                then
                    AC_MSG_RESULT(no)
                    az_matlab_root_use=false
                    AM_CONDITIONAL(MATLAB_ROOT_USE, test x"$az_matlab_root_use" = x"true")
                else
                    # $withval must be the root path then.
                    AC_SUBST([MATLAB_ROOT], ["${withval}"])
                    AC_MSG_RESULT($withval)

                    # Check for mexext
                    AC_PATH_PROG([MEXEXT],mexext,[],$MATLAB_ROOT/bin)
                    if test -z "$MEXEXT"
                    then
                        AC_MSG_ERROR(mexext not found in $MATLAB_ROOT/bin)
                    else
                      MV=`$MEXEXT`
                      AC_SUBST(MEXEXT_VALUE,["$MV"])
                    fi
                    az_matlab_root_use=true
                    AM_CONDITIONAL(MATLAB_ROOT_USE, test x"$az_matlab_root_use" = x"true")
                fi
            ],
            [
                # --with-matlab-root was not specified.
                AC_MSG_RESULT(no)
                az_matlab_use=false
                AM_CONDITIONAL(MATLAB_ROOT_USE, test x"$az_matlab_use" = x"true")
            ]
        )
    fi

])
