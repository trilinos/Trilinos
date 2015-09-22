dnl Enables export makefile specific code
dnl 
dnl The following AM_CONDITIONALS are set for makefiles to access:
dnl USING_EXPORT_MAKEFILES
dnl USING_PERL via TAC_ARG_WITH_PERL
dnl USING_GNUMAKE
dnl
dnl The following AC_DEFINES are set:
dnl HAVE_EXPORT_MAKEFILES
dnl 
dnl the following variables are set:
dnl PERL_EXE for the perl executable via TAC_ARG_WITH_PERL
dnl 
dnl This file was based on tac_arg_enable_feature.m4 by Mike Heroux
dnl @author Roger Pawlowski <rppawlo@sandia.gov>
dnl
AC_DEFUN([TAC_ARG_ENABLE_EXPORT_MAKEFILES],
[
AC_ARG_ENABLE(export-makefiles,
AC_HELP_STRING([--enable-export-makefiles],[Creates export makefiles in the install (prefix) directory.  This option requires perl to be set in your path or defined with --with-perl=<perl executable>. Note that the export makefiles are always created and used in the build directory, but will not be installable without this option to change the paths. (default is $1)]),
ac_cv_use_export_makefiles=$enableval, 
ac_cv_use_export_makefiles=$1)

AC_MSG_CHECKING(whether to build export makefiles)

if test "X$ac_cv_use_export_makefiles" != "Xno"; then

  AC_MSG_RESULT(yes)
  AC_DEFINE([HAVE_EXPORT_MAKEFILES],,[Define if you want to build export makefiles.])

else

  AC_MSG_RESULT(no)

fi

AM_CONDITIONAL(USING_EXPORT_MAKEFILES, test X${ac_cv_use_export_makefiles} = Xyes)

# Check for perl to run scripts (Required dependency)
TAC_ARG_WITH_PERL

if test "X$HAVE_PERL" != "Xyes" && 
   test "X$ac_cv_use_export_makefiles" != "Xno"; then
  AC_MSG_RESULT(no)
  AC_MSG_ERROR([Failed to find the perl executable.  The flag --enable-export-makefiles requires perl to be either in your path or explicitly defined by the flag --with-perl=<executable>.  If you do not require the export makefiles to be installed via 'make install', you can disable the export makefiles with --disable-export-makefiles.])
fi

# Check for using gnumake to clean up link lines via 
# gnumake's "shell" command. Optional dependency.
AC_DEFUN([TAC_ARG_WITH_GNUMAKE],
[
AC_ARG_WITH(gnumake,
AC_HELP_STRING([--with-gnumake],[Gnu's make has special functions we can use to eliminate redundant paths in the build and link lines. Enable this if you use gnu-make to build Trilinos. This requires that perl is in your path or that you have specified the perl executable with --with-perl=<perl executable>.  Configure will check for the existence of the perl executable and quit with an error if it is not found. (default is no)]),
ac_cv_use_gnumake=$withval, ac_cv_use_gnumake=no)

AC_MSG_CHECKING(whether gnumake specific code should be enabled)

if test "X$ac_cv_use_gnumake" != "Xno"; then
  AC_MSG_RESULT(yes)  
  AC_DEFINE([HAVE_GNUMAKE],,[Define if you are using gnumake - this will shorten your link lines.])
else
  AC_MSG_RESULT(no)
fi
AM_CONDITIONAL(USING_GNUMAKE, test "X$ac_cv_use_gnumake" = "Xyes")
])

TAC_ARG_WITH_GNUMAKE

if test "X$HAVE_PERL" != "Xyes" && 
   test "X$ac_cv_use_gnumake" != "Xno"; then
  AC_MSG_RESULT(no)
  AC_MSG_ERROR([The flag --with-gnumake requires perl to be in your path.  The perl executable can alternatively be explicitly defined by the flag --with-perl=<executable>.])
fi

])

