dnl @synopsis TAC_ARG_ENABLE_COMPILER_OPTIONS
dnl
dnl Test for --enable-debug,  --enable-optimize, --enable-profile,
dnl --enable-purify and attempt to set appropriate compiler/loader
dnl flags depending on the value of target
dnl
dnl Use this macro to facilitate definition of options to set for the
dnl compiler and loader.
dnl 
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_ENABLE_COMPILER_OPTIONS],
[
#------------------------------------------------------------------------
# Check for --enable-debug, --enable-opt, --enable-profile, and --enable-purify
#------------------------------------------------------------------------

#TAC_ARG_ENABLE_OPTION(debug, [enable debugging], DEBUG, no)

#TAC_ARG_ENABLE_OPTION(opt,[enable optimization],OPTIMIZE,no)

#TAC_ARG_ENABLE_OPTION(profile,[enable profiling],PROFILE,no)

#TAC_ARG_ENABLE_OPTION(purify,[enable purify],PURIFY,no)

#------------------------------------------------------------------------
# Check if any user overrides were defined in the above section
#------------------------------------------------------------------------
# Likely will permanently remove the concept of a custom platform when 
# another way to turn on debug, opt, profile, and purify is found.
#AC_ARG_WITH(platform,
#AC_HELP_STRING([--with-platform],[Utilize a custom platform] (must specify)),
#ac_cv_use_platform=$withval, ac_cv_use_platform=no)

#AC_MSG_CHECKING(whether to use a custom platform)

#if test "X$ac_cv_use_platform" != "Xno"; then
#  AC_MSG_RESULT(yes)  
  dnl users who define a custom platform need to add a case below
  dnl that will recognize the name of the platform and run the correct m4 macro
  # case $ac_cv_use_platform in
  dnl if your command line includes --with-platform=foo and there is a customized
  dnl version of the 'tac_set_compile_options_template' file called 
  dnl 'tac_set_compile_options_foo.m4' containing a macro named
  dnl 'TAC_SET_COMPILE_OPTIONS_FOO', use the following syntax:
  # foo*) TAC_SET_COMPILE_OPTIONS_FOO;;
  # esac
#else
#  AC_MSG_RESULT(no) 
#  if test "X$ac_cv_use_debug" != "Xno" || test "X$ac_cv_use_opt" != "Xno" ||
#     test "X$ac_cv_use_profile" != "Xno"|| test "X$ac_cv_use_purify" != "Xno"; then
#    if test "X${CXX}" != "XmpiCC"; then
#      BASE_CXX=${CXX}
#    fi
  # Check if user defined a custom platform name
  #otherwise parse info gathered
#  case $target in
#    *linux*)
#      case $BASE_CXX in
# 	 *g++*) TAC_SET_COMPILE_OPTIONS_LINUX_GCC;;
#    esac
#    ;;

#    *solaris*)
#      case $BASE_CXX in
#	*cc* | *CC*) TAC_SET_COMPILE_OPTIONS_SOLARIS_CC;;
#    esac
#    ;;
    # need to add case for cplant

    # need to add case for sgi - and sub-cases for 32 and 64

#  esac

#  fi
#fi
])

