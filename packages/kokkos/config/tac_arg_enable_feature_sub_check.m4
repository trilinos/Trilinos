dnl @synopsis TAC_ARG_ENABLE_FEATURE_SUB_CHECK(FEATURE_NAME, SUB_FEATURE_NAME, FEATURE_DESCRIPTION, HAVE_NAME)
dnl
dnl This hack gets around the fact that TAC_ARG_ENABLE_FEATURE does not support underscores
dnl in its feature names.  TAC_ARG_ENABLE_FEATURE_SUB_CHECK allows exactly one underscore.  Not great,
dnl but arguably better than supporting no underscores.
dnl
dnl TAC_ARG_ENABLE_FEATURE(feature-sub, [Configure and build feature-sub], FEATURE_SUB, yes) 
dnl   fails because tac_arg_enable_feature tests for ac_cv_use_feature-sub which gets 
dnl   rejected because the `-' is not allowed in variables.  (AC_ARG_ENABLE sets ac_cv_use_feature_sub
dnl   to avoid this problem.)  Use:
dnl 
dnl TAC_ARG_ENABLE_FEATURE_SUB_CHECK(feature, sub, [Configure and build feature-sub], FEATURE_SUB) 
dnl   instead.
dnl
dnl This macro will test for --enable-${FEATURE_NAME}-${SUB_FEATURE_NAME} when configure is run.  
dnl If it is defined and not set to "no" or not defined and --disable-${SUB_FEATURE_NAME} is not
dnl specified then HAVE_${HAVE_NAME} will be defined.
dnl
dnl *NOTE: This macro is designed for the use-case when there is an individual Trilinos package 
dnl offering fine-grained control of a Trilinos option.  This way, the individual package 
dnl option is enabled, as long as the Trilinos option is not disabled.  If the Trilinos option is
dnl disabled, then the user must enable each packages option individually.  For instance:
dnl
dnl --disable-tests --enable-teuchos-tests
dnl
dnl *NOTE: epetra, aztecoo, komplex, ifpack, and other software found in
dnl subdirectories of Trilinos/packages are "packages" in their own right.
dnl However, these packages are also "features" of the larger package
dnl "Trilinos".  Therefore, when configuring from the Trilinos directory,
dnl it is appropriate to refer to these software packages as "features".
dnl
dnl This file was based on tac_arg_enable_package.m4 by Jim Willenbring
dnl and tac_arg_enable_package_sub.m4 by Ken Stanley.
dnl
dnl @author Heidi Thornquist <hkthorn@sandia.gov>
dnl
AC_DEFUN([TAC_ARG_ENABLE_FEATURE_SUB_CHECK],
[
AC_ARG_ENABLE([$2],, ac_cv_use_$2=$enableval, ac_cv_use_$2=yes)

AC_ARG_ENABLE([$1-$2],
AC_HELP_STRING([--enable-$1-$2],[$3 (default is yes if --disable-$2 is not specified)]),
ac_cv_use_$1_$2=$enableval, ac_cv_use_$1_$2=${ac_cv_use_$2})

AC_MSG_CHECKING(whether to use [$1-$2])

if test "X$ac_cv_use_$1_$2" != "Xno"; then
  AC_MSG_RESULT(yes)  
  AC_DEFINE([HAVE_$4],,[Define if want to build $1-$2])
else
  AC_MSG_RESULT(no)
fi
])

