dnl Determine F90 vendor and version string.
AC_DEFUN([WK_FC_GET_VENDOR],
[AC_CACHE_CHECK([the compiler ID],
[wk_cv_prog_f90_version_string],
[$FC -version >conftest.log 2>&1
$FC -V >>conftest.log 2>&1
$FC --version >>conftest.log 2>&1

wk_grep_f90_NAG=`grep NAG conftest.log | head -1`
wk_grep_f90_Compaq=`grep Compaq conftest.log | head -1`
wk_grep_f90_Digital=`grep DIGITAL conftest.log | head -1`
wk_grep_f90_SGI=`grep MIPS conftest.log | head -1`
wk_grep_f90_Intel=`grep 'Intel(R)' conftest.log | head -1`
wk_grep_f90_Sun=`grep 'Sun' conftest.log | head -1`
wk_grep_f90_Lahey=`grep 'Lahey' conftest.log | head -1`
wk_grep_f90_PGI=`grep 'pgf' conftest.log | head -1`
wk_grep_f90_G95=`grep -i 'g95' conftest.log | grep -i 'gcc' | head -1`
wk_grep_f90_GFORTRAN=`grep -i 'GNU Fortran' conftest.log | head -1`
wk_grep_f90_Absoft=`grep -i 'Absoft' conftest.log | head -1`
 
if test -n "$wk_grep_f90_NAG"; then
  wk_cv_prog_f90_type="NAG"
  wk_cv_prog_f90_version_string=$wk_grep_f90_NAG
  wk_cv_prog_f90_version=[`echo $wk_cv_prog_f90_version_string | sed -e 's/.* Release \([0-9][0-9]*\.[0-9][0-9]*.*$\)/\1/'`]
  wk_cv_prog_f90_major_version=[`echo $wk_cv_prog_f90_version | sed -e 's/\([0-9][0-9]*\)\..*/\1/'`]
elif test -n "$wk_grep_f90_Compaq"; then
  wk_cv_prog_f90_type="Compaq"
  wk_cv_prog_f90_version_string=$wk_grep_f90_Compaq
elif test -n "$wk_grep_f90_Digital"; then
  wk_cv_prog_f90_type="DEC"
  wk_cv_prog_f90_version_string=$wk_grep_f90_Digital
elif test -n "$wk_grep_f90_SGI"; then
  wk_cv_prog_f90_type="SGI"
  wk_cv_prog_f90_version_string=$wk_grep_f90_SGI
elif test -n "$wk_grep_f90_Intel"; then
  wk_cv_prog_f90_type="Intel"
  wk_cv_prog_f90_version_string=$wk_grep_f90_Intel
  wk_cv_prog_f90_version=[`echo $wk_cv_prog_f90_version_string | sed -e 's/.* Version \([0-9][0-9]*\.[0-9][0-9]*\) .*/\1/'`]
  wk_cv_prog_f90_major_version=[`echo $wk_cv_prog_f90_version | sed -e 's/\([0-9][0-9]*\)\..*/\1/'`]
elif test -n "$wk_grep_f90_Sun"; then
  wk_cv_prog_f90_type="Sun"
  wk_cv_prog_f90_version_string=$wk_grep_f90_Sun
  wk_cv_prog_f90_version=[`echo $wk_cv_prog_f90_version_string | sed -e 's/.* Fortran 95 \([0-9][0-9]*\.[0-9][0-9]*\) .*/\1/'`]
  wk_cv_prog_f90_major_version=[`echo $wk_cv_prog_f90_version | sed -e 's/\([0-9][0-9]*\)\..*/\1/'`]
elif test -n "$wk_grep_f90_Lahey"; then
  wk_cv_prog_f90_type="Lahey"
  wk_cv_prog_f90_version_string=$wk_grep_f90_Lahey
elif test -n "$wk_grep_f90_PGI"; then
  wk_cv_prog_f90_type="PGI"
  wk_cv_prog_f90_version_string=$wk_grep_f90_PGI
elif test -n "$wk_grep_f90_G95"; then
  wk_cv_prog_f90_type="G95"
  wk_cv_prog_f90_version_string=$wk_grep_f90_G95
elif test -n "$wk_grep_f90_GFORTRAN"; then
  wk_cv_prog_f90_type="GNU"
  wk_cv_prog_f90_version_string=$wk_grep_f90_GFORTRAN
  wk_cv_prog_f90_version=[`echo $wk_cv_prog_f90_version_string | sed -e 's/.*\([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/'`]
  wk_cv_prog_f90_major_version=[`echo $wk_cv_prog_f90_version | sed -e 's/\([0-9][0-9]*\)\..*/\1/'`]
elif test -n "$wk_grep_f90_Absoft"; then
  wk_cv_prog_f90_type="Absoft"
  wk_cv_prog_f90_version_string=$wk_grep_f90_Absoft
else
  wk_cv_prog_f90_type="unknown"
  wk_cv_prog_f90_version_string="unknown"
fi

rm -f conftest.log

])  dnl end AC_CACHE_CHECK

dnl Vendor-specific variables:
AC_CACHE_CHECK([the compiler vendor], [wk_cv_prog_f90_type])

if test -n "$wk_cv_prog_f90_version"; then
  AC_CACHE_CHECK([the compiler version], [wk_cv_prog_f90_version])
else
  wk_cv_prog_f90_version=$wk_cv_prog_f90_version_string
fi

if test -n "$wk_cv_prog_f90_major_version"; then
  AC_CACHE_CHECK([the compiler major version], [wk_cv_prog_f90_major_version])
else
  wk_cv_prog_f90_major_version=$wk_cv_prog_f90_version
fi

FC_VERSION_STRING=$wk_cv_prog_f90_version_string
FC_VENDOR=$wk_cv_prog_f90_type
FC_VERSION=$wk_cv_prog_f90_version
FC_MAJOR_VERSION=$wk_cv_prog_f90_major_version
AC_SUBST(FC_VERSION_STRING)
AC_SUBST(FC_VENDOR)
AC_SUBST(FC_VERSION)
AC_SUBST(FC_MAJOR_VERSION)

dnl Module names: (all compilers apparently have converged to '.mod')
dnl The perl scripts need a quoted version of this
FC_MODNAME='$(1:.o=.mod)'
FC_MODNAME_Q='\$(1:.o=.mod)'
AC_SUBST(FC_MODNAME)
AC_SUBST(FC_MODNAME_Q)

])  dnl end AC_DEFUN


