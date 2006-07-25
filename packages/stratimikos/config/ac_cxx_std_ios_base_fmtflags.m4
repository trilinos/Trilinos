dnl @synopsis AC_CXX_STD_IOS_BASE_FMTFLAGS
dnl
dnl If the compiler recognizes std::ios_base::fmtflags as a format type for IO,
dnl define HAVE_STD_IOS_BASE_FMTFLAGS. Note that we assume ANSI header files 
dnl and namespaces, so this test could fail if this assumption is false, but then we
dnl should use a "long" for the format type.
dnl
dnl
AC_DEFUN([AC_CXX_STD_IOS_BASE_FMTFLAGS],
[AC_CACHE_CHECK(whether the compiler recognizes std::ios_base::fmtflags as supported IO type,
ac_cv_cxx_io_base,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
#include <iostream> using namespace std;],
[[
void dummy(ostream& os) {
const std::ios_base::fmtflags oldf = os.setf(ios::scientific,ios::floatfield);
double val = 1.0;
os << val << endl;
return;};
]], ac_cv_cxx_io_base=yes, ac_cv_cxx_io_base=no)

AC_LANG_RESTORE
])
if test "$ac_cv_cxx_io_base" = yes; then
  AC_DEFINE(HAVE_STD_IOS_BASE_FMTFLAGS,,[define if std::ios_base::fmtflags is supported type])
fi
])
