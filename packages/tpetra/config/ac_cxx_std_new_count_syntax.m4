dnl @synopsis AC_CXX_STD_NEW_COUNT_SYNTAX
dnl
dnl If the compiler recognizes std::ios_base::fmtflags as a format type for IO,
dnl define HAVE_STD_NEW_COUNT_SYNTAX.
dnl
dnl
AC_DEFUN([AC_CXX_STD_NEW_COUNT_SYNTAX],
[AC_CACHE_CHECK(whether the compiler recognizes the new form of std::count,
ac_cv_cxx_new_count,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
#ifdef HAVE_VECTOR
#include <vector>
#elif HAVE_VECTOR_H
#include <vector.h>
#endif
#ifdef HAVE_ALGORITHM
#include <algorithm>
#elif HAVE_ALGO_H
#include <algo.h>
#endif
void dummy() {
	std::vector<int> foo;
	std::count(foo.begin(), foo.end(), 0);
}],
[ dummy(); ],
ac_cv_cxx_new_count=yes, ac_cv_cxx_new_count=no)
AC_LANG_RESTORE
])
if test "$ac_cv_cxx_new_count" = yes; then
  AC_DEFINE(HAVE_STD_NEW_COUNT_SYNTAX,,[define if new form of std::count is supported])
fi
])
