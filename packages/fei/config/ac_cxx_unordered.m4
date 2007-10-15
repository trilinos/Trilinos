dnl @synopsis AC_CXX_UNORDERED
dnl
dnl If the compiler provides unordered map and set (currently known
dnl as hash_map and hash_set, but will probably be called unordered_map
dnl and unordered_set if/when they make it into the c++ standard),
dnl define HAVE_UNORDERED.
dnl
dnl @version $Id$
dnl @author Alan Williams
dnl
AC_DEFUN([AC_CXX_UNORDERED],
[AC_CACHE_CHECK(whether the compiler provides unordered associative containers,
ac_cv_cxx_unordered,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
#include <ext/hash_set>
#include <ext/hash_map>
],[__gnu_cxx::hash_set<int> hs; __gnu_cxx::hash_map<int,int> hm;],
 ac_cv_cxx_unordered=yes, ac_cv_cxx_unordered=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_unordered" = yes; then
  AC_DEFINE(HAVE_UNORDERED,,[define if the compiler provides unordered associative containers])
fi
])
