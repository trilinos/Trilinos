dnl @synopsis AC_CXX_STL_MAP_KEY_CONST
dnl
dnl Test if stl map key must be const
dnl
dnl @version $Id$
dnl @author Robert Hoekstra
dnl
AC_DEFUN([AC_CXX_STL_MAP_KEY_CONST],
[AC_CACHE_CHECK(whether stl map key must be const,
ac_cv_cxx_stl_map_key_const,
[AC_REQUIRE([AC_CXX_NAMESPACES])
 AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <map>
#ifdef HAVE_NAMESPACES
using namespace std;
#endif],
[map<int,int> x;
x.insert(pair<int,int>(0,0));
multimap<int,int> y;
y.insert(pair<int,int>(0,0));
return 0;],
ac_cv_cxx_stl_map_key_const=no,
ac_cv_cxx_stl_map_key_const=yes)
AC_LANG_RESTORE
])
if test "$ac_cv_cxx_stl_map_key_const" = yes; then
  AC_DEFINE(MUST_CONST_STL_MAP_KEY,,[define if STL map key is required to be const])
fi
])
