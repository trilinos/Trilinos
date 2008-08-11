nl @synopsis AC_CXX_STDCCBUG
dnl
dnl If the compiler by std::math functions within namespaces, define RAD_NO_USING_STDCC.
dnl
dnl @version $Id$
dnl @author David M. Gay
dnl
AC_DEFUN([AC_CXX_STDCCBUG],
[AC_CACHE_CHECK(whether the compiler has no trouble with std::math functions within namespaces<>,
ac_cv_cxx_stdccbug,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
#include <cmath>
#include <math.h>
namespace Sacado { namespace Rad {
using std::sin;
template<typename Double> class ADvari;
template<typename Double> ADvari<Double>& sin(ADvari<Double>&);
template<typename Double> class ADvari {
 protected:
	Double v;
 public:
	friend ADvari& sin<>(ADvari&);
	ADvari(){}
	ADvari(Double x): v(x) {}
	~ADvari(){}
	};
}}
namespace std { using Sacado::Rad::sin; }
typedef Sacado::Rad::ADvari<double> A;
],[A f, x; f = sin(x); return 0;], ac_cv_cxx_stdccbug=yes, ac_cv_cxx_stdccbug=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_stdccbug" = no; then
  AC_DEFINE(RAD_NO_USING_STDCC,,
            [define if the compiler is consused by std::sin, etc., within namespace Sacado::Rad])
fi
])
