dnl @synopsis A_CAN_BE_FIXEDC_CXX_COMPLEX_BLAS_PROBLEM
dnl
dnl @version
dnl @author
dnl
AC_DEFUN([AC_CXX_COMPLEX_BLAS_PROBLEM_CAN_BE_FIXED],
[AC_CACHE_CHECK(whether CDOTC and ZDOTC problems can be fixed,
ac_cv_cxx_complex_blas_problem_fixed,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
 AC_TRY_RUN([
#include  <complex>
extern "C" { 
std::complex<float> F77_BLAS_MANGLE(cdotc,CDOTC)(std::complex<float> *ret, const int* n, const  std::complex<float> x[], const int* incx, const std::complex<float> y[], const int* incy); 
}
int main() {
 const int NUM=2;
 const int INC=1;
 std::complex<float> f[NUM];
 const std::complex<float> ONE = std::complex<float>(1.0,0.0),
                           TWO = std::complex<float>(2.0,0.0);
 f[0] =  ONE; f[1] =  ONE;
 std::complex<float> ret(0.0,0.0);
 F77_BLAS_MANGLE(cdotc,CDOTC)(&ret,&NUM,f,&INC,f,&INC);
 return (ret == TWO ? 0 : 1);
}
],ac_cv_cxx_complex_blas_problem_fixed=yes,ac_cv_cxx_complex_blas_problem_fixed=no)
AC_LANG_RESTORE
LIBS="$save_LIBS"
])
if test "$ac_cv_cxx_complex_blas_problem_fixed" = yes; then
  AC_DEFINE(HAVE_FIXABLE_COMPLEX_BLAS_PROBLEM,,[define if the problem calling BLAS routines CDOTC and ZDOTC from C can be fixed])
fi
])
