dnl @synopsis AC_CXX_TEMPLATE_QUALIFIER
dnl
dnl If the compiler accepts "template" as a qualifier for templated methods
dnl where the template type cannot be deduced by the argument list,
dnl define HAVE_TEMPLATE_QUALIFIER.
dnl
AC_DEFUN([AC_CXX_TEMPLATE_QUALIFIER],
[AC_CACHE_CHECK(whether the compiler accepts the "template" qualifier,
ac_cv_cxx_template_qualifier,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
class DummyClass {
  public:
    DummyClass() : dummyint( 1 ) {}
    ~DummyClass() {}
    template<typename T>
    T get() { return( (T) dummyint ); }
  private:
    int dummyint;
};
],
[  
   DummyClass my_dummy;
   double temp = my_dummy.template get<double>(); 
],
 ac_cv_cxx_template_qualifier=yes, ac_cv_cxx_template_qualifier=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_template_qualifier" = yes; then
  AC_DEFINE(HAVE_TEMPLATE_QUALIFIER,,[define if the compiler accepts the "template" qualifier])
fi
])
