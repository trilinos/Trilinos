dnl @synopsis AC_CXX_INVALID_TEMPLATE_QUALIFIER
dnl
dnl Determine if the compiler requires "template" as a qualifier for templated 
dnl methods called from non-templated code where the template type cannot be 
dnl deduced by the argument list, even though the use of the "template" 
dnl qualifier in this case is invalid according to the ANSI C++ standard.
dnl Defines INVALID_TEMPLATE_QUALIFIER to "template" in this case, empty
dnl otherwise.
dnl
AC_DEFUN([AC_CXX_INVALID_TEMPLATE_QUALIFIER],
[AC_CACHE_CHECK(whether the compiler requires invalid use of the "template" qualifier,
ac_cv_cxx_invalid_template_qualifier,
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
   DummyClass *my_dummy_ptr = &my_dummy;
   double temp = (*my_dummy_ptr).get<double>(); 
],
 ac_cv_cxx_invalid_template_qualifier=no, ac_cv_cxx_invalid_template_qualifier=yes)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_invalid_template_qualifier" = yes; then
  AC_DEFINE(INVALID_TEMPLATE_QUALIFIER,[template],[template qualifier required for calling template methods from non-template code])
else
  AC_DEFINE(INVALID_TEMPLATE_QUALIFIER,[],[template qualifier required for calling template methods from non-template code])
fi
])
