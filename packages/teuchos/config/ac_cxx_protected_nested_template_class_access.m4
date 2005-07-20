dnl @synopsis AC_CXX_PROTECTED_NESTED_TEMPLATE_CLASS_ACCESS
dnl
dnl If the compiler supports reinterpret_cast<>, define HAVE_PROTECTED_NESTED_TEMPLATE_CLASS_ACCESS.
dnl
dnl @version $Id$
dnl @author Luc Maisonobe
dnl
AC_DEFUN([AC_CXX_PROTECTED_NESTED_TEMPLATE_CLASS_ACCESS],
[AC_CACHE_CHECK(whether the compiler supports access of protected templated nested classes in derived classes,
ac_cv_cxx_protected_nested_template_class_access,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([class BaseBase { protected: template<class T> class Boo {}; };
template<class T> class Base { public: };
template<class T> class Derived { private: BaseBase::Boo<T> boo_; };],[Derived<double> d; return 0;],
 ac_cv_cxx_protected_nested_template_class_access=yes, ac_cv_cxx_protected_nested_template_class_access=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_protected_nested_template_class_access" = yes; then
  AC_DEFINE(HAVE_PROTECTED_NESTED_TEMPLATE_CLASS_ACCESS,,
            [define if the compiler supports access of protected templated nested classes in derived classes])
fi
])
