
#ifndef F77_FUNC
#  if (defined(__IBMC__) || defined(__IBMCPP__))
#    define F77_FUNC(name,NAME) name
#  elif (defined(WIN32) || defined(ICL))
#    define F77_FUNC(name,NAME) NAME
#  else
#    define F77_FUNC(name,NAME) name ## _
#  endif
#endif

/* As F77_FUNC, but for C identifiers containing underscores. */
#ifndef F77_FUNC_
#  if (defined(__IBMC__) || defined(__IBMCPP__))
#    define F77_FUNC_(name,NAME) name
#  elif (defined(WIN32) || defined(ICL))
#    define F77_FUNC_(name,NAME) NAME
#  else
#    define F77_FUNC_(name,NAME) name ## _
#  endif
#endif

#ifndef FC_FUNC
#  if (defined(__IBMC__) || defined(__IBMCPP__))
#    define FC_FUNC(name,NAME) name
#  elif (defined(WIN32) || defined(ICL))
#    define FC_FUNC(name,NAME) NAME
#  else
#    define FC_FUNC(name,NAME) name ## _
#  endif
#endif

/* As FC_FUNC, but for C identifiers containing underscores. */
#ifndef FC_FUNC_
#  if (defined(__IBMC__) || defined(__IBMCPP__))
#    define FC_FUNC_(name,NAME) name
#  elif (defined(WIN32) || defined(ICL))
#    define FC_FUNC_(name,NAME) NAME
#  else
#    define FC_FUNC_(name,NAME) name ## _
#  endif
#endif

/* Define the Fortran name mangling to be used for the BLAS */
#define F77_BLAS_MANGLE(name,NAME) F77_FUNC(name,NAME)


 /* Mangling for Fortran global symbols without underscores. */
#ifndef FortranCInterface_GLOBAL
#  if ! (defined(__IBMC__) || defined(__IBMCPP__))
#    define FortranCInterface_GLOBAL(name,NAME) name
#  elif (defined(WIN32) || defined(ICL))
#    define FortranCInterface_GLOBAL(name,NAME) NAME
#  else
#    define FortranCInterface_GLOBAL(name,NAME) name##_
#  endif
#endif
 
 /* Mangling for Fortran global symbols with underscores. */
#ifndef FortranCInterface_GLOBAL_
#  if ! (defined(__IBMC__) || defined(__IBMCPP__))
#    define FortranCInterface_GLOBAL_(name,NAME) name
#  elif (defined(WIN32) || defined(ICL))
#    define FortranCInterface_GLOBAL_(name,NAME) NAME
#  else
#    define FortranCInterface_GLOBAL_(name,NAME) name##_
#  endif
#endif
 
 /* Mangling for Fortran module symbols without underscores. */
#ifndef FortranCInterface_MODULE
#  if (defined(__IBMC__) || defined(__IBMCPP__))
#    define FortranCInterface_MODULE(mod_name,name, mod_NAME,NAME) __##mod_name##_NMOD_##name
#  elif (defined(WIN32) || defined(ICL))
#    define FortranCInterface_MODULE(mod_name,name, mod_NAME,NAME) mod_NAME##_mp_MOD_##NAME
#  else
#    define FortranCInterface_MODULE(mod_name,name, mod_NAME,NAME) __##mod_name##_MOD_##name
#  endif
#endif
 
 /* Mangling for Fortran module symbols with underscores. */
#ifndef FortranCInterface_MODULE_
#  if (defined(__IBMC__) || defined(__IBMCPP__))
#    define FortranCInterface_MODULE_(mod_name,name, mod_NAME,NAME) __##mod_name##_NMOD_##name
#  elif (defined(WIN32) || defined(ICL))
#    define FortranCInterface_MODULE_(mod_name,name, mod_NAME,NAME) mod_NAME##_mp_MOD_##NAME
#  else
#    define FortranCInterface_MODULE_(mod_name,name, mod_NAME,NAME) __##mod_name##_MOD_##name
#  endif
#endif

#if defined(__IBMC__) || defined(__IBMCPP__)
#  ifndef TEMPLATE_FRIENDS_NOT_SUPPORTED
#    define TEMPLATE_FRIENDS_NOT_SUPPORTED
#  endif
#  ifndef TEUCHOS_PRIVIATE_DELETE_NOT_SUPPORTED
#    define TEUCHOS_PRIVIATE_DELETE_NOT_SUPPORTED
#  endif
#endif

#ifdef _AIX
#  define HAS_C99_TR1_CMATH
#endif

#if (defined(__GNUC__) && !defined(__INTEL_COMPILER)) || defined(__PGI)
#  define HAVE_RTOP_EXPLICIT_INSTANTIATION
#  define HAVE_THYRA_EXPLICIT_INSTANTIATION
#endif

#if ! (defined(WIN32) || defined(ICL))
#  ifndef __APPLE__
#    define FINITE_VALUE_HAVE_GLOBAL_ISINF
#    define FINITE_VALUE_HAVE_GLOBAL_ISNAN
#  endif
#  define FINITE_VALUE_HAVE_STD_ISINF
#  define FINITE_VALUE_HAVE_STD_ISNAN
#endif
