
/* #undef TEUCHOS_STANDALONE_PACKAGE */

/* #undef BUILD_SHARED_LIBS */

#define HAVE_TEUCHOS_PARAMETERLIST

/* Define the Fortran name mangling to be used for the BLAS */
#ifndef F77_BLAS_MANGLE
 #define F77_BLAS_MANGLE(name,NAME) name ## _
#endif

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
#ifndef F77_DUMMY_MAIN
/* #undef F77_DUMMY_MAIN */
#endif

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#ifndef F77_FUNC
 #define F77_FUNC(name,NAME) name ## _
#endif

/* As F77_FUNC, but for C identifiers containing underscores. */
#ifndef F77_FUNC_
 #define F77_FUNC_(name,NAME) name ## _
#endif

/* Define if F77 and FC dummy `main' functions are identical. */
#ifndef FC_DUMMY_MAIN_EQ_F77
/* #undef FC_DUMMY_MAIN_EQ_F77 */
#endif

/* Define to 1 if you have the <fpu_control.h> header file. */
/* #undef HAVE_FPU_CONTROL_H */

/* Define if the compiler supports abi::__cxa_demangle(...) */
#define HAVE_GCC_ABI_DEMANGLE

/* Define if the C++ compiler knows how to compile __attribute__((constructor)) */
#define HAVE_TEUCHOS_CXX_ATTRIBUTE_CONSTRUCTOR

/* Define if the C++ compiler knows how to compile __attribute__((weak)), and
   if a program can test weak functions and call them if they are not NULL. */
/* #undef HAVE_TEUCHOS_CXX_ATTRIBUTE_WEAK */

/* Define if the C++ compiler knows how to compile "#pragma weak", and
   if a program can test weak functions and call them if they are not NULL. */
/* #undef HAVE_TEUCHOS_CXX_PRAGMA_WEAK */

/* Define if building dynamic shared libraries (instead of static libraries) */
/* #undef HAVE_TEUCHOS_DYNAMIC_LIBS */

/* Define if the (Windows) compiler has intrinsic datatype __int64 */
/* #undef HAVE_TEUCHOS___INT64 */

/* Not namespaced so should be deprecated. */
/* #undef HAVE_MPI */

/* #undef HAVE_TEUCHOS_MPI */

/* detected problems with the blas and solution methods */
#define HAVE_TEUCHOS_BLASFLOAT
/* #define HAVE_TEUCHOS_BLASFLOAT_APPLE_VECLIB_BUGFIX */
/* #undef HAVE_TEUCHOS_BLASFLOAT_DOUBLE_RETURN */

#define HAVE_SLAPY2_PROBLEM
#define HAVE_SLAPY2_DOUBLE_RETURN

/* #undef HAVE_COMPLEX_BLAS */
/* #undef HAVE_COMPLEX_BLAS_PROBLEM */
/* #undef HAVE_FIXABLE_COMPLEX_BLAS_PROBLEM */
/* #undef HAVE_VECLIB_COMPLEX_BLAS */

/* define if the compiler supports access of protected templated nested
   classes in derived classes */
/* #undef HAVE_PROTECTED_NESTED_TEMPLATE_CLASS_ACCESS */

/* #undef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK */

/* #undef HAVE_TEUCHOS_THREAD_SAFE */

/* #undef HAVE_TEUCHOS_LAPACKLARND */

/* Deprecated */
/* #undef HAVE_TEUCHOS_BOOST */

/* Deprecated */
/* #undef HAVE_TEUCHOS_QT */

/* #undef HAVE_TEUCHOS_QD */

/* #undef HAVE_TEUCHOSNUMERICS_EIGEN */

/* #undef HAVE_TEUCHOS_DOUBLE_TO_QD */

/* #undef HAVE_TEUCHOS_ARPREC */

/* #undef HAVE_TEUCHOS_DOUBLE_TO_ARPREC */

#define HAVE_TEUCHOS_ADD_TIME_MONITOR_TO_STACKED_TIMER

/* #undef HAVE_TEUCHOS_COMM_TIMERS */

/* #undef HAVE_TEUCHOS_FLOAT */

/* #undef HAVE_TEUCHOS_LONG_DOUBLE */

#define TEUCHOS_ORDINAL_TYPE ptrdiff_t

/* #undef HAVE_TEUCHOS_COMPLEX */

/* #undef HAVE_TEUCHOS_INST_FLOAT */

/* #undef HAVE_TEUCHOS_INST_COMPLEX_FLOAT */

/* #undef HAVE_TEUCHOS_INST_COMPLEX_DOUBLE */

// This exists only for backwards compatibility.
#define HAVE_TEUCHOS_LONG_LONG_INT 1

/* #undef HAVE_TEUCHOS_DEBUG */

/* #undef HAVE_TEUCHOS_DEBUG_RCP_NODE_TRACING */

#define HAVE_TEUCHOS_DEMANGLE

/* #undef HAVE_TEUCHOS_EXPAT */

#define HAVE_TEUCHOS_EXPLICIT_INSTANTIATION

#define HAVE_TEUCHOS_EXTENDED

/* #undef HAVE_TEUCHOS_GNU_MP */

/* #undef HAVE_TEUCHOS_LIBXML2 */

/* #undef HAVE_TEUCHOS_C_EXCEPTIONS */

/* #undef HAVE_TEUCHOS_LINK */

/* #undef HAVE_TEUCHOS_BFD */

/* #undef HAVE_TEUCHOS_STACKTRACE */

/* #undef HAVE_TEUCHOS_DEFAULT_STACKTRACE */

/* #undef HAVE_TEUCHOS_GLOBALLY_REDUCE_UNITTEST_RESULTS */

/* #undef HAVE_TEUCHOS_TIME_MASSIF_SNAPSHOTS */

/* #undef HAVE_TEUCHOS_KOKKOS_PROFILING */

/* template qualifier required for calling template methods from non-template
   code */
#define INVALID_TEMPLATE_QUALIFIER 

#ifndef TEUCHOS_DEPRECATED
#  if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#    define TEUCHOS_DEPRECATED  __attribute__((__deprecated__))
#  else
#    define TEUCHOS_DEPRECATED
#  endif
#endif

#ifndef TEUCHOS_DEPRECATED_MSG
#  if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 5))
#    define TEUCHOS_DEPRECATED_MSG(MSG)  __attribute__((__deprecated__ (#MSG) ))
#  elif (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#    define TEUCHOS_DEPRECATED_MSG(MSG)  __attribute__((__deprecated__))
#  else
#    define TEUCHOS_DEPRECATED_MSG(MSG)
#  endif
#endif

