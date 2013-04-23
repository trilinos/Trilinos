
/* #undef BUILD_SHARED_LIBS */

#define HAVE_TEUCHOS_PARAMETERLIST

/* Define the Fortran name mangling to be used for the BLAS */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */


/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to 1 if you have the <fpu_control.h> header file. */
/* #undef HAVE_FPU_CONTROL_H */

/* define if the compiler supports abi::__cxa_demangle(...) */
#define HAVE_GCC_ABI_DEMANGLE

/* Define if the (Windows) compiler has intrinsic datatype __int64 */
/* #undef HAVE_TEUCHOS___INT64 */

/* Not namespaced so should be depreciated. */
#define HAVE_MPI

#define HAVE_TEUCHOS_MPI

/* detected problems with the blas and solution methods */
#define HAVE_TEUCHOS_BLASFLOAT
/* #undef HAVE_TEUCHOS_BLASFLOAT_APPLE_VECLIB_BUGFIX */
#define HAVE_TEUCHOS_BLASFLOAT_DOUBLE_RETURN

#define HAVE_SLAPY2_PROBLEM
#define HAVE_SLAPY2_DOUBLE_RETURN

#define HAVE_COMPLEX_BLAS
#define HAVE_COMPLEX_BLAS_PROBLEM
#define HAVE_FIXABLE_COMPLEX_BLAS_PROBLEM
/* #undef HAVE_VECLIB_COMPLEX_BLAS */

/* define if the compiler supports access of protected templated nested
   classes in derived classes */
/* #undef HAVE_PROTECTED_NESTED_TEMPLATE_CLASS_ACCESS */

/* #undef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK */

/* #undef HAVE_TEUCHOS_LAPACKLARND */

/* #undef HAVE_TEUCHOSCORE_BOOST */

/* Deprecated */
/* #undef HAVE_TEUCHOS_BOOST */

/* #undef HAVE_TEUCHOSCORE_QT */

/* Deprecated */
/* #undef HAVE_TEUCHOS_QT */

/* #undef HAVE_TEUCHOS_QD */

/* #undef HAVE_TEUCHOSNUMERICS_EIGEN */

/* Deprecated */
/* #undef HAVE_TEUCHOSCORE_QD */

/* #undef HAVE_TEUCHOS_DOUBLE_TO_QD */

/* #undef HAVE_TEUCHOS_ARPREC */

/* Deprecated */
/* #undef HAVE_TEUCHOSCORE_ARPREC */

/* #undef HAVE_TEUCHOS_DOUBLE_TO_ARPREC */

/* #undef HAVE_TEUCHOS_COMM_TIMERS */

/* #undef HAVE_TEUCHOS_FLOAT */

#define TEUCHOS_ORDINAL_TYPE ptrdiff_t

#define HAVE_TEUCHOS_COMPLEX

/* #undef HAVE_TEUCHOS_LONG_LONG_INT */

/* #undef HAVE_TEUCHOS_DEBUG */

/* #undef HAVE_TEUCHOS_DEBUG_RCP_NODE_TRACING */

/* #undef HAS_TEUCHOS_BOOST_IS_POLYMORPHIC */

#define HAVE_TEUCHOS_DEMANGLE

/* #undef HAVE_TEUCHOS_EXPAT */

/* #undef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION */

#define HAVE_TEUCHOS_EXTENDED

/* #undef HAVE_TEUCHOS_GNU_MP */

/* #undef HAVE_TEUCHOS_LIBXML2 */

/* #undef HAVE_TEUCHOS_C_EXCEPTIONS */

/* #undef HAVE_TEUCHOS_LINK */

/* #undef HAVE_TEUCHOS_BFD */

/* #undef HAVE_TEUCHOS_STACKTRACE */

/* #undef HAVE_TEUCHOS_DEFAULT_STACKTRACE */

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

