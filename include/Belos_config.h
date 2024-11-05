/* Define if want to build with Amesos enabled */
/* #undef HAVE_BELOS_AMESOS */

/* Define if want to build with AztecOO enabled */
#define HAVE_BELOS_AZTECOO

/* Define if want to build with IFPACK enabled */
#define HAVE_BELOS_IFPACK

/* Define if want to build with Epetra enabled */
#define HAVE_BELOS_EPETRA

/* Define if want to build with EpetraExt enabled */
/* #undef HAVE_BELOS_EPETRAEXT */

/* Define if want to build with KokkosKernels enabled */
/* #undef HAVE_BELOS_KOKKOSKERNELS */

/* Define if want to build with ML enabled */
#define HAVE_BELOS_ML

/* Define if want to build with tpetra enabled */
/* #undef HAVE_BELOS_TPETRA */

/* Define if want to build with xpetra enabled */
/* #undef HAVE_BELOS_XPETRA */

/* Define whether we want to time Tpetra MultiVector operations */
/* #undef HAVE_BELOS_TPETRA_TIMERS */

/* Define if want to build with thyra enabled */
/* #undef HAVE_BELOS_THYRA */

/* Define if want to build with triutils enabled */
#define HAVE_BELOS_TRIUTILS

/* Define if want to build belos-examples */
/* #undef HAVE_BELOS_EXAMPLES */

/* Define if want to build examples */
/* #undef HAVE_EXAMPLES */

/* Define if want to build belos-tests */
/* #undef HAVE_BELOS_TESTS */

/* Define if want to build tests */
/* #undef HAVE_TESTS */

/* define if we want to use MPI */
#define HAVE_MPI

/* define if the compiler supports the mutable keyword */
#define HAVE_MUTABLE

/* define if the compiler implements namespaces */
#define HAVE_NAMESPACES

/* define if the compiler accepts the new for scoping rules */
#define HAVE_NEW_FOR_SCOPING

/* Define if want to build teuchos-complex */
/* #undef HAVE_TEUCHOS_COMPLEX */

/* Define if we are building Belos TSQR with TSQR support */
/* #undef HAVE_BELOS_TSQR */

/* Define if we are building with Teuchos TimeMonitors enabled */
#define BELOS_TEUCHOS_TIME_MONITOR

/* Define if we are building with experimental code enabled */
/* #undef HAVE_BELOS_EXPERIMENTAL */

/* Add macros for declaring functions deprecated */
#ifndef BELOS_DEPRECATED
#  if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#    define BELOS_DEPRECATED  __attribute__((__deprecated__))
#  else
#    define BELOS_DEPRECATED
#  endif
#endif

#ifndef BELOS_DEPRECATED_MSG
#  if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 5))
#    define BELOS_DEPRECATED_MSG(MSG)  __attribute__((__deprecated__ (#MSG) ))
#  elif (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#    define BELOS_DEPRECATED_MSG(MSG)  __attribute__((__deprecated__))
#  else
#    define BELOS_DEPRECATED_MSG(MSG)
#  endif
#endif


/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif
