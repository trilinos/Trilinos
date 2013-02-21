/* src/AztecOO_config.h.in.  Generated from configure.ac by autoheader.  */

/* Define the Fortran name mangling to be used for the BLAS */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */


/* Define to enable capture-matrix feature */
/* #undef AZ_ENABLE_CAPTURE_MATRIX */

/* Define to enable Teuchos TimeMonitors within Aztec solvers */
/* #undef AZ_ENABLE_TIMEMONITOR */

/* Define to 1 if you have the <assert.h> header file. */
/* #undef HAVE_ASSERT_H */

/* Define if you want to build AZ_lu */
#define HAVE_AZLU

/* Define if want to build with aztecoo enabled */
/* #undef HAVE_AZTECOO_EPETRAEXT */

/* Define if want to build aztecoo-examples */
/* #undef HAVE_AZTECOO_EXAMPLES */

/* Define if want to build with aztecoo enabled */
/* #undef HAVE_AZTECOO_IFPACK */

/* Define if want to build aztecoo-tests */
/* #undef HAVE_AZTECOO_TESTS */

/* Define if want to build with aztecoo enabled */
#define HAVE_AZTECOO_TEUCHOS

/* Define if you have a BLAS library. */
/* #undef HAVE_BLAS */

/* define if bool is a built-in type */
#define HAVE_BOOL

/* Define to 1 if you have the <cfloat> header file. */
/* #undef HAVE_CFLOAT */

/* Define to 1 if you have the <ctime> header file. */
/* #undef HAVE_CTIME */

/* Define if want to build examples */
/* #undef HAVE_EXAMPLES */

/* Define if you want to build export makefiles. */
/* #undef HAVE_EXPORT_MAKEFILES */

/* Define to 1 if you have the <float.h> header file. */
/* #undef HAVE_FLOAT_H */

/* Define to 1 if you have the `floor' function. */
/* #undef HAVE_FLOOR */

/* Define if want to build with fortran enabled */
#define HAVE_FORTRAN_SUPPORT

/* Define if you are using gnumake - this will shorten your link lines. */
/* #undef HAVE_GNUMAKE */

/* Define if you have LAPACK library. */
/* #undef HAVE_LAPACK */

/* Define if want to build libcheck */
/* #undef HAVE_LIBCHECK */

/* Define to 1 if your system has a GNU libc compatible `malloc' function, and
   to 0 otherwise. */
/* #undef HAVE_MALLOC */

/* Define to 1 if you have the <math.h> header file. */
/* #undef HAVE_MATH_H */

/* Define to 1 if you have the <memory.h> header file. */
/* #undef HAVE_MEMORY_H */

/* define if we want to use MPI */
#define HAVE_MPI

/* define if the compiler supports the mutable keyword */
#define HAVE_MUTABLE

/* define if the compiler implements namespaces */
#define HAVE_NAMESPACES

/* define if the compiler accepts the new for scoping rules */
#define HAVE_NEW_FOR_SCOPING

/* Define to 1 if you have the `pow' function. */
#define HAVE_POW

/* Define to 1 if you have the `sqrt' function. */
/* #undef HAVE_SQRT */

/* Define to 1 if you have the <stdio.h> header file. */
/* #undef HAVE_STDIO_H */

/* Define to 1 if you have the <stdlib.h> header file. */
/* #undef HAVE_STDLIB_H */

/* define if std::sprintf is supported */
/* #undef HAVE_STD_SPRINTF */

/* define if the compiler supports Standard Template Library */
#define HAVE_STL

/* Define to 1 if you have the `strchr' function. */
/* #undef HAVE_STRCHR */

/* Define to 1 if you have the <string> header file. */
#define HAVE_STRING

/* Define to 1 if you have the <strings.h> header file. */
/* #undef HAVE_STRINGS_H */

/* Define to 1 if you have the <string.h> header file. */
/* #undef HAVE_STRING_H */

/* Define to 1 if you have the <sys/time.h> header file. */
#ifndef HAVE_SYS_TIME_H
#define HAVE_SYS_TIME_H
#endif

/* Define if want to build tests */
/* #undef HAVE_TESTS */

/* Define to 1 if you have the `uname' function. */
/* #undef HAVE_UNAME */

/* Define to rpl_malloc if the replacement function should be used. */
/* #undef malloc */
