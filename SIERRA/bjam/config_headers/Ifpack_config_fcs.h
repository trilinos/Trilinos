/* src/Ifpack_config.h.in.  Generated from configure.ac by autoheader.  */

/* Define if want to build with ifpack enabled */
#define HAVE_IFPACK_AMESOS

/* Define if want to build with ifpack enabled */
#define HAVE_IFPACK_AZTECOO

/* Define if want to build with ifpack enabled */
#define HAVE_IFPACK_EPETRAEXT

/* Define if want to build with ifpack enabled */
#define HAVE_IFPACK_GALERI

/* Define if you want to build Ifpack's ParameterList interface */
/* #undef HAVE_IFPACK_METIS */

/* define if we want to use MPI */
#define HAVE_MPI

/* Define if want to build with hypre enabled */
/* #undef HAVE_HYPRE */

/* Define if we want to use the euclid preconditioner */
/* #undef HAVE_EUCLID */

/* Define if we want to use the HIPS preconditioner */
/* #undef HAVE_IFPACK_HIPS */

/* Define if we want to build with SuperLU enabled */
/* #undef HAVE_IFPACK_SUPERLU */

/* Define if we want to build with SPARSKIT enabled */
/* #undef HAVE_IFPACK_SPARSKIT */

/* Define if we want Teuchos Time Monitors enabled */
#define IFPACK_TEUCHOS_TIME_MONITOR

/* Define if we want Ifpack Internal FlopCounters enabled */
/* #undef IFPACK_FLOPCOUNTERS */

/* Add macros for declaring functions deprecated */
#ifndef IFPACK_DEPRECATED
#  if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#    define IFPACK_DEPRECATED  __attribute__((__deprecated__))
#  else
#    define IFPACK_DEPRECATED
#  endif
#endif


