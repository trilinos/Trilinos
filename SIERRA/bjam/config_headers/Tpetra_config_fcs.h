/* src/Tpetra_config.h.in.  Generated from configure.ac by autoheader.  */

/* define if new form of std::count is supported */
/* #undef HAVE_STD_NEW_COUNT_SYNTAX */

/* Define if want to build with epetra enabled */
#define HAVE_TPETRA_EPETRA

/* Define if want to build teuchos-debug */
/* #undef HAVE_TPETRA_DEBUG */

/* #undef HAVE_TPETRA_ENABLE_SS_TESTING */

/* Define if want to build tpetra with OpenMP */
/* #undef HAVE_TPETRA_OPENMP */

/* Define if we have MPI */
#define HAVE_TPETRA_MPI

/* Determine if we have CUDA, Thrust, QD */
/* #undef HAVE_TPETRA_CUDA */
/* #undef HAVE_TPETRA_THRUST */
/* #undef HAVE_TPETRA_QD */

/* Define if want to build tpetra-throw_warnings */
/* #undef HAVE_TPETRA_THROW_WARNINGS */

/* Define if want to build tpetra-print_warnings */
/* #undef HAVE_TPETRA_PRINT_WARNINGS */

/* Define if want to build tpetra-throw_efficiency_warnings */
/* #undef HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS */

/* Define if want to build tpetra-print_efficiency_warnings */
/* #undef HAVE_TPETRA_PRINT_EFFICIENCY_WARNINGS */

/* Define if want to build tpetra-throw_abuse_warnings */
/* #undef HAVE_TPETRA_THROW_ABUSE_WARNINGS */

/* Define if want to build tpetra-print_abuse_warnings */
/* #undef HAVE_TPETRA_PRINT_ABUSE_WARNINGS */

/* Define when building Tpetra with TSQR support */
#define HAVE_TPETRA_TSQR

/* Define when enabling the Murmur hash function in Tpetra */
/* #undef TPETRA_USE_MURMUR_HASH */

/* Define if user requested explicit instantiation of classes into libtpetra */
/* #undef HAVE_TPETRA_EXPLICIT_INSTANTIATION */

/* Define if user requested explicit instantiation over ordinal pair <int,long> into libtpetra */
#define HAVE_TPETRA_INST_INT_LONG

/* #undef HAVE_TPETRA_INST_FLOAT */
#define HAVE_TPETRA_INST_DOUBLE
/* #undef HAVE_TPETRA_INST_COMPLEX_FLOAT */
/* #undef HAVE_TPETRA_INST_COMPLEX_DOUBLE */
/* #undef HAVE_TPETRA_INST_DD_REAL */
/* #undef HAVE_TPETRA_INST_QD_REAL */

/* #undef HAVE_TPETRA_THREADED_MKL */

/* #undef HAVE_TPETRA_RTI */

#ifndef TPETRA_DEPRECATED
#  if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#    define TPETRA_DEPRECATED  __attribute__((__deprecated__))
#  else
#    define TPETRA_DEPRECATED
#  endif
#endif

