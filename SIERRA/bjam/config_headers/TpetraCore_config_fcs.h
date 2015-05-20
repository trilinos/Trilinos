#ifndef TPETRACORE_CONFIG_H
#define TPETRACORE_CONFIG_H
/* CMake uses this file to generate TpetraCore_config.h automatically */

/* define if new form of std::count is supported */
/* #undef HAVE_STD_NEW_COUNT_SYNTAX */

/* Define if want to build with epetra enabled */
#define HAVE_TPETRACORE_EPETRA
#ifdef HAVE_TPETRACORE_EPETRA
#  define HAVE_TPETRA_EPETRA
#endif // HAVE_TPETRACORE_EPETRA

/* Define if want to build teuchos-debug */
/* #undef HAVE_TPETRA_DEBUG */

#define HAVE_TPETRA_ENABLE_SS_TESTING

/* Define if we have MPI */
#define HAVE_TPETRACORE_MPI
#ifdef HAVE_TPETRACORE_MPI
#  define HAVE_TPETRA_MPI
#endif // HAVE_TPETRACORE_MPI

/* Determine if we have CUDA */
/* #undef HAVE_TPETRACORE_CUDA */
#ifdef HAVE_TPETRACORE_CUDA
#  define HAVE_TPETRA_CUDA
#endif // HAVE_TPETRACORE_CUDA

/* Determine if we have the quadmath TPL */
/* #undef HAVE_TPETRACORE_QUADMATH */

/* Determine if we have QD */
/* #undef HAVE_TPETRACORE_QD */
#ifdef HAVE_TPETRACORE_QD
#  define HAVE_TPETRA_QD
#endif // HAVE_TPETRACORE_QD

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

/* Define when building Tpetra with the TpetraTSQR subpackage enabled.
   The HAVE_TPETRA_TSQR macro tells you whether it is safe to use TSQR
   in Tpetra (and downstream packages).  TSQR is enabled by default in
   Tpetra if KokkosTSQR is enabled, but users may turn off TSQR by
   setting Tpetra_ENABLE_TSQR to OFF. */
#define HAVE_TPETRACORE_TPETRATSQR
#ifdef HAVE_TPETRACORE_TPETRATSQR
#  define HAVE_TPETRA_KOKKOSTSQR
#  define HAVE_TPETRA_TPETRATSQR
#endif // HAVE_TPETRACORE_TPETRATSQR

/* Define when building TpetraCore with TSQR enabled.  TSQR is enabled
   by default in TpetraCore if the TpetraTSQR subpackage is enabled.
   Users may turn off TSQR even if TpetraTSQR is enabled, by setting
   Tpetra_ENABLE_TSQR to OFF. */
#define HAVE_TPETRA_TSQR

/* Define when enabling the Murmur hash function in Tpetra */
/* #undef TPETRA_USE_MURMUR_HASH */

/* Define when enabling KokkosContainers in Tpetra */
/* #undef HAVE_TPETRA_MMM_TIMINGS */

/* Define when enabling KokkosCore in Tpetra */
#define HAVE_TPETRACORE_KOKKOSCORE
#ifdef HAVE_TPETRACORE_KOKKOSCORE
// For backwards compatibility
#  define HAVE_TPETRA_KOKKOSCORE
#endif // HAVE_TPETRACORE_KOKKOSCORE

/* Define when enabling KokkosCompat in TpetraCore */
#define HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT
#ifdef HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT
// For backwards compatibility
#  define HAVE_TPETRACORE_KOKKOSCOMPAT
#  define HAVE_TPETRA_TEUCHOSKOKKOSCOMPAT
#  define HAVE_TPETRA_KOKKOSCOMPAT
#endif // HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT

/* Define when enabling TpetraKernels in TpetraCore */
/* #undef HAVE_TPETRACORE_TPETRAKERNELS */
#ifdef HAVE_TPETRACORE_TPETRAKERNELS
// For backwards compatibility
#  define HAVE_TPETRA_TPETRAKERNELS
#  define HAVE_TPETRA_KOKKOSLINALG
#endif // HAVE_TPETRACORE_TPETRAKERNELS

/* Define when enabling KokkosContainers in TpetraCore */
#define HAVE_TPETRACORE_KOKKOSCONTAINERS
#ifdef HAVE_TPETRACORE_KOKKOSCONTAINERS
// For backwards compatibility
#  define HAVE_TPETRA_KOKKOSCONTAINERS
#endif // HAVE_TPETRACORE_KOKKOSCONTAINERS

/* Define when enabling KokkosComm in TpetraCore */
#define HAVE_TPETRACORE_TEUCHOSKOKKOSCOMM
#ifdef HAVE_TPETRACORE_TEUCHOSKOKKOSCOMM
// For backwards compatibility
#  define HAVE_TPETRACORE_KOKKOSMPICOMM
#  define HAVE_TPETRA_KOKKOSMPICOMM
#endif // HAVE_TPETRACORE_TEUCHOSKOKKOSCOMM

/* Define when enabling Kokkos::View DistObject in Tpetra */
/* #undef TPETRA_ENABLE_KOKKOS_DISTOBJECT */

/* Define when enabling RDMA support for MPI communication for CUDA GPUs */
/* #undef TPETRA_ENABLE_MPI_CUDA_RDMA */

/* Define when enabling Tpetra for refactoring to new KokkosCore interface */
#define TPETRA_HAVE_KOKKOS_REFACTOR

/* Define if you want to use the new Kokkos refactor version of Map */
/* #undef TPETRA_USE_KOKKOS_REFACTOR_MAP */

/* Define if user requested explicit instantiation of classes into libtpetra */
/* #undef HAVE_TPETRA_EXPLICIT_INSTANTIATION */

/* Define if user requested explicit instantiation over ordinal pair <int,long> into libtpetra */
/* #undef HAVE_TPETRA_INST_INT_INT */
/* #undef HAVE_TPETRA_INST_INT_LONG */
/* #undef HAVE_TPETRA_INST_INT_LONG_LONG */
/* #undef HAVE_TPETRA_INST_INT_UNSIGNED */

/* #undef HAVE_TPETRA_INST_FLOAT */
/* #undef HAVE_TPETRA_INST_DOUBLE */
/* #undef HAVE_TPETRA_INST_COMPLEX_FLOAT */
/* #undef HAVE_TPETRA_INST_COMPLEX_DOUBLE */
/* #undef HAVE_TPETRA_INST_DD_REAL */
/* #undef HAVE_TPETRA_INST_QD_REAL */

/* #undef HAVE_TPETRA_INST_SERIALCLASSIC */
/* #undef HAVE_TPETRA_INST_SERIAL */
/* #undef HAVE_TPETRA_INST_PTHREAD */
/* #undef HAVE_TPETRA_INST_OPENMP */
/* #undef HAVE_TPETRA_INST_CUDA */

#define HAVE_TPETRA_INT_INT
#define HAVE_TPETRA_INT_LONG
/* #undef HAVE_TPETRA_INT_LONG_LONG */
#define HAVE_TPETRA_INT_UNSIGNED

#define HAVE_TPETRA_FLOAT
#define HAVE_TPETRA_DOUBLE
#define HAVE_TPETRA_COMPLEX_FLOAT
#define HAVE_TPETRA_COMPLEX_DOUBLE
/* #undef HAVE_TPETRA_DD_REAL */
/* #undef HAVE_TPETRA_QD_REAL */

#define HAVE_TPETRA_SERIALCLASSIC
#define HAVE_TPETRA_SERIAL
#define HAVE_TPETRA_PTHREAD
/* #undef HAVE_TPETRA_OPENMP */
/* #undef HAVE_TPETRA_CUDA */

/* #undef HAVE_TPETRA_THREADED_MKL */

/* #undef HAVE_TPETRACORE_RTI */
#ifdef HAVE_TPETRACORE_RTI
#  define HAVE_TPETRA_RTI
#endif // HAVE_TPETRACORE_RTI

#define TPETRA_DEPRECATED
#define TPETRA_DEPRECATED_MSG(MSG)


#endif // TPETRACORE_CONFIG_H
