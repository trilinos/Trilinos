/* src/Tpetra_config.h.in.  Generated from configure.ac by autoheader.  */

/* define if new form of std::count is supported */
/* #undef HAVE_STD_NEW_COUNT_SYNTAX */

/* Define if want to build with epetra enabled */
#define HAVE_TPETRA_EPETRA

/* Define if want to build teuchos-debug */
/* #undef HAVE_TPETRA_DEBUG */

#define HAVE_TPETRA_ENABLE_SS_TESTING

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

/* Define when building Tpetra with the KokkosTSQR subpackage enabled.
   The HAVE_TPETRA_TSQR macro tells you whether it is safe to use TSQR
   in Tpetra (and downstream packages).  TSQR is enabled by default in
   Tpetra if KokkosTSQR is enabled, but users may turn off TSQR by
   setting Tpetra_ENABLE_TSQR to OFF. */
/* #undef HAVE_TPETRA_KOKKOSTSQR */

/* Define when building Tpetra with TSQR enabled.  TSQR is enabled by
   default in Tpetra if the KokkosTSQR subpackage is enabled.  Users
   may turn off TSQR even if KokkosTSQR is enabled, by setting
   Tpetra_ENABLE_TSQR to OFF. */
/* #undef HAVE_TPETRA_TSQR */

/* Define when the variable-block-size classes VbrMatrix, BlockMap,
   BlockCrsGraph, and BlockMultiVector (in the Tpetra namespace) are
   enabled. */
#define HAVE_TPETRA_CLASSIC_VBR

/* Define when enabling the Murmur hash function in Tpetra */
/* #undef TPETRA_USE_MURMUR_HASH */

/* Define when enabling KokkosContainers in Tpetra */
/* #undef HAVE_TPETRA_MMM_TIMINGS */

/* Define when enabling KokkosCore in Tpetra */
#define HAVE_TPETRA_KOKKOSCORE

/* Define when enabling KokkosCompat in Tpetra */
/* #undef HAVE_TPETRA_KOKKOSCOMPAT */

/* Define when enabling KokkosLinAlg in Tpetra */
/* #undef HAVE_TPETRA_KOKKOSLINALG */

/* Define when enabling KokkosContainers in Tpetra */
#define HAVE_TPETRA_KOKKOSCONTAINERS

/* Define when enabling KokkosContainers in Tpetra */
/* #undef HAVE_TPETRA_KOKKOSMPICOMM */

/* Define when enabling Kokkos::View DistObject in Tpetra */
/* #undef TPETRA_ENABLE_KOKKOS_DISTOBJECT */

/* Define when enabling RDMA support for MPI communication for CUDA GPUs */
/* #undef TPETRA_ENABLE_MPI_CUDA_RDMA */

/* Define when enabling Tpetra for refactoring to new KokkosCore interface */
/* #undef TPETRA_HAVE_KOKKOS_REFACTOR */

/* Define if you want to use the new Kokkos refactor version of Map */
/* #undef TPETRA_USE_KOKKOS_REFACTOR_MAP */

/* Define if user requested explicit instantiation of classes into libtpetra */
/* #undef HAVE_TPETRA_EXPLICIT_INSTANTIATION */

/* Define if user requested explicit instantiation over ordinal pair <int,long> into libtpetra */
#define HAVE_TPETRA_INST_INT_LONG
/* #undef HAVE_TPETRA_INST_INT_LONG_LONG */

#define HAVE_TPETRA_INST_FLOAT
#define HAVE_TPETRA_INST_DOUBLE
#define HAVE_TPETRA_INST_COMPLEX_FLOAT
#define HAVE_TPETRA_INST_COMPLEX_DOUBLE
/* #undef HAVE_TPETRA_INST_DD_REAL */
/* #undef HAVE_TPETRA_INST_QD_REAL */

/* #undef HAVE_TPETRA_THREADED_MKL */

/* #undef HAVE_TPETRA_RTI */

#define TPETRA_DEPRECATED
