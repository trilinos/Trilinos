/* src/KokkosClassic_config.h.in.  Used by CMake to generate KokkosClassic_config.h. */

/* Define if you want to build kokkos-examples */
/* #undef HAVE_KOKKOSCLASSIC_EXAMPLES */

/* Define if you want to build kokkos-debug */
/* #undef HAVE_KOKKOSCLASSIC_DEBUG */

/* Define if you want to build Kokkos with OpenMP */
/* #undef HAVE_KOKKOSCLASSIC_OPENMP */

/* Define if the KokkosCompat subpackage is available in KokkosClassic */
/* #undef HAVE_KOKKOSCLASSIC_KOKKOSCOMPAT */

/* Define if you want to build kokkos-experimental */
/* #undef HAVE_KOKKOSCLASSIC_EXPERIMENTAL */

/* Enable float type for CUDA nodes. */
#define HAVE_KOKKOSCLASSIC_CUDA_FLOAT

/* Enable double type for CUDA nodes. */
/* #undef HAVE_KOKKOSCLASSIC_CUDA_DOUBLE */

/* #undef HAVE_KOKKOSCLASSIC_CUDA_NODE_MEMORY_PROFILING */
/* #undef HAVE_KOKKOSCLASSIC_CUDA_NODE_MEMORY_TRACE */
/* #undef HAVE_KOKKOSCLASSIC_TREAT_SERIALNODE_AS_DEVICE */

/* Enable complex<float> type for CUDA nodes. */
/* #undef HAVE_KOKKOSCLASSIC_CUDA_COMPLEX_FLOAT */

/* Enable complex<double> type for CUDA nodes. */
/* #undef HAVE_KOKKOSCLASSIC_CUDA_COMPLEX_DOUBLE */

/* Define if you want to build kokkos-tests */
/* #undef HAVE_KOKKOSCLASSIC_TESTS */

/* Other package dependencies */
/* #undef HAVE_KOKKOSCLASSIC_TEUCHOS */
#define HAVE_KOKKOSCLASSIC_THREADPOOL
/* #undef HAVE_KOKKOSCLASSIC_KOKKOSCORE */

/* TPL Dependencies */
#define HAVE_KOKKOSCLASSIC_SERIAL
/* #undef HAVE_KOKKOSCLASSIC_TBB */
/* #undef HAVE_KOKKOSCLASSIC_CUDA */
/* #undef HAVE_KOKKOSCLASSIC_CUSPARSE */
/* #undef HAVE_KOKKOSCLASSIC_THRUST */
/* #undef HAVE_KOKKOSCLASSIC_CUSP */
/* #undef HAVE_KOKKOSCLASSIC_MKL */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif
