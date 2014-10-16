#ifndef KOKKOS_CORE_CONFIG_H
#define KOKKOS_CORE_CONFIG_H

/* The trivial 'src/build_common.sh' creates a config
 * that must stay in sync with this file.
 */
#define KOKKOS_FOR_SIERRA

#if !defined( KOKKOS_FOR_SIERRA )

#define KOKKOS_HAVE_MPI
/* #undef KOKKOS_HAVE_CUDA */
/* #undef KOKKOS_USE_CUDA_UVM */
#define KOKKOS_HAVE_PTHREAD
/* #undef KOKKOS_HAVE_QTHREAD */
/* #undef KOKKOS_HAVE_Winthread */
#define KOKKOS_HAVE_OPENMP
/* #undef KOKKOS_HAVE_HWLOC */
/* #undef KOKKOS_EXPRESSION_CHECK */
/* #undef KOKKOS_HAVE_CXX11 */

#endif

#endif
