// Ensure that if CUDA and KokkosCompat are enabled, then only the .cu
// version of this file will actually be compiled.
#include <Tpetra_config.h>
#ifdef HAVE_TPETRA_KOKKOSCOMPAT
#  include <KokkosCore_config.h>
#  ifdef KOKKOS_HAVE_CUDA
#    define KOKKOS_USE_CUDA_BUILD
#  include "MultiVector/rcb.cpp"
#  undef KOKKOS_USE_CUDA_BUILD
#  endif
#endif

