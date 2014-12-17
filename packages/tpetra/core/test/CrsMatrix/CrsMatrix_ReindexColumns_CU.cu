// Ensure that if CUDA and KokkosCompat are enabled, then
// only the .cu version of this file is actually compiled.
#include <Tpetra_config.h>
#ifdef HAVE_TPETRA_KOKKOSCOMPAT
#  include <KokkosCore_config.h>
#  ifdef KOKKOS_HAVE_CUDA
#    define KOKKOS_USE_CUDA_BUILD
#    include "CrsMatrix/CrsMatrix_ReindexColumns.cpp"
#    undef KOKKOS_USE_CUDA_BUILD
#  endif // KOKKOS_HAVE_CUDA
#endif // HAVE_TPETRA_KOKKOSCOMPAT
