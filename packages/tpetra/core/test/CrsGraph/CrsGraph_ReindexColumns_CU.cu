#include <Tpetra_ConfigDefs.hpp>
#ifdef HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT
#  include <KokkosCore_config.h>
#  ifdef KOKKOS_HAVE_CUDA
#    define KOKKOS_USE_CUDA_BUILD
#    include "CrsGraph/CrsGraph_ReindexColumns.cpp"
#    undef KOKKOS_USE_CUDA_BUILD
#  endif // KOKKOS_HAVE_CUDA
#endif // HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT
