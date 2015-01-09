#include <Tpetra_ConfigDefs.hpp>
#ifdef HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT
#include <KokkosCore_config.h>
#ifdef KOKKOS_HAVE_CUDA
#define KOKKOS_USE_CUDA_BUILD
#include "inout/MatrixMarket_MV_Output_Perm.cpp"
#undef KOKKOS_USE_CUDA_BUILD
#endif
#endif
