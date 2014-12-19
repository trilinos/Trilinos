#include <Tpetra_ConfigDefs.hpp>
#ifdef HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT
#include <KokkosCore_config.h>
#ifdef KOKKOS_HAVE_CUDA
#define KOKKOS_USE_CUDA_BUILD
#include "inout/MatrixMarket_Tpetra_MultiVector_InOutTest.cpp"
#undef KOKKOS_USE_CUDA_BUILD
#endif
#endif
