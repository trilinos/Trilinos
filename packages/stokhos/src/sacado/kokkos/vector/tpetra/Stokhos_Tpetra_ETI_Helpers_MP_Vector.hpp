// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// MP::Vector includes
#include "Stokhos_Tpetra_MP_Vector.hpp"
#include "TpetraCore_ETIHelperMacros.h"

#define INSTANTIATE_MP_VECTOR_STORAGE(INSTMACRO, STORAGE, LO, GO, N)      \
  INSTMACRO( Sacado::MP::Vector<STORAGE>, LO, GO, N )

#define INSTANTIATE_MP_VECTOR_SFS_SLND(INSTMACRO, S, L, NUM, D, LO, GO, N) \
  typedef Stokhos::StaticFixedStorage<L,S,NUM,D::execution_space> SFS_ ## L ## _ ## S ## _ ## NUM ## _ ## D; \
  INSTANTIATE_MP_VECTOR_STORAGE(INSTMACRO, SFS_ ## L ## _ ## S ## _ ## NUM ## _ ## D, LO, GO, N)

#if defined(__MIC__)
// For MIC (Xeon Phi) -- vector width = 8 (double precision)
#define INSTANTIATE_MP_VECTOR_SFS_SLD_CPU(INSTMACRO, S, L, D, LO, GO, N) \
  INSTANTIATE_MP_VECTOR_SFS_SLND(INSTMACRO, S, L,  8, D, LO, GO, N)      \
  INSTANTIATE_MP_VECTOR_SFS_SLND(INSTMACRO, S, L, 16, D, LO, GO, N)      \
  INSTANTIATE_MP_VECTOR_SFS_SLND(INSTMACRO, S, L, 32, D, LO, GO, N)
#else
// For CPU with AVX instructions -- vector width = 4 (double precision)
#define INSTANTIATE_MP_VECTOR_SFS_SLD_CPU(INSTMACRO, S, L, D, LO, GO, N) \
  INSTANTIATE_MP_VECTOR_SFS_SLND(INSTMACRO, S, L,  4, D, LO, GO, N)      \
  INSTANTIATE_MP_VECTOR_SFS_SLND(INSTMACRO, S, L,  8, D, LO, GO, N)      \
  INSTANTIATE_MP_VECTOR_SFS_SLND(INSTMACRO, S, L, 16, D, LO, GO, N)      \
  INSTANTIATE_MP_VECTOR_SFS_SLND(INSTMACRO, S, L, 32, D, LO, GO, N)
#endif

// For CUDA GPU -- warp size = 32
#define INSTANTIATE_MP_VECTOR_SFS_SLD_GPU(INSTMACRO, S, L, D, LO, GO, N) \
  INSTANTIATE_MP_VECTOR_SFS_SLND(INSTMACRO, S, L, 16, D, LO, GO, N)      \
  INSTANTIATE_MP_VECTOR_SFS_SLND(INSTMACRO, S, L, 32, D, LO, GO, N)

#define INSTANTIATE_MP_VECTOR_DS_SLD(INSTMACRO, S, L, D, LO, GO, N)       \
  typedef Stokhos::DynamicStorage<L,S,D> DS_ ## L ## _ ## S ## _ ## _ ## D; \
  INSTANTIATE_MP_VECTOR_STORAGE(INSTMACRO, DS_ ## L ## _ ## S ## _ ## _ ## D, LO, GO, N)

#define INSTANTIATE_MP_VECTOR_S_D_CPU(INSTMACRO, D, LO, GO, N) \
  INSTANTIATE_MP_VECTOR_SFS_SLD_CPU(INSTMACRO, double, int, D, LO, GO, N)
#define INSTANTIATE_MP_VECTOR_S_D_GPU(INSTMACRO, D, LO, GO, N) \
  INSTANTIATE_MP_VECTOR_SFS_SLD_GPU(INSTMACRO, double, int, D, LO, GO, N)

// Disabling dynamic storage ETI -- we don't really need it
//  INSTANTIATE_MP_VECTOR_DS_SLD(INSTMACRO, double, int, D, LO, GO, N)

#define INSTANTIATE_MP_VECTOR_S_CPU(INSTMACRO, LO, GO, N) \
  typedef Stokhos::DeviceForNode<N>::type DFN_CPU_ ## LO ## _ ## GO ## _ ## N; \
  INSTANTIATE_MP_VECTOR_S_D_CPU(INSTMACRO, DFN_CPU_ ## LO ## _ ## GO ## _ ## N, LO, GO, N)
#define INSTANTIATE_MP_VECTOR_S_GPU(INSTMACRO, LO, GO, N) \
  typedef Stokhos::DeviceForNode<N>::type DFN_GPU_ ## LO ## _ ## GO ## _ ## N; \
  INSTANTIATE_MP_VECTOR_S_D_GPU(INSTMACRO, DFN_GPU_ ## LO ## _ ## GO ## _ ## N, LO, GO, N)

#if defined(HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT) && defined(HAVE_TPETRA_INST_SERIAL)
#define INSTANTIATE_TPETRA_MP_VECTOR_SERIAL(INSTMACRO) \
  using default_local_ordinal_type = Tpetra::Map<>::local_ordinal_type; \
  using default_global_ordinal_type = Tpetra::Map<>::global_ordinal_type; \
  INSTANTIATE_MP_VECTOR_S_CPU(INSTMACRO, default_local_ordinal_type, default_global_ordinal_type, Tpetra_KokkosCompat_KokkosSerialWrapperNode)
#else
#define INSTANTIATE_TPETRA_MP_VECTOR_SERIAL(INSTMACRO)
#endif


#if defined(HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT) && defined(HAVE_TPETRA_INST_PTHREAD)
#define INSTANTIATE_TPETRA_MP_VECTOR_THREADS(INSTMACRO) \
  using default_local_ordinal_type = Tpetra::Map<>::local_ordinal_type; \
  using default_global_ordinal_type = Tpetra::Map<>::global_ordinal_type; \
  INSTANTIATE_MP_VECTOR_S_CPU(INSTMACRO, default_local_ordinal_type, default_global_ordinal_type, Tpetra_KokkosCompat_KokkosThreadsWrapperNode)
#else
#define INSTANTIATE_TPETRA_MP_VECTOR_THREADS(INSTMACRO)
#endif

#if defined(HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT) && defined(HAVE_TPETRA_INST_OPENMP)
#define INSTANTIATE_TPETRA_MP_VECTOR_OPENMP(INSTMACRO) \
  using default_local_ordinal_type = Tpetra::Map<>::local_ordinal_type; \
  using default_global_ordinal_type = Tpetra::Map<>::global_ordinal_type; \
  INSTANTIATE_MP_VECTOR_S_CPU(INSTMACRO, default_local_ordinal_type, default_global_ordinal_type, Tpetra_KokkosCompat_KokkosOpenMPWrapperNode)
#else
#define INSTANTIATE_TPETRA_MP_VECTOR_OPENMP(INSTMACRO)
#endif

#if defined(HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT) && defined(HAVE_TPETRA_INST_CUDA)
#define INSTANTIATE_TPETRA_MP_VECTOR_CUDA(INSTMACRO) \
  using default_local_ordinal_type = Tpetra::Map<>::local_ordinal_type; \
  using default_global_ordinal_type = Tpetra::Map<>::global_ordinal_type; \
  INSTANTIATE_MP_VECTOR_S_GPU(INSTMACRO, default_local_ordinal_type, default_global_ordinal_type, Tpetra_KokkosCompat_KokkosCudaWrapperNode)
#else
#define INSTANTIATE_TPETRA_MP_VECTOR_CUDA(INSTMACRO)
#endif

#define INSTANTIATE_TPETRA_MP_VECTOR_WRAPPER_NODES(INSTMACRO) \
  INSTANTIATE_TPETRA_MP_VECTOR_THREADS(INSTMACRO)             \
  INSTANTIATE_TPETRA_MP_VECTOR_OPENMP(INSTMACRO)              \
  INSTANTIATE_TPETRA_MP_VECTOR_CUDA(INSTMACRO)

#define INSTANTIATE_TPETRA_MP_VECTOR(INSTMACRO)                 \
  namespace Tpetra {                                            \
                                                                \
  TPETRA_ETI_MANGLING_TYPEDEFS()                                \
                                                                \
  INSTANTIATE_TPETRA_MP_VECTOR_WRAPPER_NODES(INSTMACRO)         \
                                                                \
}
