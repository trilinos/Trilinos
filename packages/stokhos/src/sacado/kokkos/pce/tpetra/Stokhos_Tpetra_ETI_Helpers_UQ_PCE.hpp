// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// UQ::PCE includes
#include "Stokhos_Tpetra_UQ_PCE.hpp"
#include "TpetraCore_ETIHelperMacros.h"

#define INSTANTIATE_UQ_PCE_STORAGE(INSTMACRO, STORAGE, LO, GO, N)      \
  INSTMACRO( Sacado::UQ::PCE<STORAGE>, LO, GO, N )

#define INSTANTIATE_UQ_PCE_DS_SLD(INSTMACRO, S, L, D, LO, GO, N)       \
  typedef Stokhos::DynamicStorage<L,S,D::execution_space> DS_ ## L ## _ ## S ## _ ## _ ## D; \
  INSTANTIATE_UQ_PCE_STORAGE(INSTMACRO, DS_ ## L ## _ ## S ## _ ## _ ## D, LO, GO, N)

#define INSTANTIATE_UQ_PCE_S_D(INSTMACRO, D, LO, GO, N) \
  INSTANTIATE_UQ_PCE_DS_SLD(INSTMACRO, double, int, D, LO, GO, N)

#define INSTANTIATE_UQ_PCE_S(INSTMACRO, LO, GO, N) \
  typedef Stokhos::DeviceForNode2<N>::type DFN_ ## LO ## _ ## GO ## _ ## N; \
  INSTANTIATE_UQ_PCE_S_D(INSTMACRO, DFN_ ## LO ## _ ## GO ## _ ## N, LO, GO, N)

#define INSTANTIATE_UQ_PCE(INSTMACRO, LO, GO, N) \
  INSTANTIATE_UQ_PCE_S(INSTMACRO, LO, GO, N)

#define INSTANTIATE_TPETRA_UQ_PCE_N(INSTMACRO, N)  \
  using default_local_ordinal_type = Tpetra::Map<>::local_ordinal_type; \
  using default_global_ordinal_type = Tpetra::Map<>::global_ordinal_type; \
  INSTANTIATE_UQ_PCE_S(INSTMACRO, default_local_ordinal_type, default_global_ordinal_type, N)

#if defined(HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT) && defined(HAVE_TPETRA_INST_SERIAL)
#define INSTANTIATE_TPETRA_UQ_PCE_SERIAL(INSTMACRO) \
  INSTMACRO(Tpetra_KokkosCompat_KokkosSerialWrapperNode)
#else
#define INSTANTIATE_TPETRA_UQ_PCE_SERIAL(INSTMACRO)
#endif

#if defined(HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT) && defined(HAVE_TPETRA_INST_PTHREAD)
#define INSTANTIATE_TPETRA_UQ_PCE_THREADS(INSTMACRO) \
  INSTMACRO(Tpetra_KokkosCompat_KokkosThreadsWrapperNode)
#else
#define INSTANTIATE_TPETRA_UQ_PCE_THREADS(INSTMACRO)
#endif

#if defined(HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT) && defined(HAVE_TPETRA_INST_OPENMP)
#define INSTANTIATE_TPETRA_UQ_PCE_OPENMP(INSTMACRO) \
  INSTMACRO(Tpetra_KokkosCompat_KokkosOpenMPWrapperNode)
#else
#define INSTANTIATE_TPETRA_UQ_PCE_OPENMP(INSTMACRO)
#endif

#if defined(HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT) && defined(HAVE_TPETRA_INST_CUDA)
#define INSTANTIATE_TPETRA_UQ_PCE_CUDA(INSTMACRO) \
  INSTMACRO(Tpetra_KokkosCompat_KokkosCudaWrapperNode)
#else
#define INSTANTIATE_TPETRA_UQ_PCE_CUDA(INSTMACRO)
#endif

#define INSTANTIATE_TPETRA_UQ_PCE_WRAPPER_NODES(INSTMACRO) \
  INSTANTIATE_TPETRA_UQ_PCE_THREADS(INSTMACRO)             \
  INSTANTIATE_TPETRA_UQ_PCE_OPENMP(INSTMACRO)              \
  INSTANTIATE_TPETRA_UQ_PCE_CUDA(INSTMACRO)

#define INSTANTIATE_TPETRA_UQ_PCE(INSTMACRO)                    \
  namespace Tpetra {                                            \
                                                                \
  TPETRA_ETI_MANGLING_TYPEDEFS()                                \
                                                                \
  INSTANTIATE_TPETRA_UQ_PCE_WRAPPER_NODES(INSTMACRO)            \
                                                                \
}
