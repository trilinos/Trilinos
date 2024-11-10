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

#define INSTANTIATE_UQ_PCE_STORAGE_SD(INSTMACRO, STORAGE, N)      \
  INSTMACRO( Sacado::UQ::PCE<STORAGE>, N )

#define INSTANTIATE_UQ_PCE_DS_SLD(INSTMACRO, S, L, D, LO, GO, N)       \
  typedef Stokhos::DynamicStorage<L,S,D::execution_space> DS_ ## L ## _ ## S ## _ ## _ ## D; \
  INSTANTIATE_UQ_PCE_STORAGE(INSTMACRO, DS_ ## L ## _ ## S ## _ ## _ ## D, LO, GO, N)

#define INSTANTIATE_UQ_PCE_DS_SLD_SD(INSTMACRO, S, L, D, N)       \
  typedef Stokhos::DynamicStorage<L,S,D::execution_space> DS_ ## L ## _ ## S ## _ ## _ ## D ## _2; \
  INSTANTIATE_UQ_PCE_STORAGE_SD(INSTMACRO, DS_ ## L ## _ ## S ## _ ## _ ## D ## _2, N)

#define INSTANTIATE_UQ_PCE_S_D(INSTMACRO, D, LO, GO, N) \
  INSTANTIATE_UQ_PCE_DS_SLD(INSTMACRO, double, int, D, LO, GO, N)

#define INSTANTIATE_UQ_PCE_S_D_SD(INSTMACRO, D, N) \
  INSTANTIATE_UQ_PCE_DS_SLD_SD(INSTMACRO, double, int, D, N)

#define INSTANTIATE_UQ_PCE_S(INSTMACRO, LO, GO, N) \
  typedef Stokhos::DeviceForNode2<N>::type DFN_ ## LO ## _ ## GO ## _ ## N; \
  INSTANTIATE_UQ_PCE_S_D(INSTMACRO, DFN_ ## LO ## _ ## GO ## _ ## N, LO, GO, N)

#if @IS_DEVICE_NODE@

// Add instantiation on HostMirror for device nodes.
#define INSTANTIATE_UQ_PCE_S_SD(INSTMACRO, N) \
  typedef Stokhos::DeviceForNode2<N>::type DFN_ ## N; \
  typedef Kokkos::View<double*, N::device_type>::HostMirror::device_type host_device_type_##N; \
  INSTANTIATE_UQ_PCE_S_D_SD(INSTMACRO, DFN_ ## N, N) \
  INSTANTIATE_UQ_PCE_S_D_SD(INSTMACRO, DFN_ ## N, host_device_type_##N)

#else

#define INSTANTIATE_UQ_PCE_S_SD(INSTMACRO, N) \
  typedef Stokhos::DeviceForNode2<N>::type DFN_ ## N; \
  INSTANTIATE_UQ_PCE_S_D_SD(INSTMACRO, DFN_ ## N, N)

#endif

#define INSTANTIATE_UQ_PCE(INSTMACRO, LO, GO, N) \
  INSTANTIATE_UQ_PCE_S(INSTMACRO, LO, GO, N)

#define INSTANTIATE_UQ_PCE_SD(INSTMACRO, N) \
  INSTANTIATE_UQ_PCE_S_SD(INSTMACRO, N)

#define INSTANTIATE_TPETRA_UQ_PCE_N(INSTMACRO, N)  \
  using default_local_ordinal_type = Tpetra::Map<>::local_ordinal_type; \
  using default_global_ordinal_type = Tpetra::Map<>::global_ordinal_type; \
  INSTANTIATE_UQ_PCE_S(INSTMACRO, default_local_ordinal_type, default_global_ordinal_type, N)

#define INSTANTIATE_TPETRA_UQ_PCE_N_SD(INSTMACRO, N)  \
  INSTANTIATE_UQ_PCE_S_SD(INSTMACRO, N)

#define INSTANTIATE_TPETRA_UQ_PCE_WRAPPER_NODES(INSTMACRO) \
  INSTMACRO(Tpetra_KokkosCompat_Kokkos@DEVICE@WrapperNode)

#define INSTANTIATE_TPETRA_UQ_PCE(INSTMACRO)                    \
  namespace Tpetra {                                            \
                                                                \
  TPETRA_ETI_MANGLING_TYPEDEFS()                                \
                                                                \
  INSTANTIATE_TPETRA_UQ_PCE_WRAPPER_NODES(INSTMACRO)            \
                                                                \
}
