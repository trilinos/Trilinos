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

#define INSTANTIATE_MP_VECTOR_STORAGE_SD(INSTMACRO, STORAGE, N)      \
  INSTMACRO( Sacado::MP::Vector<STORAGE>, N )

#define INSTANTIATE_MP_VECTOR_SFS_SLND(INSTMACRO, S, L, NUM, D, LO, GO, N) \
  typedef Stokhos::StaticFixedStorage<L,S,NUM,D::execution_space> SFS_ ## L ## _ ## S ## _ ## NUM ## _ ## D; \
  INSTANTIATE_MP_VECTOR_STORAGE(INSTMACRO, SFS_ ## L ## _ ## S ## _ ## NUM ## _ ## D, LO, GO, N)

#define INSTANTIATE_MP_VECTOR_SFS_SLND_SD(INSTMACRO, S, L, NUM, D, N) \
  typedef Stokhos::StaticFixedStorage<L,S,NUM,D::execution_space> SFS_ ## L ## _ ## S ## _ ## NUM ## _ ## D ## _2; \
  INSTANTIATE_MP_VECTOR_STORAGE_SD(INSTMACRO, SFS_ ## L ## _ ## S ## _ ## NUM ## _ ## D ## _2, N)

#define INSTANTIATE_MP_VECTOR_SFS_SLD(INSTMACRO, S, L, D, LO, GO, N) \
  INSTANTIATE_MP_VECTOR_SFS_SLND(INSTMACRO, S, L,  @ENSEMBLE_SIZE@, D, LO, GO, N)

#define INSTANTIATE_MP_VECTOR_SFS_SLD_SD(INSTMACRO, S, L, D, N) \
  INSTANTIATE_MP_VECTOR_SFS_SLND_SD(INSTMACRO, S, L,  @ENSEMBLE_SIZE@, D, N)

#define INSTANTIATE_MP_VECTOR_S_D(INSTMACRO, D, LO, GO, N) \
  INSTANTIATE_MP_VECTOR_SFS_SLD(INSTMACRO, double, int, D, LO, GO, N)

#define INSTANTIATE_MP_VECTOR_S_D_SD(INSTMACRO, D, N) \
  INSTANTIATE_MP_VECTOR_SFS_SLD_SD(INSTMACRO, double, int, D, N)

#define INSTANTIATE_MP_VECTOR_S(INSTMACRO, LO, GO, N) \
  typedef Stokhos::DeviceForNode<N>::type DFN_CPU_ ## LO ## _ ## GO ## _ ## N; \
  INSTANTIATE_MP_VECTOR_S_D(INSTMACRO, DFN_CPU_ ## LO ## _ ## GO ## _ ## N, LO, GO, N)

#if @IS_DEVICE_NODE@

// Add instantiation on HostMirror for device nodes.
#define INSTANTIATE_MP_VECTOR_S_SD(INSTMACRO, N) \
  typedef Stokhos::DeviceForNode<N>::type DFN_CPU_ ## N; \
  typedef Kokkos::View<double*, N::device_type>::HostMirror::device_type host_device_type_##N; \
  INSTANTIATE_MP_VECTOR_S_D_SD(INSTMACRO, DFN_CPU_ ## N, N) \
  INSTANTIATE_MP_VECTOR_S_D_SD(INSTMACRO, DFN_CPU_ ## N, host_device_type_##N)

#else

#define INSTANTIATE_MP_VECTOR_S_SD(INSTMACRO, N) \
  typedef Stokhos::DeviceForNode<N>::type DFN_CPU_ ## N; \
  INSTANTIATE_MP_VECTOR_S_D_SD(INSTMACRO, DFN_CPU_ ## N, N)

#endif

#define INSTANTIATE_TPETRA_MP_VECTOR_WRAPPER_NODES(INSTMACRO)           \
  using default_local_ordinal_type = Tpetra::Map<>::local_ordinal_type; \
  using default_global_ordinal_type = Tpetra::Map<>::global_ordinal_type; \
  INSTANTIATE_MP_VECTOR_S(INSTMACRO, default_local_ordinal_type, default_global_ordinal_type, Tpetra_KokkosCompat_Kokkos@DEVICE@WrapperNode)

#define INSTANTIATE_TPETRA_MP_VECTOR_WRAPPER_NODES_SD(INSTMACRO)           \
  INSTANTIATE_MP_VECTOR_S_SD(INSTMACRO, Tpetra_KokkosCompat_Kokkos@DEVICE@WrapperNode)

#define INSTANTIATE_TPETRA_MP_VECTOR(INSTMACRO)                 \
  namespace Tpetra {                                            \
                                                                \
  TPETRA_ETI_MANGLING_TYPEDEFS()                                \
                                                                \
  INSTANTIATE_TPETRA_MP_VECTOR_WRAPPER_NODES(INSTMACRO)         \
                                                                \
}

#define INSTANTIATE_TPETRA_MP_VECTOR_SD(INSTMACRO)              \
  namespace Tpetra {                                            \
                                                                \
  TPETRA_ETI_MANGLING_TYPEDEFS()                                \
                                                                \
  INSTANTIATE_TPETRA_MP_VECTOR_WRAPPER_NODES_SD(INSTMACRO)      \
                                                                \
}
