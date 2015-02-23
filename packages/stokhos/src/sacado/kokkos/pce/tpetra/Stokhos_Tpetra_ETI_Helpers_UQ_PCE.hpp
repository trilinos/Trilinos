// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

// UQ::PCE includes
#include "Stokhos_Tpetra_UQ_PCE.hpp"
#include "TpetraCore_ETIHelperMacros.h"

#define INSTANTIATE_UQ_PCE_STORAGE(INSTMACRO, STORAGE, LO, GO, N)      \
  INSTMACRO( Sacado::UQ::PCE<STORAGE>, LO, GO, N )

#define INSTANTIATE_UQ_PCE_DS_SLD(INSTMACRO, S, L, D, LO, GO, N)       \
  typedef Stokhos::DynamicStorage<L,S,typename D::execution_space> DS_ ## L ## _ ## S ## _ ## _ ## D; \
  INSTANTIATE_UQ_PCE_STORAGE(INSTMACRO, DS_ ## L ## _ ## S ## _ ## _ ## D, LO, GO, N)

#define INSTANTIATE_UQ_PCE_S_D(INSTMACRO, D, LO, GO, N) \
  INSTANTIATE_UQ_PCE_DS_SLD(INSTMACRO, double, int, D, LO, GO, N)

#define INSTANTIATE_UQ_PCE_S(INSTMACRO, LO, GO, N) \
  typedef Stokhos::DeviceForNode2<N>::type DFN_ ## LO ## _ ## GO ## _ ## N; \
  INSTANTIATE_UQ_PCE_S_D(INSTMACRO, DFN_ ## LO ## _ ## GO ## _ ## N, LO, GO, N)

#define INSTANTIATE_UQ_PCE(INSTMACRO, LO, GO, N) \
  INSTANTIATE_UQ_PCE_S(INSTMACRO, LO, GO, N)

#define INSTANTIATE_TPETRA_UQ_PCE_N(INSTMACRO, N)  \
  INSTANTIATE_UQ_PCE_S(INSTMACRO, int, int, N)

#if defined(HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT) && defined(HAVE_TPETRA_INST_SERIAL)
#define INSTANTIATE_TPETRA_UQ_PCE_SERIAL(INSTMACRO) \
  INSTMACRO(Kokkos_Compat_KokkosSerialWrapperNode)
#else
#define INSTANTIATE_TPETRA_UQ_PCE_SERIAL(INSTMACRO)
#endif

#if defined(HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT) && defined(HAVE_TPETRA_INST_PTHREAD)
#define INSTANTIATE_TPETRA_UQ_PCE_THREADS(INSTMACRO) \
  INSTMACRO(Kokkos_Compat_KokkosThreadsWrapperNode)
#else
#define INSTANTIATE_TPETRA_UQ_PCE_THREADS(INSTMACRO)
#endif

#if defined(HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT) && defined(HAVE_TPETRA_INST_OPENMP)
#define INSTANTIATE_TPETRA_UQ_PCE_OPENMP(INSTMACRO) \
  INSTMACRO(Kokkos_Compat_KokkosOpenMPWrapperNode)
#else
#define INSTANTIATE_TPETRA_UQ_PCE_OPENMP(INSTMACRO)
#endif

#if defined(HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT) && defined(HAVE_TPETRA_INST_CUDA)
#define INSTANTIATE_TPETRA_UQ_PCE_CUDA(INSTMACRO) \
  INSTMACRO(Kokkos_Compat_KokkosCudaWrapperNode)
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
