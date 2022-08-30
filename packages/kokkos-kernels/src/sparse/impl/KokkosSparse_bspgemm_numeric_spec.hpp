/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#ifndef KOKKOSSPARSE_IMPL_BSPGEMM_NUMERIC_SPEC_HPP_
#define KOKKOSSPARSE_IMPL_BSPGEMM_NUMERIC_SPEC_HPP_

#include <KokkosKernels_config.h>

#include <Kokkos_Core.hpp>
//#include <Kokkos_ArithTraits.hpp>
#include "KokkosKernels_Handle.hpp"
// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//#include "KokkosSparse_spgemm_symbolic.hpp"
#include "KokkosSparse_spgemm_cuSPARSE_impl.hpp"
#include "KokkosSparse_spgemm_CUSP_impl.hpp"
#include "KokkosSparse_bspgemm_impl.hpp"
#include "KokkosSparse_bspgemm_impl_seq.hpp"
#include "KokkosSparse_spgemm_mkl_impl.hpp"
#include "KokkosSparse_spgemm_mkl2phase_impl.hpp"
#include "KokkosSparse_spgemm_viennaCL_impl.hpp"
#endif

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class KernelHandle, class a_size_view_t_, class a_lno_view_t,
          class a_scalar_view_t, class b_size_view_t_, class b_lno_view_t,
          class b_scalar_view_t, class c_size_view_t_, class c_lno_view_t,
          class c_scalar_view_t>
struct bspgemm_numeric_eti_spec_avail {
  enum : bool { value = false };
};

}  // namespace Impl
}  // namespace KokkosSparse

#define KOKKOSSPARSE_BSPGEMM_NUMERIC_ETI_SPEC_AVAIL(                      \
    SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, \
    FAST_MEM_SPACE_TYPE, SLOW_MEM_SPACE_TYPE)                             \
  template <>                                                             \
  struct bspgemm_numeric_eti_spec_avail<                                  \
      KokkosKernels::Experimental::KokkosKernelsHandle<                   \
          const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE,       \
          EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE, SLOW_MEM_SPACE_TYPE>,     \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE,                     \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE,                     \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE,                           \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE,                            \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {          \
    enum : bool { value = true };                                         \
  };                                                                      \
                                                                          \
  template <>                                                             \
  struct bspgemm_numeric_eti_spec_avail<                                  \
      KokkosKernels::Experimental::KokkosKernelsHandle<                   \
          const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE,       \
          EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE, SLOW_MEM_SPACE_TYPE>,     \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE,                     \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE,                     \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE,                            \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE,                           \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE,                            \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {          \
    enum : bool { value = true };                                         \
  };

// Include the actual specialization declarations
//#include <KokkosSparse_bspgemm_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosSparse_bspgemm_numeric_eti_spec_avail.hpp>

namespace KokkosSparse {
namespace Impl {

// For future use (when TPL with block SpGEMM numeric phase is encountered)
template <class KernelHandle, class a_size_view_t_, class a_lno_view_t,
          class a_scalar_view_t, class b_size_view_t_, class b_lno_view_t,
          class b_scalar_view_t, class c_size_view_t_, class c_lno_view_t,
          class c_scalar_view_t>
struct bspgemm_numeric_tpl_spec_avail {
  enum : bool { value = false };
};

// Unification layer
/// \brief Implementation of BSR sparse block matrix - matrix multiplication

template <class KernelHandle, class a_size_view_t_, class a_lno_view_t,
          class a_scalar_view_t, class b_size_view_t_, class b_lno_view_t,
          class b_scalar_view_t, class c_size_view_t_, class c_lno_view_t,
          class c_scalar_view_t,
          bool tpl_spec_avail = bspgemm_numeric_tpl_spec_avail<
              KernelHandle, a_size_view_t_, a_lno_view_t, a_scalar_view_t,
              b_size_view_t_, b_lno_view_t, b_scalar_view_t, c_size_view_t_,
              c_lno_view_t, c_scalar_view_t>::value,
          bool eti_spec_avail = bspgemm_numeric_eti_spec_avail<
              KernelHandle, a_size_view_t_, a_lno_view_t, a_scalar_view_t,
              b_size_view_t_, b_lno_view_t, b_scalar_view_t, c_size_view_t_,
              c_lno_view_t, c_scalar_view_t>::value>
struct BSPGEMM_NUMERIC {
  static void bspgemm_numeric(KernelHandle *handle,
                              typename KernelHandle::const_nnz_lno_t m,
                              typename KernelHandle::const_nnz_lno_t n,
                              typename KernelHandle::const_nnz_lno_t k,
                              typename KernelHandle::const_nnz_lno_t blockDim,
                              a_size_view_t_ row_mapA, a_lno_view_t entriesA,
                              a_scalar_view_t valuesA,

                              bool transposeA, b_size_view_t_ row_mapB,
                              b_lno_view_t entriesB, b_scalar_view_t valuesB,
                              bool transposeB, c_size_view_t_ row_mapC,
                              c_lno_view_t &entriesC, c_scalar_view_t &valuesC);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY

//! Full specialization of block spgemm
// Unification layer
template <class KernelHandle, class a_size_view_t_, class a_lno_view_t,
          class a_scalar_view_t, class b_size_view_t_, class b_lno_view_t,
          class b_scalar_view_t, class c_size_view_t_, class c_lno_view_t,
          class c_scalar_view_t>
struct BSPGEMM_NUMERIC<
    KernelHandle, a_size_view_t_, a_lno_view_t, a_scalar_view_t, b_size_view_t_,
    b_lno_view_t, b_scalar_view_t, c_size_view_t_, c_lno_view_t,
    c_scalar_view_t, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void bspgemm_numeric(
      KernelHandle *handle, typename KernelHandle::nnz_lno_t m,
      typename KernelHandle::nnz_lno_t n, typename KernelHandle::nnz_lno_t k,
      typename KernelHandle::const_nnz_lno_t blockDim, a_size_view_t_ row_mapA,
      a_lno_view_t entriesA, a_scalar_view_t valuesA,

      bool transposeA, b_size_view_t_ row_mapB, b_lno_view_t entriesB,
      b_scalar_view_t valuesB, bool transposeB, c_size_view_t_ row_mapC,
      c_lno_view_t &entriesC, c_scalar_view_t &valuesC) {
    typedef typename KernelHandle::SPGEMMHandleType spgemmHandleType;
    spgemmHandleType *sh = handle->get_spgemm_handle();
    if (!sh->is_symbolic_called()) {
      throw std::runtime_error(
          "Call spgemm symbolic before calling SpGEMM numeric");
    }

    switch (sh->get_algorithm_type()) {
      case SPGEMM_CUSPARSE:
        throw std::runtime_error(
            "cuSPARSE implementation for block SpGEMM is not available");
      case SPGEMM_CUSP:
        throw std::runtime_error(
            "CUSP implementation for block SpGEMM is not available");
      case SPGEMM_MKL:
      case SPGEMM_MKL2PHASE:
        throw std::runtime_error(
            "MKL implementation available for block SpGEMM is not available");
      case SPGEMM_VIENNA:
        throw std::runtime_error(
            "Vienna implementation available for block SpGEMM is not "
            "available");

      default:

      {
        KokkosBSPGEMM<KernelHandle, a_size_view_t_, a_lno_view_t,
                      a_scalar_view_t, b_size_view_t_, b_lno_view_t,
                      b_scalar_view_t>
            kbspgemm(handle, m, n, k, blockDim, row_mapA, entriesA, valuesA,
                     transposeA, row_mapB, entriesB, valuesB, transposeB);
        kbspgemm.KokkosBSPGEMM_numeric(row_mapC, entriesC, valuesC);
      } break;
      case SPGEMM_SERIAL:
      case SPGEMM_DEBUG:
        bspgemm_debug_numeric(handle, m, n, k, blockDim, row_mapA, entriesA,
                              valuesA, transposeA, row_mapB, entriesB, valuesB,
                              transposeB, row_mapC, entriesC, valuesC);
        break;
    }
  }
};

#endif

}  // namespace Impl
}  // namespace KokkosSparse

#define KOKKOSSPARSE_BSPGEMM_NUMERIC_ETI_SPEC_DECL(                       \
    SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, \
    FAST_MEM_SPACE_TYPE, SLOW_MEM_SPACE_TYPE)                             \
  extern template struct BSPGEMM_NUMERIC<                                 \
      typename KokkosKernels::Experimental::KokkosKernelsHandle<          \
          const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE,       \
          EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE, SLOW_MEM_SPACE_TYPE>,     \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE,                     \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE,                     \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE,                           \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE,                            \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      false, true>;                                                       \
                                                                          \
  extern template struct BSPGEMM_NUMERIC<                                 \
      typename KokkosKernels::Experimental::KokkosKernelsHandle<          \
          const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE,       \
          EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE, SLOW_MEM_SPACE_TYPE>,     \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE,                     \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE,                     \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE,                            \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE,                           \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE,                            \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      false, true>;

#define KOKKOSSPARSE_BSPGEMM_NUMERIC_ETI_SPEC_INST(                       \
    SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, \
    FAST_MEM_SPACE_TYPE, SLOW_MEM_SPACE_TYPE)                             \
  template struct BSPGEMM_NUMERIC<                                        \
      KokkosKernels::Experimental::KokkosKernelsHandle<                   \
          const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE,       \
          EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE, SLOW_MEM_SPACE_TYPE>,     \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE,                     \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE,                     \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE,                           \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE,                            \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      false, true>;                                                       \
                                                                          \
  template struct BSPGEMM_NUMERIC<                                        \
      KokkosKernels::Experimental::KokkosKernelsHandle<                   \
          const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE,       \
          EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE, SLOW_MEM_SPACE_TYPE>,     \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE,                     \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE,                     \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE,                            \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE,                           \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE,                            \
                   Kokkos::Device<EXEC_SPACE_TYPE, FAST_MEM_SPACE_TYPE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      false, true>;

//#include <KokkosSparse_spgemm_tpl_spec_decl.hpp>
#include <generated_specializations_hpp/KokkosSparse_bspgemm_numeric_eti_spec_decl.hpp>

#endif  // KOKKOS_BLAS1_MV_IMPL_DOT_HPP_
