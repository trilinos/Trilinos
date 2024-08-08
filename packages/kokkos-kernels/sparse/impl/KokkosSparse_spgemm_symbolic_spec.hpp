//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
#ifndef KOKKOSSPARSE_IMPL_SPGEMM_SYMBOLIC_SPEC_HPP_
#define KOKKOSSPARSE_IMPL_SPGEMM_SYMBOLIC_SPEC_HPP_

#include <KokkosKernels_config.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include "KokkosKernels_Handle.hpp"
// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include "KokkosSparse_spgemm_impl.hpp"
#include "KokkosSparse_spgemm_impl_seq.hpp"
#endif

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class KernelHandle, class a_size_view_t_, class a_lno_view_t, class b_size_view_t_, class b_lno_view_t,
          class c_size_view_t_>
struct spgemm_symbolic_eti_spec_avail {
  enum : bool { value = false };
};

}  // namespace Impl
}  // namespace KokkosSparse

#define KOKKOSSPARSE_SPGEMM_SYMBOLIC_ETI_SPEC_AVAIL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,         \
                                                    EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                             \
  template <>                                                                                                    \
  struct spgemm_symbolic_eti_spec_avail<                                                                         \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                                                 \
    enum : bool { value = true };                                                                                \
  };

// Include the actual specialization declarations
#include <KokkosSparse_spgemm_symbolic_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosSparse_spgemm_symbolic_eti_spec_avail.hpp>

namespace KokkosSparse {
namespace Impl {

// Unification layer
/// \brief Implementation of KokkosSparse::spgemm (sparse matrix - sparse
///   matrix multiply)
///
template <class KernelHandle, class a_size_view_t_, class a_lno_view_t, class b_size_view_t_, class b_lno_view_t,
          class c_size_view_t_,
          bool tpl_spec_avail = spgemm_symbolic_tpl_spec_avail<KernelHandle, a_size_view_t_, a_lno_view_t,
                                                               b_size_view_t_, b_lno_view_t, c_size_view_t_>::value,
          bool eti_spec_avail = spgemm_symbolic_eti_spec_avail<KernelHandle, a_size_view_t_, a_lno_view_t,
                                                               b_size_view_t_, b_lno_view_t, c_size_view_t_>::value>
struct SPGEMM_SYMBOLIC {
  static void spgemm_symbolic(KernelHandle *handle, typename KernelHandle::nnz_lno_t m,
                              typename KernelHandle::nnz_lno_t n, typename KernelHandle::nnz_lno_t k,
                              a_size_view_t_ row_mapA, a_lno_view_t entriesA, bool transposeA, b_size_view_t_ row_mapB,
                              b_lno_view_t entriesB, bool transposeB, c_size_view_t_ row_mapC, bool computeRowptrs);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//! Full specialization of spgemm for single vectors (1-D Views).
// Unification layer
template <class KernelHandle, class a_size_view_t_, class a_lno_view_t, class b_size_view_t_, class b_lno_view_t,
          class c_size_view_t_>
struct SPGEMM_SYMBOLIC<KernelHandle, a_size_view_t_, a_lno_view_t, b_size_view_t_, b_lno_view_t, c_size_view_t_, false,
                       KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void spgemm_symbolic(KernelHandle *handle, typename KernelHandle::nnz_lno_t m,
                              typename KernelHandle::nnz_lno_t n, typename KernelHandle::nnz_lno_t k,
                              a_size_view_t_ row_mapA, a_lno_view_t entriesA, bool transposeA, b_size_view_t_ row_mapB,
                              b_lno_view_t entriesB, bool transposeB, c_size_view_t_ row_mapC,
                              bool /* computeRowptrs */) {
    typedef typename KernelHandle::SPGEMMHandleType spgemmHandleType;
    spgemmHandleType *sh = handle->get_spgemm_handle();
    if (sh->is_symbolic_called() && sh->are_rowptrs_computed()) return;
    if (m == 0 || n == 0 || k == 0 || !entriesA.extent(0) || !entriesB.extent(0)) {
      sh->set_computed_rowptrs();
      sh->set_call_symbolic();
      sh->set_c_nnz(0);
      Kokkos::deep_copy(typename spgemmHandleType::HandleExecSpace(), row_mapC,
                        typename c_size_view_t_::non_const_value_type(0));
      return;
    }
    switch (sh->get_algorithm_type()) {
      case SPGEMM_SERIAL:
      case SPGEMM_DEBUG:
        spgemm_debug_symbolic(handle, m, n, k, row_mapA, entriesA,

                              transposeA, row_mapB, entriesB, transposeB, row_mapC);
        break;
      default: {
        KokkosSPGEMM<KernelHandle, a_size_view_t_, a_lno_view_t, typename KernelHandle::in_scalar_nnz_view_t,
                     b_size_view_t_, b_lno_view_t, typename KernelHandle::in_scalar_nnz_view_t>
            kspgemm(handle, m, n, k, row_mapA, entriesA, transposeA, row_mapB, entriesB, transposeB);
        kspgemm.KokkosSPGEMM_symbolic(row_mapC);
      } break;
    }
    sh->set_call_symbolic();
    // The KokkosKernels implementation of symbolic always populates rowptrs.
    sh->set_computed_rowptrs();
  }
};

#endif

}  // namespace Impl
}  // namespace KokkosSparse

//
// Macro for declaration of full specialization of
// KokkosSparse::Impl::SPGEMM_SYMBOLIC.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
//
#define KOKKOSSPARSE_SPGEMM_SYMBOLIC_ETI_SPEC_DECL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,          \
                                                   EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                              \
  extern template struct SPGEMM_SYMBOLIC<                                                                        \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      false, true>;

#define KOKKOSSPARSE_SPGEMM_SYMBOLIC_ETI_SPEC_INST(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,          \
                                                   EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                              \
  template struct SPGEMM_SYMBOLIC<                                                                               \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      false, true>;

#include <KokkosSparse_spgemm_symbolic_tpl_spec_decl.hpp>

#endif  // KOKKOS_BLAS1_MV_IMPL_DOT_HPP_
