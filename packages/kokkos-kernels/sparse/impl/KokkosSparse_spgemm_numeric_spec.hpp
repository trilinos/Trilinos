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
#ifndef KOKKOSSPARSE_IMPL_SPGEMM_NUMERIC_SPEC_HPP_
#define KOKKOSSPARSE_IMPL_SPGEMM_NUMERIC_SPEC_HPP_

#include <KokkosKernels_config.h>

#include <Kokkos_Core.hpp>
// #include <Kokkos_ArithTraits.hpp>
#include "KokkosKernels_Handle.hpp"
// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include "KokkosSparse_spgemm_symbolic.hpp"
#include "KokkosSparse_spgemm_impl.hpp"
#include "KokkosSparse_spgemm_impl_seq.hpp"
#include "KokkosSparse_SortCrs.hpp"
#endif

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class KernelHandle, class a_size_view_t_, class a_lno_view_t, class a_scalar_view_t, class b_size_view_t_,
          class b_lno_view_t, class b_scalar_view_t, class c_size_view_t_, class c_lno_view_t, class c_scalar_view_t>
struct spgemm_numeric_eti_spec_avail {
  enum : bool { value = false };
};

}  // namespace Impl
}  // namespace KokkosSparse

#define KOKKOSSPARSE_SPGEMM_NUMERIC_ETI_SPEC_AVAIL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,          \
                                                   EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                              \
  template <>                                                                                                    \
  struct spgemm_numeric_eti_spec_avail<                                                                          \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                                                 \
    enum : bool { value = true };                                                                                \
  };

// Include the actual specialization declarations
#include <KokkosSparse_spgemm_numeric_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosSparse_spgemm_numeric_eti_spec_avail.hpp>

namespace KokkosSparse {
namespace Impl {

// Unification layer
/// \brief Implementation of KokkosBlas::spgemm (sparse matrix - dense
///   vector multiply) for multiple vectors at a time (multivectors)
///   and possibly multiple coefficients at a time.

template <class KernelHandle, class a_size_view_t_, class a_lno_view_t, class a_scalar_view_t, class b_size_view_t_,
          class b_lno_view_t, class b_scalar_view_t, class c_size_view_t_, class c_lno_view_t, class c_scalar_view_t,
          bool tpl_spec_avail = spgemm_numeric_tpl_spec_avail<
              KernelHandle, a_size_view_t_, a_lno_view_t, a_scalar_view_t, b_size_view_t_, b_lno_view_t,
              b_scalar_view_t, c_size_view_t_, c_lno_view_t, c_scalar_view_t>::value,
          bool eti_spec_avail = spgemm_numeric_eti_spec_avail<
              KernelHandle, a_size_view_t_, a_lno_view_t, a_scalar_view_t, b_size_view_t_, b_lno_view_t,
              b_scalar_view_t, c_size_view_t_, c_lno_view_t, c_scalar_view_t>::value>
struct SPGEMM_NUMERIC {
  static void spgemm_numeric(KernelHandle *handle, typename KernelHandle::nnz_lno_t m,
                             typename KernelHandle::nnz_lno_t n, typename KernelHandle::nnz_lno_t k,
                             a_size_view_t_ row_mapA, a_lno_view_t entriesA, a_scalar_view_t valuesA,

                             bool transposeA, b_size_view_t_ row_mapB, b_lno_view_t entriesB, b_scalar_view_t valuesB,
                             bool transposeB, c_size_view_t_ row_mapC, c_lno_view_t entriesC, c_scalar_view_t valuesC);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY

//! Full specialization of spgemm_mv for single vectors (2-D Views).
// Unification layer
template <class KernelHandle, class a_size_view_t_, class a_lno_view_t, class a_scalar_view_t, class b_size_view_t_,
          class b_lno_view_t, class b_scalar_view_t, class c_size_view_t_, class c_lno_view_t, class c_scalar_view_t>
struct SPGEMM_NUMERIC<KernelHandle, a_size_view_t_, a_lno_view_t, a_scalar_view_t, b_size_view_t_, b_lno_view_t,
                      b_scalar_view_t, c_size_view_t_, c_lno_view_t, c_scalar_view_t, false,
                      KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void spgemm_numeric(KernelHandle *handle, typename KernelHandle::nnz_lno_t m,
                             typename KernelHandle::nnz_lno_t n, typename KernelHandle::nnz_lno_t k,
                             a_size_view_t_ row_mapA, a_lno_view_t entriesA, a_scalar_view_t valuesA, bool transposeA,
                             b_size_view_t_ row_mapB, b_lno_view_t entriesB, b_scalar_view_t valuesB, bool transposeB,
                             c_size_view_t_ row_mapC, c_lno_view_t entriesC, c_scalar_view_t valuesC) {
    typedef typename KernelHandle::SPGEMMHandleType spgemmHandleType;
    spgemmHandleType *sh = handle->get_spgemm_handle();
    if (!sh->is_symbolic_called()) {
      throw std::runtime_error("Call spgemm symbolic before calling SpGEMM numeric");
    }
    if (m == 0 || n == 0 || k == 0 || !entriesA.extent(0) || !entriesB.extent(0)) {
      sh->set_call_numeric();
      sh->set_computed_entries();
      return;
    }
    switch (sh->get_algorithm_type()) {
      case SPGEMM_SERIAL:
      case SPGEMM_DEBUG:
        spgemm_debug_numeric(handle, m, n, k, row_mapA, entriesA, valuesA,

                             transposeA, row_mapB, entriesB, valuesB, transposeB, row_mapC, entriesC, valuesC);
        break;
      default: {
        KokkosSPGEMM<KernelHandle, a_size_view_t_, a_lno_view_t, a_scalar_view_t, b_size_view_t_, b_lno_view_t,
                     b_scalar_view_t>
            kspgemm(handle, m, n, k, row_mapA, entriesA, valuesA, transposeA, row_mapB, entriesB, valuesB, transposeB);
        kspgemm.KokkosSPGEMM_numeric(row_mapC, entriesC, valuesC);
      } break;
    }
    // Current implementation does not produce sorted matrix
    // TODO: remove this call when impl sorts
    KokkosSparse::sort_crs_matrix<typename KernelHandle::HandleExecSpace>(row_mapC, entriesC, valuesC);
    sh->set_call_numeric();
    sh->set_computed_entries();
  }
};

#endif

}  // namespace Impl
}  // namespace KokkosSparse

#define KOKKOSSPARSE_SPGEMM_NUMERIC_ETI_SPEC_DECL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,                \
                                                  EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                                    \
  extern template struct SPGEMM_NUMERIC<                                                                              \
      typename KokkosKernels::Experimental::KokkosKernelsHandle<                                                      \
          const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>, \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      false, true>;

#define KOKKOSSPARSE_SPGEMM_NUMERIC_ETI_SPEC_INST(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,           \
                                                  EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                               \
  template struct SPGEMM_NUMERIC<                                                                                \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      false, true>;

#include <KokkosSparse_spgemm_numeric_tpl_spec_decl.hpp>

#endif  // KOKKOS_BLAS1_MV_IMPL_DOT_HPP_
