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
#ifndef KOKKOSSPARSE_IMPL_SPGEMM_JACOBI_SPEC_HPP_
#define KOKKOSSPARSE_IMPL_SPGEMM_JACOBI_SPEC_HPP_

#include <KokkosKernels_config.h>

#include <Kokkos_Core.hpp>

#include "KokkosKernels_Handle.hpp"
// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include "KokkosSparse_spgemm_symbolic.hpp"
#include "KokkosSparse_spgemm_impl.hpp"
#include "KokkosSparse_spgemm_jacobi_denseacc_impl.hpp"
#include "KokkosSparse_spgemm_jacobi_sparseacc_impl.hpp"
#include "KokkosSparse_spgemm_jacobi_seq_impl.hpp"
#include "KokkosSparse_SortCrs.hpp"
#endif

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class KernelHandle, class a_size_view_t_, class a_lno_view_t, class a_scalar_view_t, class b_size_view_t_,
          class b_lno_view_t, class b_scalar_view_t, class c_size_view_t_, class c_lno_view_t, class c_scalar_view_t,
          class dinv_scalar_view_t>
struct spgemm_jacobi_eti_spec_avail {
  enum : bool { value = false };
};

}  // namespace Impl
}  // namespace KokkosSparse

#define KOKKOSSPARSE_SPGEMM_JACOBI_ETI_SPEC_AVAIL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,           \
                                                  EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                               \
  template <>                                                                                                    \
  struct spgemm_jacobi_eti_spec_avail<                                                                           \
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
      Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                                                 \
    enum : bool { value = true };                                                                                \
  };

// Include the actual specialization declarations
#include <KokkosSparse_spgemm_jacobi_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosSparse_spgemm_jacobi_eti_spec_avail.hpp>

namespace KokkosSparse {
namespace Impl {

// Unification layer
template <class KernelHandle, class a_size_view_t_, class a_lno_view_t, class a_scalar_view_t, class b_size_view_t_,
          class b_lno_view_t, class b_scalar_view_t, class c_size_view_t_, class c_lno_view_t, class c_scalar_view_t,
          class dinv_scalar_view_t,
          bool tpl_spec_avail = spgemm_jacobi_tpl_spec_avail<
              KernelHandle, a_size_view_t_, a_lno_view_t, a_scalar_view_t, b_size_view_t_, b_lno_view_t,
              b_scalar_view_t, c_size_view_t_, c_lno_view_t, c_scalar_view_t, dinv_scalar_view_t>::value,
          bool eti_spec_avail = spgemm_jacobi_eti_spec_avail<
              KernelHandle, a_size_view_t_, a_lno_view_t, a_scalar_view_t, b_size_view_t_, b_lno_view_t,
              b_scalar_view_t, c_size_view_t_, c_lno_view_t, c_scalar_view_t, dinv_scalar_view_t>::value>

struct SPGEMM_JACOBI {
  static void spgemm_jacobi(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t m,
                            typename KernelHandle::const_nnz_lno_t n, typename KernelHandle::const_nnz_lno_t k,

                            a_size_view_t_ row_mapA, a_lno_view_t entriesA, a_scalar_view_t valuesA,

                            bool transposeA, b_size_view_t_ row_mapB, b_lno_view_t entriesB, b_scalar_view_t valuesB,

                            bool transposeB, c_size_view_t_ row_mapC, c_lno_view_t &entriesC, c_scalar_view_t &valuesC,

                            typename a_scalar_view_t::const_value_type omega, dinv_scalar_view_t dinv

  );
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY

//! Full specialization of spgemm jacobi
// Unification layer
template <class KernelHandle, class a_size_view_t_, class a_lno_view_t, class a_scalar_view_t, class b_size_view_t_,
          class b_lno_view_t, class b_scalar_view_t, class c_size_view_t_, class c_lno_view_t, class c_scalar_view_t,
          class dinv_scalar_view_t>
struct SPGEMM_JACOBI<KernelHandle, a_size_view_t_, a_lno_view_t, a_scalar_view_t, b_size_view_t_, b_lno_view_t,
                     b_scalar_view_t, c_size_view_t_, c_lno_view_t, c_scalar_view_t, dinv_scalar_view_t, false,
                     KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void spgemm_jacobi(KernelHandle *handle, typename KernelHandle::nnz_lno_t m,
                            typename KernelHandle::nnz_lno_t n, typename KernelHandle::nnz_lno_t k,

                            a_size_view_t_ row_mapA, a_lno_view_t entriesA, a_scalar_view_t valuesA,

                            bool transposeA, b_size_view_t_ row_mapB, b_lno_view_t entriesB, b_scalar_view_t valuesB,

                            bool transposeB, c_size_view_t_ row_mapC, c_lno_view_t &entriesC, c_scalar_view_t &valuesC,

                            typename c_scalar_view_t::const_value_type omega, dinv_scalar_view_t dinv) {
    typedef typename KernelHandle::SPGEMMHandleType spgemmHandleType;
    spgemmHandleType *sh = handle->get_spgemm_handle();
    if (!sh->is_symbolic_called()) {
      throw std::runtime_error(
          "KokkosSparse::spgemm_jacobi: must first call spgemm_symbolic with "
          "the same handle.");
    }
    if (!sh->are_rowptrs_computed()) {
      // Call symbolic, and make sure rowptrs are populated. This will not
      // duplicate any work if the user has already called symbolic. Must also
      // cast away constness of row_mapC.
      using c_size_view_t_nc = typename c_size_view_t_::non_const_type;
      using c_size_type      = typename c_size_view_t_::non_const_value_type;
      c_size_view_t_nc row_mapC_nc(const_cast<c_size_type *>(row_mapC.data()), row_mapC.extent(0));
      KokkosSparse::Experimental::spgemm_symbolic(handle, m, n, k, row_mapA, entriesA, transposeA, row_mapB, entriesB,
                                                  transposeB, row_mapC_nc, true);
    }
    if (!sh->are_rowflops_computed()) {
      KokkosSPGEMM<KernelHandle, a_size_view_t_, a_lno_view_t, a_scalar_view_t, b_size_view_t_, b_lno_view_t,
                   b_scalar_view_t>
          kspgemm(handle, m, n, k, row_mapA, entriesA, valuesA, transposeA, row_mapB, entriesB, valuesB, transposeB);
      kspgemm.compute_row_flops();
    }

    if (sh->get_algorithm_type() == SPGEMM_SERIAL) {
      spgemm_jacobi_seq(handle, m, n, k, row_mapA, entriesA, valuesA, transposeA, row_mapB, entriesB, valuesB,
                        transposeB, row_mapC, entriesC, valuesC, omega, dinv);
    } else {
      KokkosSPGEMM<KernelHandle, a_size_view_t_, a_lno_view_t, a_scalar_view_t, b_size_view_t_, b_lno_view_t,
                   b_scalar_view_t>
          kspgemm(handle, m, n, k, row_mapA, entriesA, valuesA, transposeA, row_mapB, entriesB, valuesB, transposeB);
      KokkosKernels::Impl::ExecSpaceType myExecSpace =
          KokkosKernels::Impl::get_exec_space_type<typename KernelHandle::HandleExecSpace>();

      kspgemm.KokkosSPGEMM_jacobi_sparseacc(row_mapC, entriesC, valuesC, omega, dinv, myExecSpace);
    }
    // Current implementation does not produce sorted matrix
    // TODO: remove this call when impl sorts
    KokkosSparse::sort_crs_matrix<typename KernelHandle::HandleExecSpace>(row_mapC, entriesC, valuesC);
  }
};

#endif

}  // namespace Impl
}  // namespace KokkosSparse

#define KOKKOSSPARSE_SPGEMM_JACOBI_ETI_SPEC_DECL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, \
                                                 MEM_SPACE_TYPE)                                                       \
  extern template struct SPGEMM_JACOBI<                                                                                \
      typename KokkosKernels::Experimental::KokkosKernelsHandle<                                                       \
          const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,  \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      false, true>;

#define KOKKOSSPARSE_SPGEMM_JACOBI_ETI_SPEC_INST(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, \
                                                 MEM_SPACE_TYPE)                                                       \
  template struct SPGEMM_JACOBI<                                                                                       \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE,       \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,               \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      false, true>;

#include <KokkosSparse_spgemm_jacobi_tpl_spec_decl.hpp>

#endif
