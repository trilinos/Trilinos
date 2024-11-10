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
#ifndef KOKKOSSPARSE_IMPL_SPADD_NUMERIC_SPEC_HPP_
#define KOKKOSSPARSE_IMPL_SPADD_NUMERIC_SPEC_HPP_

#include <KokkosKernels_config.h>

#include <Kokkos_Core.hpp>
#include "KokkosKernels_Handle.hpp"
// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include "KokkosSparse_spadd_numeric_impl.hpp"
#endif

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class ExecSpace, class KernelHandle, class a_size_view_t_, class a_lno_view_t, class a_scalar_view_t,
          class b_size_view_t_, class b_lno_view_t, class b_scalar_view_t, class c_size_view_t_, class c_lno_view_t,
          class c_scalar_view_t>
struct spadd_numeric_eti_spec_avail {
  enum : bool { value = false };
};

}  // namespace Impl
}  // namespace KokkosSparse

#define KOKKOSSPARSE_SPADD_NUMERIC_ETI_SPEC_AVAIL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,           \
                                                  EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                               \
  template <>                                                                                                    \
  struct spadd_numeric_eti_spec_avail<                                                                           \
      EXEC_SPACE_TYPE,                                                                                           \
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
#include <KokkosSparse_spadd_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosSparse_spadd_numeric_eti_spec_avail.hpp>

namespace KokkosSparse {
namespace Impl {

// Unification layer
/// \brief Implementation of KokkosBlas::spadd (sparse-sparse matrix addition)

template <class ExecSpace, class KernelHandle, class a_size_view_t, class a_lno_view_t, class a_scalar_view_t,
          class b_size_view_t, class b_lno_view_t, class b_scalar_view_t, class c_size_view_t, class c_lno_view_t,
          class c_scalar_view_t,
          bool tpl_spec_avail = spadd_numeric_tpl_spec_avail<
              ExecSpace, KernelHandle, a_size_view_t, a_lno_view_t, a_scalar_view_t, b_size_view_t, b_lno_view_t,
              b_scalar_view_t, c_size_view_t, c_lno_view_t, c_scalar_view_t>::value,
          bool eti_spec_avail = spadd_numeric_eti_spec_avail<
              ExecSpace, KernelHandle, a_size_view_t, a_lno_view_t, a_scalar_view_t, b_size_view_t, b_lno_view_t,
              b_scalar_view_t, c_size_view_t, c_lno_view_t, c_scalar_view_t>::value>
struct SPADD_NUMERIC {
  static void spadd_numeric(const ExecSpace &exec, KernelHandle *handle, typename KernelHandle::const_nnz_lno_t m,
                            typename KernelHandle::const_nnz_lno_t n, typename a_scalar_view_t::const_value_type alpha,
                            a_size_view_t row_mapA, a_lno_view_t entriesA, a_scalar_view_t valuesA,
                            typename b_scalar_view_t::const_value_type beta, b_size_view_t row_mapB,
                            b_lno_view_t entriesB, b_scalar_view_t valuesB, c_size_view_t row_mapC,
                            c_lno_view_t entriesC, c_scalar_view_t valuesC);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY

template <class ExecSpace, class KernelHandle, class a_size_view_t, class a_lno_view_t, class a_scalar_view_t,
          class b_size_view_t, class b_lno_view_t, class b_scalar_view_t, class c_size_view_t, class c_lno_view_t,
          class c_scalar_view_t>
struct SPADD_NUMERIC<ExecSpace, KernelHandle, a_size_view_t, a_lno_view_t, a_scalar_view_t, b_size_view_t, b_lno_view_t,
                     b_scalar_view_t, c_size_view_t, c_lno_view_t, c_scalar_view_t, false,
                     KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void spadd_numeric(const ExecSpace &exec, KernelHandle *handle, typename KernelHandle::const_nnz_lno_t /* m */,
                            typename KernelHandle::const_nnz_lno_t /* n */,
                            typename a_scalar_view_t::const_value_type alpha, a_size_view_t row_mapA,
                            a_lno_view_t entriesA, a_scalar_view_t valuesA,
                            typename b_scalar_view_t::const_value_type beta, b_size_view_t row_mapB,
                            b_lno_view_t entriesB, b_scalar_view_t valuesB, c_size_view_t row_mapC,
                            c_lno_view_t entriesC, c_scalar_view_t valuesC) {
    spadd_numeric_impl(exec, handle, row_mapA, entriesA, valuesA, alpha, row_mapB, entriesB, valuesB, beta, row_mapC,
                       entriesC, valuesC);
  }
};

#endif

}  // namespace Impl
}  // namespace KokkosSparse

#define KOKKOSSPARSE_SPADD_NUMERIC_ETI_SPEC_DECL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, \
                                                 MEM_SPACE_TYPE)                                                       \
  extern template struct SPADD_NUMERIC<                                                                                \
      EXEC_SPACE_TYPE,                                                                                                 \
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
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      false, true>;

#define KOKKOSSPARSE_SPADD_NUMERIC_ETI_SPEC_INST(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, \
                                                 MEM_SPACE_TYPE)                                                       \
  template struct SPADD_NUMERIC<                                                                                       \
      EXEC_SPACE_TYPE,                                                                                                 \
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
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      false, true>;

#include <KokkosSparse_spadd_numeric_tpl_spec_decl.hpp>

#endif
