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
#ifndef KOKKOSSPARSE_IMPL_COLOR_D1_SPEC_HPP_
#define KOKKOSSPARSE_IMPL_COLOR_D1_SPEC_HPP_

#include <KokkosKernels_config.h>

#include <Kokkos_Core.hpp>
#include "KokkosKernels_Handle.hpp"
// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include "KokkosGraph_Distance1Color_impl.hpp"
#endif

namespace KokkosGraph {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class KernelHandle, class size_view_t_, class lno_view_t>
struct color_d1_eti_spec_avail {
  enum : bool { value = false };
};

}  // namespace Impl
}  // namespace KokkosGraph

#define KOKKOSGRAPH_COLOR_D1_ETI_SPEC_AVAIL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, \
                                            MEM_SPACE_TYPE)                                                       \
  template <>                                                                                                     \
  struct color_d1_eti_spec_avail<                                                                                 \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE,  \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,          \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,             \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                      \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {                                                    \
    enum : bool { value = true };                                                                                 \
  };

// Include the actual specialization declarations
#include <generated_specializations_hpp/KokkosGraph_color_d1_eti_spec_avail.hpp>

namespace KokkosGraph {
namespace Impl {

// Unification layer
/// \brief Implementation of KokkosGraph::graph_color (distance-1 greedy
/// coloring)

template <class KernelHandle, class size_view_t, class lno_view_t, bool tpl_spec_avail = false,
          bool eti_spec_avail = color_d1_eti_spec_avail<KernelHandle, size_view_t, lno_view_t>::value>
struct COLOR_D1 {
  static void color_d1(KernelHandle *handle, typename lno_view_t::non_const_value_type num_rows, size_view_t rowmap,
                       lno_view_t entries);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY

template <class KernelHandle, class size_view_t, class lno_view_t>
struct COLOR_D1<KernelHandle, size_view_t, lno_view_t, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void color_d1(KernelHandle *handle, typename lno_view_t::non_const_value_type num_rows, size_view_t rowmap,
                       lno_view_t entries) {
    KokkosGraph::Impl::graph_color_impl(handle, num_rows, rowmap, entries);
  }
};

#endif

}  // namespace Impl
}  // namespace KokkosGraph

#define KOKKOSGRAPH_COLOR_D1_ETI_SPEC_DECL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE,      \
                                           MEM_SPACE_TYPE)                                                            \
  extern template struct COLOR_D1<                                                                                    \
      typename KokkosKernels::Experimental::KokkosKernelsHandle<                                                      \
          const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>, \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                          \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                          \
      false, true>;

#define KOKKOSGRAPH_COLOR_D1_ETI_SPEC_INST(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, \
                                           MEM_SPACE_TYPE)                                                       \
  template struct COLOR_D1<                                                                                      \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                     \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                     \
      false, true>;

#endif
