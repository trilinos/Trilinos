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
#ifndef KOKKOSSPARSE_IMPL_SPTRSV_SOLVE_SPEC_HPP_
#define KOKKOSSPARSE_IMPL_SPTRSV_SOLVE_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosKernels_Handle.hpp"

// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosSparse_sptrsv_solve_impl.hpp>
#include <KokkosSparse_sptrsv_symbolic_impl.hpp>
#endif

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class ExecutionSpace, class KernelHandle, class RowMapType, class EntriesType, class ValuesType, class BType,
          class XType>
struct sptrsv_solve_eti_spec_avail {
  enum : bool { value = false };
};

}  // namespace Impl
}  // namespace KokkosSparse

#define KOKKOSSPARSE_SPTRSV_SOLVE_ETI_SPEC_AVAIL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, \
                                                 MEM_SPACE_TYPE)                                                       \
  template <>                                                                                                          \
  struct sptrsv_solve_eti_spec_avail<                                                                                  \
      EXEC_SPACE_TYPE,                                                                                                 \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE,       \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,               \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                                   \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                                   \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                                   \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                                   \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                                                       \
    enum : bool { value = true };                                                                                      \
  };

// Include the actual specialization declarations
#include <KokkosSparse_sptrsv_solve_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosSparse_sptrsv_solve_eti_spec_avail.hpp>

namespace KokkosSparse {
namespace Impl {

#if defined(KOKKOS_ENABLE_CUDA) && 10000 < CUDA_VERSION && defined(KOKKOSKERNELS_ENABLE_EXP_CUDAGRAPH)
#define KOKKOSKERNELS_SPTRSV_CUDAGRAPHSUPPORT
#endif

// Unification layer
/// \brief Implementations of KokkosSparse::sptrsv_solve and
/// \brief KokkosSparse::sptrsv_solve_streams

template <class ExecutionSpace, class KernelHandle, class RowMapType, class EntriesType, class ValuesType, class BType,
          class XType,
          bool tpl_spec_avail = sptrsv_solve_tpl_spec_avail<ExecutionSpace, KernelHandle, RowMapType, EntriesType,
                                                            ValuesType, BType, XType>::value,
          bool eti_spec_avail = sptrsv_solve_eti_spec_avail<ExecutionSpace, KernelHandle, RowMapType, EntriesType,
                                                            ValuesType, BType, XType>::value>
struct SPTRSV_SOLVE {
  static void sptrsv_solve(ExecutionSpace &space, KernelHandle *handle, const RowMapType row_map,
                           const EntriesType entries, const ValuesType values, BType b, XType x);

  static void sptrsv_solve_streams(const std::vector<ExecutionSpace> &execspace_v, std::vector<KernelHandle> &handle_v,
                                   const std::vector<RowMapType> &row_map_v, const std::vector<EntriesType> &entries_v,
                                   const std::vector<ValuesType> &values_v, const std::vector<BType> &b_v,
                                   std::vector<XType> &x_v);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//! Full specialization of sptrsv_solve and sptrsv_solve_streams
// Unification layer
template <class ExecutionSpace, class KernelHandle, class RowMapType, class EntriesType, class ValuesType, class BType,
          class XType>
struct SPTRSV_SOLVE<ExecutionSpace, KernelHandle, RowMapType, EntriesType, ValuesType, BType, XType, false,
                    KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void sptrsv_solve(ExecutionSpace &space, KernelHandle *handle, const RowMapType row_map,
                           const EntriesType entries, const ValuesType values, BType b, XType x) {
    using Sptrsv = Experimental::SptrsvWrap<typename KernelHandle::SPTRSVHandleType>;

    // Call specific algorithm type
    auto sptrsv_handle       = handle->get_sptrsv_handle();
    const auto block_enabled = sptrsv_handle->is_block_enabled();
    Kokkos::Profiling::pushRegion(sptrsv_handle->is_lower_tri() ? "KokkosSparse_sptrsv[lower]"
                                                                : "KokkosSparse_sptrsv[upper]");
    if (sptrsv_handle->is_lower_tri()) {
      if (sptrsv_handle->is_symbolic_complete() == false) {
        Experimental::lower_tri_symbolic(space, *sptrsv_handle, row_map, entries);
      }
      if (sptrsv_handle->get_algorithm() == KokkosSparse::Experimental::SPTRSVAlgorithm::SEQLVLSCHD_TP1CHAIN) {
        Sptrsv::template tri_solve_chain<true>(space, *sptrsv_handle, row_map, entries, values, b, x);
      } else {
#ifdef KOKKOSKERNELS_SPTRSV_CUDAGRAPHSUPPORT
        using ExecSpace = typename RowMapType::memory_space::execution_space;
        if (std::is_same<ExecSpace, Kokkos::Cuda>::value)
          // TODO: set stream in thandle's sptrsvCudaGraph
          Sptrsv::tri_solve_cg<true>(*sptrsv_handle, row_map, entries, values, b, x);
        else
#endif
        {
          if (block_enabled) {
            Sptrsv::template lower_tri_solve<true>(space, *sptrsv_handle, row_map, entries, values, b, x);
          } else {
            Sptrsv::template lower_tri_solve<false>(space, *sptrsv_handle, row_map, entries, values, b, x);
          }
        }
      }
    } else {
      if (sptrsv_handle->is_symbolic_complete() == false) {
        Experimental::upper_tri_symbolic(space, *sptrsv_handle, row_map, entries);
      }
      if (sptrsv_handle->get_algorithm() == KokkosSparse::Experimental::SPTRSVAlgorithm::SEQLVLSCHD_TP1CHAIN) {
        Sptrsv::template tri_solve_chain<false>(space, *sptrsv_handle, row_map, entries, values, b, x);
      } else {
#ifdef KOKKOSKERNELS_SPTRSV_CUDAGRAPHSUPPORT
        using ExecSpace = typename RowMapType::memory_space::execution_space;
        if (std::is_same<ExecSpace, Kokkos::Cuda>::value)
          // TODO: set stream in thandle's sptrsvCudaGraph
          Sptrsv::tri_solve_cg<false>(*sptrsv_handle, row_map, entries, values, b, x);
        else
#endif
        {
          if (block_enabled) {
            Sptrsv::template upper_tri_solve<true>(space, *sptrsv_handle, row_map, entries, values, b, x);
          } else {
            Sptrsv::template upper_tri_solve<false>(space, *sptrsv_handle, row_map, entries, values, b, x);
          }
        }
      }
    }
    Kokkos::Profiling::popRegion();
  }

  static void sptrsv_solve_streams(const std::vector<ExecutionSpace> &execspace_v, std::vector<KernelHandle> &handle_v,
                                   const std::vector<RowMapType> &row_map_v, const std::vector<EntriesType> &entries_v,
                                   const std::vector<ValuesType> &values_v, const std::vector<BType> &b_v,
                                   std::vector<XType> &x_v) {
    using Sptrsv = Experimental::SptrsvWrap<typename KernelHandle::SPTRSVHandleType>;
    // Call specific algorithm type
    // NOTE: Only support SEQLVLSCHD_RP and SEQLVLSCHD_TP1 at this moment
    //       Assume streams have the same either lower or upper matrix type
    std::vector<typename KernelHandle::SPTRSVHandleType *> sptrsv_handle_v(execspace_v.size());
    for (int i = 0; i < static_cast<int>(execspace_v.size()); i++) {
      sptrsv_handle_v[i] = handle_v[i].get_sptrsv_handle();
    }
    Kokkos::Profiling::pushRegion(sptrsv_handle_v[0]->is_lower_tri() ? "KokkosSparse_sptrsv[lower]"
                                                                     : "KokkosSparse_sptrsv[upper]");
    if (sptrsv_handle_v[0]->is_lower_tri()) {
      for (int i = 0; i < static_cast<int>(execspace_v.size()); i++) {
        if (sptrsv_handle_v[i]->is_symbolic_complete() == false) {
          Experimental::lower_tri_symbolic(execspace_v[i], *(sptrsv_handle_v[i]), row_map_v[i], entries_v[i]);
        }
      }
      Sptrsv::template tri_solve_streams<true>(execspace_v, sptrsv_handle_v, row_map_v, entries_v, values_v, b_v, x_v);
    } else {
      for (int i = 0; i < static_cast<int>(execspace_v.size()); i++) {
        if (sptrsv_handle_v[i]->is_symbolic_complete() == false) {
          Experimental::upper_tri_symbolic(execspace_v[i], *(sptrsv_handle_v[i]), row_map_v[i], entries_v[i]);
        }
      }
      Sptrsv::template tri_solve_streams<false>(execspace_v, sptrsv_handle_v, row_map_v, entries_v, values_v, b_v, x_v);
    }
    Kokkos::Profiling::popRegion();
  }
};

#endif
}  // namespace Impl
}  // namespace KokkosSparse

//
// Macro for declaration of full specialization of
// This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSSPARSE_SPTRSV_SOLVE_ETI_SPEC_DECL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, \
                                                MEM_SPACE_TYPE)                                                       \
  extern template struct SPTRSV_SOLVE<                                                                                \
      EXEC_SPACE_TYPE,                                                                                                \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE,      \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,              \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                                  \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                                  \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                                  \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                                  \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      false, true>;

#define KOKKOSSPARSE_SPTRSV_SOLVE_ETI_SPEC_INST(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, \
                                                MEM_SPACE_TYPE)                                                       \
  template struct SPTRSV_SOLVE<                                                                                       \
      EXEC_SPACE_TYPE,                                                                                                \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE,      \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,              \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                                  \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                                  \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                                  \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                                  \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      false, true>;

#include <KokkosSparse_sptrsv_solve_tpl_spec_decl.hpp>

#endif
