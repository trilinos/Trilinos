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
#ifndef __KOKKOSBATCHED_HOSTLEVEL_GEMM_SPEC_HPP__
#define __KOKKOSBATCHED_HOSTLEVEL_GEMM_SPEC_HPP__

#include <Kokkos_Core.hpp>
#include <KokkosKernels_config.h>
#include <KokkosBatched_HostLevel_Gemm_Handle.hpp>  // BatchedGemmHandle

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include "KokkosBatched_HostLevel_Gemm_Impl.hpp"
#include "KokkosKernels_ExecSpaceUtils.hpp"
#endif

namespace KokkosBatched {
namespace Impl {
// Specialization struct which defines whether a specialization exists
// This struct is currently never specialized.
template <class ArgTransA, class ArgTransB, class ArgBatchSzDim, class BatchedGemmHandleType, class ScalarType,
          class AViewType, class BViewType, class CViewType>
struct batched_gemm_tpl_spec_avail {
  enum : bool { value = false };
};

// Specialization struct which defines whether a specialization exists
template <class ArgTransA, class ArgTransB, class ArgBatchSzDim, class BatchedGemmHandleType, class ScalarType,
          class AViewType, class BViewType, class CViewType>
struct batched_gemm_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBatched

// ETI specalization macros, consumed by generated *_eti_spec_avail.hpp files
#define KOKKOSBATCHED_GEMM_ETI_SPEC_AVAIL_INNER(ARG_TRANS_A, ARG_TRANS_B, ARG_BATCH_LAYOUT, SCALAR, LAYOUT,  \
                                                EXEC_SPACE, MEM_SPACE)                                       \
  template <>                                                                                                \
  struct batched_gemm_eti_spec_avail<ARG_TRANS_A, ARG_TRANS_B, ARG_BATCH_LAYOUT, BatchedGemmHandle, SCALAR,  \
                                     Kokkos::View<SCALAR ***, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                                                  Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                  \
                                     Kokkos::View<SCALAR ***, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                                                  Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                  \
                                     Kokkos::View<SCALAR ***, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                                                  Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {                \
    enum : bool { value = true };                                                                            \
  };

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
#define KOKKOSBATCHED_GEMM_BLL_ETI_SPEC_AVAIL_INNER(ARG_TRANS_A, ARG_TRANS_B, ARG_BATCH_LAYOUT, SCALAR, LAYOUT,    \
                                                    EXEC_SPACE, MEM_SPACE)                                         \
  KOKKOSBATCHED_GEMM_ETI_SPEC_AVAIL_INNER(ARG_TRANS_A, ARG_TRANS_B, ARG_BATCH_LAYOUT, SCALAR, Kokkos::LayoutRight, \
                                          EXEC_SPACE, MEM_SPACE)
#else
#define KOKKOSBATCHED_GEMM_BLL_ETI_SPEC_AVAIL_INNER(ARG_TRANS_A, ARG_TRANS_B, ARG_BATCH_LAYOUT, SCALAR, LAYOUT, \
                                                    EXEC_SPACE, MEM_SPACE)
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
#define KOKKOSBATCHED_GEMM_BLR_ETI_SPEC_AVAIL_INNER(ARG_TRANS_A, ARG_TRANS_B, ARG_BATCH_LAYOUT, SCALAR, LAYOUT,   \
                                                    EXEC_SPACE, MEM_SPACE)                                        \
  KOKKOSBATCHED_GEMM_ETI_SPEC_AVAIL_INNER(ARG_TRANS_A, ARG_TRANS_B, ARG_BATCH_LAYOUT, SCALAR, Kokkos::LayoutLeft, \
                                          EXEC_SPACE, MEM_SPACE)
#else
#define KOKKOSBATCHED_GEMM_BLR_ETI_SPEC_AVAIL_INNER(ARG_TRANS_A, ARG_TRANS_B, ARG_BATCH_LAYOUT, SCALAR, LAYOUT, \
                                                    EXEC_SPACE, MEM_SPACE)
#endif

///////////////// BatchLayout::Left Permutations /////////////////
#define KOKKOSBATCHED_GEMM_NT_NT_BLL_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                       \
  KOKKOSBATCHED_GEMM_BLL_ETI_SPEC_AVAIL_INNER(Trans::NoTranspose, Trans::NoTranspose, BatchLayout::Left, SCALAR, \
                                              LAYOUT, EXEC_SPACE, MEM_SPACE)

#define KOKKOSBATCHED_GEMM_NT_T_BLL_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                              \
  KOKKOSBATCHED_GEMM_BLL_ETI_SPEC_AVAIL_INNER(Trans::NoTranspose, Trans::Transpose, BatchLayout::Left, SCALAR, LAYOUT, \
                                              EXEC_SPACE, MEM_SPACE)

#define KOKKOSBATCHED_GEMM_T_NT_BLL_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                              \
  KOKKOSBATCHED_GEMM_BLL_ETI_SPEC_AVAIL_INNER(Trans::Transpose, Trans::NoTranspose, BatchLayout::Left, SCALAR, LAYOUT, \
                                              EXEC_SPACE, MEM_SPACE)

#define KOKKOSBATCHED_GEMM_T_T_BLL_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                             \
  KOKKOSBATCHED_GEMM_BLL_ETI_SPEC_AVAIL_INNER(Trans::Transpose, Trans::Transpose, BatchLayout::Left, SCALAR, LAYOUT, \
                                              EXEC_SPACE, MEM_SPACE)

// Include the BLL ETI specalizations
#include <generated_specializations_hpp/KokkosBatched_Gemm_nt_nt_bll_eti_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBatched_Gemm_nt_t_bll_eti_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBatched_Gemm_t_nt_bll_eti_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBatched_Gemm_t_t_bll_eti_spec_avail.hpp>

///////////////// BatchLayout::Right Permutations /////////////////
#define KOKKOSBATCHED_GEMM_NT_NT_BLR_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                        \
  KOKKOSBATCHED_GEMM_BLR_ETI_SPEC_AVAIL_INNER(Trans::NoTranspose, Trans::NoTranspose, BatchLayout::Right, SCALAR, \
                                              LAYOUT, EXEC_SPACE, MEM_SPACE)

#define KOKKOSBATCHED_GEMM_NT_T_BLR_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                       \
  KOKKOSBATCHED_GEMM_BLR_ETI_SPEC_AVAIL_INNER(Trans::NoTranspose, Trans::Transpose, BatchLayout::Right, SCALAR, \
                                              LAYOUT, EXEC_SPACE, MEM_SPACE)

#define KOKKOSBATCHED_GEMM_T_NT_BLR_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                       \
  KOKKOSBATCHED_GEMM_BLR_ETI_SPEC_AVAIL_INNER(Trans::Transpose, Trans::NoTranspose, BatchLayout::Right, SCALAR, \
                                              LAYOUT, EXEC_SPACE, MEM_SPACE)

#define KOKKOSBATCHED_GEMM_T_T_BLR_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                              \
  KOKKOSBATCHED_GEMM_BLR_ETI_SPEC_AVAIL_INNER(Trans::Transpose, Trans::Transpose, BatchLayout::Right, SCALAR, LAYOUT, \
                                              EXEC_SPACE, MEM_SPACE)

// Include the BLR ETI specalizations
#include <generated_specializations_hpp/KokkosBatched_Gemm_nt_nt_blr_eti_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBatched_Gemm_nt_t_blr_eti_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBatched_Gemm_t_nt_blr_eti_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBatched_Gemm_t_t_blr_eti_spec_avail.hpp>

namespace KokkosBatched {
namespace Impl {
template <class ArgTransA, class ArgTransB, class ArgBatchSzDim, class BatchedGemmHandleType, class ScalarType,
          class AViewType, class BViewType, class CViewType,
          bool tpl_spec_avail = batched_gemm_tpl_spec_avail<ArgTransA, ArgTransB, ArgBatchSzDim, BatchedGemmHandleType,
                                                            ScalarType, AViewType, BViewType, CViewType>::value,
          bool eti_spec_avail = batched_gemm_eti_spec_avail<ArgTransA, ArgTransB, ArgBatchSzDim, BatchedGemmHandleType,
                                                            ScalarType, AViewType, BViewType, CViewType>::value>
struct BatchedGemmSpec {
  static int run(BatchedGemmHandleType *const handle, const ScalarType alpha, const AViewType &A, const BViewType &B,
                 const ScalarType beta, const CViewType &C)
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
  {
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
#if KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
    printf(
        "KokkosBatched::BatchedGemm<> ETI specialization for < %s, %s, %s, "
        "%s, %s, %s, %s, %s >\n",
        typeid(ArgTransA).name(), typeid(ArgTransB).name(), typeid(ArgBatchSzDim).name(),
        typeid(BatchedGemmHandleType).name(), typeid(ScalarType).name(), typeid(AViewType).name(),
        typeid(BViewType).name(), typeid(CViewType).name());
#else
    printf(
        "KokkosBatched::BatchedGemm<> non-ETI specialization for < %s, %s, "
        "%s, %s, %s, %s, %s, %s >\n",
        typeid(ArgTransA).name(), typeid(ArgTransB).name(), typeid(ArgBatchSzDim).name(),
        typeid(BatchedGemmHandleType).name(), typeid(ScalarType).name(), typeid(AViewType).name(),
        typeid(BViewType).name(), typeid(CViewType).name());
#endif  // KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#endif  // KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    return KokkosBatched::Impl::BatchedGemmImpl<ArgTransA, ArgTransB, ArgBatchSzDim, BatchedGemmHandleType, ScalarType,
                                                AViewType, BViewType, CViewType>(handle, alpha, A, B, beta, C);
  }
#else
      ;
#endif  // !defined(KOKKOSKERNELS_ETI_ONLY) ||
        // KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
};
}  // namespace Impl
}  // namespace KokkosBatched

// ETI instantiation macros, consumed by *.cpp.in files
#define KOKKOSBATCHED_GEMM_ETI_SPEC_INST_INNER(ARG_TRANS_A, ARG_TRANS_B, ARG_BATCH_LAYOUT, SCALAR, LAYOUT, EXEC_SPACE, \
                                               MEM_SPACE)                                                              \
  template struct BatchedGemmSpec<ARG_TRANS_A, ARG_TRANS_B, ARG_BATCH_LAYOUT, BatchedGemmHandle, SCALAR,               \
                                  Kokkos::View<SCALAR ***, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,              \
                                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                               \
                                  Kokkos::View<SCALAR ***, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,              \
                                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                               \
                                  Kokkos::View<SCALAR ***, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,              \
                                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                               \
                                  false, true>;

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
#define KOKKOSBATCHED_GEMM_BLL_ETI_SPEC_INST_INNER(ARG_TRANS_A, ARG_TRANS_B, ARG_BATCH_LAYOUT, SCALAR, LAYOUT,    \
                                                   EXEC_SPACE, MEM_SPACE)                                         \
  KOKKOSBATCHED_GEMM_ETI_SPEC_INST_INNER(ARG_TRANS_A, ARG_TRANS_B, ARG_BATCH_LAYOUT, SCALAR, Kokkos::LayoutRight, \
                                         EXEC_SPACE, MEM_SPACE)
#else
#define KOKKOSBATCHED_GEMM_BLL_ETI_SPEC_INST_INNER(ARG_TRANS_A, ARG_TRANS_B, ARG_BATCH_LAYOUT, SCALAR, LAYOUT, \
                                                   EXEC_SPACE, MEM_SPACE)
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
#define KOKKOSBATCHED_GEMM_BLR_ETI_SPEC_INST_INNER(ARG_TRANS_A, ARG_TRANS_B, ARG_BATCH_LAYOUT, SCALAR, LAYOUT,   \
                                                   EXEC_SPACE, MEM_SPACE)                                        \
  KOKKOSBATCHED_GEMM_ETI_SPEC_INST_INNER(ARG_TRANS_A, ARG_TRANS_B, ARG_BATCH_LAYOUT, SCALAR, Kokkos::LayoutLeft, \
                                         EXEC_SPACE, MEM_SPACE)
#else
#define KOKKOSBATCHED_GEMM_BLR_ETI_SPEC_INST_INNER(ARG_TRANS_A, ARG_TRANS_B, ARG_BATCH_LAYOUT, SCALAR, LAYOUT, \
                                                   EXEC_SPACE, MEM_SPACE)
#endif

///////////////// BatchLayout::Left Permutations /////////////////
#define KOKKOSBATCHED_GEMM_NT_NT_BLL_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                       \
  KOKKOSBATCHED_GEMM_BLL_ETI_SPEC_INST_INNER(Trans::NoTranspose, Trans::NoTranspose, BatchLayout::Left, SCALAR, \
                                             LAYOUT, EXEC_SPACE, MEM_SPACE)

#define KOKKOSBATCHED_GEMM_NT_T_BLL_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                              \
  KOKKOSBATCHED_GEMM_BLL_ETI_SPEC_INST_INNER(Trans::NoTranspose, Trans::Transpose, BatchLayout::Left, SCALAR, LAYOUT, \
                                             EXEC_SPACE, MEM_SPACE)

#define KOKKOSBATCHED_GEMM_T_NT_BLL_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                              \
  KOKKOSBATCHED_GEMM_BLL_ETI_SPEC_INST_INNER(Trans::Transpose, Trans::NoTranspose, BatchLayout::Left, SCALAR, LAYOUT, \
                                             EXEC_SPACE, MEM_SPACE)

#define KOKKOSBATCHED_GEMM_T_T_BLL_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                             \
  KOKKOSBATCHED_GEMM_BLL_ETI_SPEC_INST_INNER(Trans::Transpose, Trans::Transpose, BatchLayout::Left, SCALAR, LAYOUT, \
                                             EXEC_SPACE, MEM_SPACE)

///////////////// BatchLayout::Right Permutations /////////////////
#define KOKKOSBATCHED_GEMM_NT_NT_BLR_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                        \
  KOKKOSBATCHED_GEMM_BLR_ETI_SPEC_INST_INNER(Trans::NoTranspose, Trans::NoTranspose, BatchLayout::Right, SCALAR, \
                                             LAYOUT, EXEC_SPACE, MEM_SPACE)

#define KOKKOSBATCHED_GEMM_NT_T_BLR_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                               \
  KOKKOSBATCHED_GEMM_BLR_ETI_SPEC_INST_INNER(Trans::NoTranspose, Trans::Transpose, BatchLayout::Right, SCALAR, LAYOUT, \
                                             EXEC_SPACE, MEM_SPACE)

#define KOKKOSBATCHED_GEMM_T_NT_BLR_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                               \
  KOKKOSBATCHED_GEMM_BLR_ETI_SPEC_INST_INNER(Trans::Transpose, Trans::NoTranspose, BatchLayout::Right, SCALAR, LAYOUT, \
                                             EXEC_SPACE, MEM_SPACE)

#define KOKKOSBATCHED_GEMM_T_T_BLR_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                              \
  KOKKOSBATCHED_GEMM_BLR_ETI_SPEC_INST_INNER(Trans::Transpose, Trans::Transpose, BatchLayout::Right, SCALAR, LAYOUT, \
                                             EXEC_SPACE, MEM_SPACE)
#endif  // __KOKKOSBATCHED_HOSTLEVEL_GEMM_SPEC_HPP__
