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
#ifndef KOKKOSBLAS3_GEMM_SPEC_HPP_
#define KOKKOSBLAS3_GEMM_SPEC_HPP_

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_InnerProductSpaceTraits.hpp"

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include "KokkosBlas3_gemm_impl.hpp"
#include "KokkosBlas3_gemm_dotbased_impl.hpp"
#include "KokkosKernels_ExecSpaceUtils.hpp"
#endif

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class execution_space, class AVT, class BVT, class CVT>
struct gemm_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::GEMM.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS3_GEMM_ETI_SPEC_AVAIL_LAYOUT(SCALAR, LAYOUTA, LAYOUTB, LAYOUTC, EXEC_SPACE, MEM_SPACE)  \
  template <>                                                                                             \
  struct gemm_eti_spec_avail<EXEC_SPACE,                                                                  \
                             Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                      \
                             Kokkos::View<const SCALAR**, LAYOUTB, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                      \
                             Kokkos::View<SCALAR**, LAYOUTC, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,       \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                   \
    enum : bool { value = true };                                                                         \
  };

#define KOKKOSBLAS3_GEMM_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                 \
  KOKKOSBLAS3_GEMM_ETI_SPEC_AVAIL_LAYOUT(SCALAR, Kokkos::LayoutLeft, Kokkos::LayoutLeft, LAYOUT, EXEC_SPACE,   \
                                         MEM_SPACE)                                                            \
  KOKKOSBLAS3_GEMM_ETI_SPEC_AVAIL_LAYOUT(SCALAR, Kokkos::LayoutLeft, Kokkos::LayoutRight, LAYOUT, EXEC_SPACE,  \
                                         MEM_SPACE)                                                            \
  KOKKOSBLAS3_GEMM_ETI_SPEC_AVAIL_LAYOUT(SCALAR, Kokkos::LayoutRight, Kokkos::LayoutLeft, LAYOUT, EXEC_SPACE,  \
                                         MEM_SPACE)                                                            \
  KOKKOSBLAS3_GEMM_ETI_SPEC_AVAIL_LAYOUT(SCALAR, Kokkos::LayoutRight, Kokkos::LayoutRight, LAYOUT, EXEC_SPACE, \
                                         MEM_SPACE)

// Include the actual specialization declarations
#include <KokkosBlas3_gemm_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBlas3_gemm_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

//
// gemm
//

// Implementation of KokkosBlas::gemm.
template <class execution_space, class AViewType, class BViewType, class CViewType,
          bool tpl_spec_avail = gemm_tpl_spec_avail<execution_space, AViewType, BViewType, CViewType>::value,
          bool eti_spec_avail = gemm_eti_spec_avail<execution_space, AViewType, BViewType, CViewType>::value>
struct GEMM {
  static void gemm(const execution_space& space, const char transA[], const char transB[],
                   typename AViewType::const_value_type& alpha, const AViewType& A, const BViewType& B,
                   typename CViewType::const_value_type& beta, const CViewType& C)
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
  {
    static_assert(Kokkos::is_view<AViewType>::value, "AViewType must be a Kokkos::View.");
    static_assert(Kokkos::is_view<BViewType>::value, "BViewType must be a Kokkos::View.");
    static_assert(Kokkos::is_view<CViewType>::value, "CViewType must be a Kokkos::View.");
    static_assert(static_cast<int>(AViewType::rank) == 2, "AViewType must have rank 2.");
    static_assert(static_cast<int>(BViewType::rank) == 2, "BViewType must have rank 2.");
    static_assert(static_cast<int>(CViewType::rank) == 2, "CViewType must have rank 2.");

    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY ? "KokkosBlas::gemm[ETI]"
                                                                     : "KokkosBlas::gemm[noETI]");
    // Figure out Scalar Types
    typedef typename AViewType::non_const_value_type ScalarA;
    typedef typename BViewType::non_const_value_type ScalarB;
    typedef typename CViewType::non_const_value_type ScalarC;

    // Figure out whether to use DotBased implementation
    const int M = static_cast<int>(C.extent(0));
    const int N = static_cast<int>(C.extent(1));

    const bool is_device_space = KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>();
    const bool A_is_lr         = std::is_same<Kokkos::LayoutRight, typename AViewType::array_layout>::value;
    const bool A_is_tr         = ((transA[0] == 'T') || (transA[0] == 't') || (transA[0] == 'C') || (transA[0] == 'c'));
    const bool B_is_tr         = ((transB[0] == 'T') || (transB[0] == 't') || (transB[0] == 'C') || (transB[0] == 'c'));

    // NOTE: these thresholds were copied from TPL CUBLAS, and may need to be
    // retuned
    constexpr int numDotsLayoutLeftThreshold  = 1600;
    constexpr int numDotsLayoutRightThreshold = 100;
    if (((!A_is_lr && A_is_tr && !B_is_tr && M * N < numDotsLayoutLeftThreshold) ||
         (A_is_lr && A_is_tr && !B_is_tr && M * N < numDotsLayoutRightThreshold)) &&
        is_device_space) {
      // call dot-based GEMM, only for C := beta * C + alpha * A^T * B, on
      // device
      bool A_is_conj = ((transA[0] == 'C') || (transA[0] == 'c'));
      DotBasedGEMM<execution_space, AViewType, BViewType, CViewType> dotBasedGemm(alpha, A, B, beta, C);
      dotBasedGemm.run(space, A_is_conj);

    } else {
      // Define Blocking sizes (this will be used for scratch spaces)
      static constexpr int blockA0 = 24;
      static constexpr int blockB1 = 64;
      static constexpr int blockA1 =
          (sizeof(ScalarA) * blockA0 * 16 + sizeof(ScalarB) * 16 * blockB1 + sizeof(ScalarC) * blockA0 * blockB1 <
           24000)
              ? 16
          : (sizeof(ScalarA) * blockA0 * 8 + sizeof(ScalarB) * 8 * blockB1 + sizeof(ScalarC) * blockA0 * blockB1 <
             24000)
              ? 8
          : (sizeof(ScalarA) * blockA0 * 4 + sizeof(ScalarB) * 4 * blockB1 + sizeof(ScalarC) * blockA0 * blockB1 <
             24000)
              ? 4
              : 16;
      int vector_length     = blockB1 / 4;
      int max_vector_length = KokkosKernels::Impl::kk_get_max_vector_size<execution_space>();
      if (vector_length > max_vector_length) vector_length = max_vector_length;

      // Compute scratch space size
      typedef KokkosBlas::Impl::GEMMImpl<execution_space, AViewType, BViewType, CViewType, blockA0, blockA1, blockB1, 0,
                                         0>
          gemm_dummy_type;
      const int scratch_memory_size = gemm_dummy_type::ViewTypeAScratch::required_allocation_size() +
                                      gemm_dummy_type::ViewTypeBScratch::required_allocation_size() +
                                      gemm_dummy_type::ViewTypeCScratch::required_allocation_size();
      const int scratch_level = scratch_memory_size < 24000 ? 0 : 1;

      // Figure out Team Sizes
      int team_size = 1;
#if defined(KOKKOS_ENABLE_CUDA)
      if (std::is_same<execution_space, Kokkos::Cuda>::value) team_size = blockA0;
#endif
#if defined(KOKKOS_ENABLE_HIP)
      if (std::is_same<execution_space, Kokkos::HIP>::value) team_size = blockA0;
#endif
#if defined(KOKKOS_ENABLE_ROCM)
      if (std::is_same<execution_space, Kokkos::ROCm>::value) team_size = blockA0;
#endif
#if defined(KOKKOS_ENABLE_SYCL)
      if (std::is_same<execution_space, Kokkos::Experimental::SYCL>::value) team_size = blockA0;
#endif

      // Call the correct kernel
      if ((transA[0] == 'N' || transA[0] == 'n') && (transB[0] == 'N' || transB[0] == 'n')) {
        KokkosBlas::Impl::GEMMImpl<execution_space, AViewType, BViewType, CViewType, blockA0, blockA1, blockB1, 0, 0>
            gemm(alpha, A, B, beta, C);
        gemm.run(space, team_size, vector_length, scratch_level);
      }
      if ((transA[0] == 'T' || transA[0] == 't') && (transB[0] == 'N' || transB[0] == 'n')) {
        KokkosBlas::Impl::GEMMImpl<execution_space, AViewType, BViewType, CViewType, blockA0, blockA1, blockB1, 1, 0>
            gemm(alpha, A, B, beta, C);
        gemm.run(space, team_size, vector_length, scratch_level);
      }
      if ((transA[0] == 'C' || transA[0] == 'c') && (transB[0] == 'N' || transB[0] == 'n')) {
        KokkosBlas::Impl::GEMMImpl<execution_space, AViewType, BViewType, CViewType, blockA0, blockA1, blockB1, 2, 0>
            gemm(alpha, A, B, beta, C);
        gemm.run(space, team_size, vector_length, scratch_level);
      }
      if ((transA[0] == 'N' || transA[0] == 'n') && (transB[0] == 'T' || transB[0] == 't')) {
        KokkosBlas::Impl::GEMMImpl<execution_space, AViewType, BViewType, CViewType, blockA0, blockA1, blockB1, 0, 1>
            gemm(alpha, A, B, beta, C);
        gemm.run(space, team_size, vector_length, scratch_level);
      }
      if ((transA[0] == 'T' || transA[0] == 't') && (transB[0] == 'T' || transB[0] == 't')) {
        KokkosBlas::Impl::GEMMImpl<execution_space, AViewType, BViewType, CViewType, blockA0, blockA1, blockB1, 1, 1>
            gemm(alpha, A, B, beta, C);
        gemm.run(space, team_size, vector_length, scratch_level);
      }
      if ((transA[0] == 'C' || transA[0] == 'c') && (transB[0] == 'T' || transB[0] == 't')) {
        KokkosBlas::Impl::GEMMImpl<execution_space, AViewType, BViewType, CViewType, blockA0, blockA1, blockB1, 2, 1>
            gemm(alpha, A, B, beta, C);
        gemm.run(space, team_size, vector_length, scratch_level);
      }
      if ((transA[0] == 'N' || transA[0] == 'n') && (transB[0] == 'C' || transB[0] == 'c')) {
        KokkosBlas::Impl::GEMMImpl<execution_space, AViewType, BViewType, CViewType, blockA0, blockA1, blockB1, 0, 2>
            gemm(alpha, A, B, beta, C);
        gemm.run(space, team_size, vector_length, scratch_level);
      }
      if ((transA[0] == 'T' || transA[0] == 't') && (transB[0] == 'C' || transB[0] == 'c')) {
        KokkosBlas::Impl::GEMMImpl<execution_space, AViewType, BViewType, CViewType, blockA0, blockA1, blockB1, 1, 2>
            gemm(alpha, A, B, beta, C);
        gemm.run(space, team_size, vector_length, scratch_level);
      }
      if ((transA[0] == 'C' || transA[0] == 'c') && (transB[0] == 'C' || transB[0] == 'c')) {
        KokkosBlas::Impl::GEMMImpl<execution_space, AViewType, BViewType, CViewType, blockA0, blockA1, blockB1, 2, 2>
            gemm(alpha, A, B, beta, C);
        gemm.run(space, team_size, vector_length, scratch_level);
      }
    }
    Kokkos::Profiling::popRegion();
  }
#else
      ;
#endif  //! defined(KOKKOSKERNELS_ETI_ONLY) ||
        //! KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
};

}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::GEMM.  This is NOT for users!!!
// All the declarations of full specializations go in this header
// file.  We may spread out definitions (see _DEF macro below) across
// one or more .cpp files.
//

#define KOKKOSBLAS3_GEMM_ETI_SPEC_DECL_LAYOUTS(SCALAR, LAYOUTA, LAYOUTB, LAYOUTC, EXEC_SPACE, MEM_SPACE)   \
  extern template struct GEMM<EXEC_SPACE,                                                                  \
                              Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                      \
                              Kokkos::View<const SCALAR**, LAYOUTB, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                      \
                              Kokkos::View<SCALAR**, LAYOUTC, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,       \
                                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                      \
                              false, true>;

#define KOKKOSBLAS3_GEMM_ETI_SPEC_INST_LAYOUTS(SCALAR, LAYOUTA, LAYOUTB, LAYOUTC, EXEC_SPACE, MEM_SPACE) \
  template struct GEMM<EXEC_SPACE,                                                                       \
                       Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,      \
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                           \
                       Kokkos::View<const SCALAR**, LAYOUTB, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,      \
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                           \
                       Kokkos::View<SCALAR**, LAYOUTC, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,            \
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                           \
                       false, true>;

#define KOKKOSBLAS3_GEMM_ETI_SPEC_DECL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                  \
  KOKKOSBLAS3_GEMM_ETI_SPEC_DECL_LAYOUTS(SCALAR, Kokkos::LayoutLeft, Kokkos::LayoutLeft, LAYOUT, EXEC_SPACE,   \
                                         MEM_SPACE)                                                            \
  KOKKOSBLAS3_GEMM_ETI_SPEC_DECL_LAYOUTS(SCALAR, Kokkos::LayoutLeft, Kokkos::LayoutRight, LAYOUT, EXEC_SPACE,  \
                                         MEM_SPACE)                                                            \
  KOKKOSBLAS3_GEMM_ETI_SPEC_DECL_LAYOUTS(SCALAR, Kokkos::LayoutRight, Kokkos::LayoutLeft, LAYOUT, EXEC_SPACE,  \
                                         MEM_SPACE)                                                            \
  KOKKOSBLAS3_GEMM_ETI_SPEC_DECL_LAYOUTS(SCALAR, Kokkos::LayoutRight, Kokkos::LayoutRight, LAYOUT, EXEC_SPACE, \
                                         MEM_SPACE)

#define KOKKOSBLAS3_GEMM_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                  \
  KOKKOSBLAS3_GEMM_ETI_SPEC_INST_LAYOUTS(SCALAR, Kokkos::LayoutLeft, Kokkos::LayoutLeft, LAYOUT, EXEC_SPACE,   \
                                         MEM_SPACE)                                                            \
  KOKKOSBLAS3_GEMM_ETI_SPEC_INST_LAYOUTS(SCALAR, Kokkos::LayoutLeft, Kokkos::LayoutRight, LAYOUT, EXEC_SPACE,  \
                                         MEM_SPACE)                                                            \
  KOKKOSBLAS3_GEMM_ETI_SPEC_INST_LAYOUTS(SCALAR, Kokkos::LayoutRight, Kokkos::LayoutLeft, LAYOUT, EXEC_SPACE,  \
                                         MEM_SPACE)                                                            \
  KOKKOSBLAS3_GEMM_ETI_SPEC_INST_LAYOUTS(SCALAR, Kokkos::LayoutRight, Kokkos::LayoutRight, LAYOUT, EXEC_SPACE, \
                                         MEM_SPACE)

#include <KokkosBlas3_gemm_tpl_spec_decl.hpp>

#endif  // KOKKOSBLAS3_GEMM_SPEC_HPP_
