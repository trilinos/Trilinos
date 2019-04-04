/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#ifndef KOKKOSBLAS3_GEMM_SPEC_HPP_
#define KOKKOSBLAS3_GEMM_SPEC_HPP_

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_InnerProductSpaceTraits.hpp"

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include<KokkosBlas3_gemm_impl.hpp>
#endif

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template<class AVT, class BVT, class CVT>
struct gemm_eti_spec_avail {
  enum : bool { value = false };
};
}
}


//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::GEMM.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS3_GEMM_ETI_SPEC_AVAIL_LAYOUT( SCALAR, LAYOUTA, LAYOUTB, LAYOUTC, EXEC_SPACE, MEM_SPACE ) \
    template<> \
    struct gemm_eti_spec_avail< \
         Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         Kokkos::View<const SCALAR**, LAYOUTB, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         Kokkos::View<SCALAR**, LAYOUTC, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> > \
         > { enum : bool { value = true }; };

#define KOKKOSBLAS3_GEMM_ETI_SPEC_AVAIL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
    KOKKOSBLAS3_GEMM_ETI_SPEC_AVAIL_LAYOUT( SCALAR, Kokkos::LayoutLeft, Kokkos::LayoutLeft, LAYOUT, EXEC_SPACE, MEM_SPACE) \
    KOKKOSBLAS3_GEMM_ETI_SPEC_AVAIL_LAYOUT( SCALAR, Kokkos::LayoutLeft, Kokkos::LayoutRight, LAYOUT, EXEC_SPACE, MEM_SPACE) \
    KOKKOSBLAS3_GEMM_ETI_SPEC_AVAIL_LAYOUT( SCALAR, Kokkos::LayoutRight, Kokkos::LayoutLeft, LAYOUT, EXEC_SPACE, MEM_SPACE) \
    KOKKOSBLAS3_GEMM_ETI_SPEC_AVAIL_LAYOUT( SCALAR, Kokkos::LayoutRight, Kokkos::LayoutRight, LAYOUT, EXEC_SPACE, MEM_SPACE)

// Include the actual specialization declarations
#include<KokkosBlas3_gemm_tpl_spec_avail.hpp>
#include<generated_specializations_hpp/KokkosBlas3_gemm_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

//
// gemm
//

// Implementation of KokkosBlas::gemm.
template<class AViewType,
         class BViewType,
         class CViewType,
         bool tpl_spec_avail = gemm_tpl_spec_avail<AViewType, BViewType, CViewType>::value,
         bool eti_spec_avail = gemm_eti_spec_avail<AViewType, BViewType, CViewType>::value
         >
struct GEMM {
  static void
  gemm (const char transA[],
        const char transB[],
        typename AViewType::const_value_type& alpha,
        const AViewType& A,
        const BViewType& B,
        typename CViewType::const_value_type& beta,
        const CViewType& C)
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
{
  static_assert (Kokkos::Impl::is_view<AViewType>::value,
                 "AViewType must be a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<BViewType>::value,
                 "BViewType must be a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<CViewType>::value,
                 "CViewType must be a Kokkos::View.");
  static_assert (static_cast<int> (AViewType::rank) == 2,
                 "AViewType must have rank 2.");
  static_assert (static_cast<int> (BViewType::rank) == 2,
                 "BViewType must have rank 2.");
  static_assert (static_cast<int> (CViewType::rank) == 2,
                 "CViewType must have rank 2.");

  Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY?"KokkosBlas::gemm[ETI]":"KokkosBlas::gemm[noETI]");
  // Figure out Scalar Types
  typedef typename AViewType::non_const_value_type ScalarA;
  typedef typename BViewType::non_const_value_type ScalarB;
  typedef typename CViewType::non_const_value_type ScalarC;

  // Define Blocking sizes (this will be used for scratch spaces)
  static constexpr int blockA0 = 24;
  static constexpr int blockB1 = 64;
  static constexpr int blockA1 = (sizeof(ScalarA)*blockA0*16 + sizeof(ScalarB)*16*blockB1 + sizeof(ScalarC)*blockA0*blockB1 < 24000) ? 16 :
                                 (sizeof(ScalarA)*blockA0*8 + sizeof(ScalarB)*8*blockB1 + sizeof(ScalarC)*blockA0*blockB1 < 24000) ? 8 :
                                 (sizeof(ScalarA)*blockA0*4 + sizeof(ScalarB)*4*blockB1 + sizeof(ScalarC)*blockA0*blockB1 < 24000) ? 4 : 16 ;
  static constexpr int vector_length = blockB1/4;

  // Compute scratch space size
  typedef KokkosBlas::Impl::GEMMImpl<typename CViewType::execution_space,AViewType,BViewType,CViewType,blockA0,blockA1,blockB1,0,0> gemm_dummy_type;
  const int scratch_memory_size =
        gemm_dummy_type::ViewTypeAScratch::required_allocation_size() +
        gemm_dummy_type::ViewTypeBScratch::required_allocation_size() +
        gemm_dummy_type::ViewTypeCScratch::required_allocation_size();
  const int scratch_level = scratch_memory_size < 24000 ? 0 : 1;

  // Figure out Team Sizes
  int team_size = 1;
  #if defined(KOKKOS_ENABLE_CUDA)
  if(std::is_same<typename CViewType::execution_space,Kokkos::Cuda>::value)
    team_size = blockA0;
  #endif
  #if defined(KOKKOS_ENABLE_ROCM)
  if(std::is_same<typename CViewType::execution_space,Kokkos::ROCm>::value)
    team_size = blockA0;
  #endif

  // Call the correct kernel
  if((transA[0]=='N' || transA[0]=='n') && (transB[0]=='N' || transB[0]=='n')) {
    KokkosBlas::Impl::GEMMImpl<typename CViewType::execution_space,AViewType,BViewType,CViewType,blockA0,blockA1,blockB1,0,0> gemm(alpha,A,B,beta,C);
    gemm.run(team_size,vector_length,scratch_level);
  }
  if((transA[0]=='T' || transA[0]=='t') && (transB[0]=='N' || transB[0]=='n')) {
    KokkosBlas::Impl::GEMMImpl<typename CViewType::execution_space,AViewType,BViewType,CViewType,blockA0,blockA1,blockB1,1,0> gemm(alpha,A,B,beta,C);
    gemm.run(team_size,vector_length,scratch_level);
  }
  if((transA[0]=='C' || transA[0]=='c') && (transB[0]=='N' || transB[0]=='n')) {
    KokkosBlas::Impl::GEMMImpl<typename CViewType::execution_space,AViewType,BViewType,CViewType,blockA0,blockA1,blockB1,2,0> gemm(alpha,A,B,beta,C);
    gemm.run(team_size,vector_length,scratch_level);
  }
  if((transA[0]=='N' || transA[0]=='n') && (transB[0]=='T' || transB[0]=='t')) {
    KokkosBlas::Impl::GEMMImpl<typename CViewType::execution_space,AViewType,BViewType,CViewType,blockA0,blockA1,blockB1,0,1> gemm(alpha,A,B,beta,C);
    gemm.run(team_size,vector_length,scratch_level);
  }
  if((transA[0]=='T' || transA[0]=='t') && (transB[0]=='T' || transB[0]=='t')) {
    KokkosBlas::Impl::GEMMImpl<typename CViewType::execution_space,AViewType,BViewType,CViewType,blockA0,blockA1,blockB1,1,1> gemm(alpha,A,B,beta,C);
    gemm.run(team_size,vector_length,scratch_level);
  }
  if((transA[0]=='C' || transA[0]=='c') && (transB[0]=='T' || transB[0]=='t')) {
    KokkosBlas::Impl::GEMMImpl<typename CViewType::execution_space,AViewType,BViewType,CViewType,blockA0,blockA1,blockB1,2,1> gemm(alpha,A,B,beta,C);
    gemm.run(team_size,vector_length,scratch_level);
  }
  if((transA[0]=='N' || transA[0]=='n') && (transB[0]=='C' || transB[0]=='c')) {
    KokkosBlas::Impl::GEMMImpl<typename CViewType::execution_space,AViewType,BViewType,CViewType,blockA0,blockA1,blockB1,0,2> gemm(alpha,A,B,beta,C);
    gemm.run(team_size,vector_length,scratch_level);
  }
  if((transA[0]=='T' || transA[0]=='t') && (transB[0]=='C' || transB[0]=='c')) {
    KokkosBlas::Impl::GEMMImpl<typename CViewType::execution_space,AViewType,BViewType,CViewType,blockA0,blockA1,blockB1,1,2> gemm(alpha,A,B,beta,C);
    gemm.run(team_size,vector_length,scratch_level);
  }
  if((transA[0]=='C' || transA[0]=='c') && (transB[0]=='C' || transB[0]=='c')) {
    KokkosBlas::Impl::GEMMImpl<typename CViewType::execution_space,AViewType,BViewType,CViewType,blockA0,blockA1,blockB1,2,2> gemm(alpha,A,B,beta,C);
    gemm.run(team_size,vector_length,scratch_level);
  }
  Kokkos::Profiling::popRegion();
}
#else
;
#endif //!defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY

};




} // namespace Impl
} // namespace KokkosBlas


//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::GEMM.  This is NOT for users!!!
// All the declarations of full specializations go in this header
// file.  We may spread out definitions (see _DEF macro below) across
// one or more .cpp files.
//

#define KOKKOSBLAS3_GEMM_ETI_SPEC_DECL_LAYOUTS( SCALAR, LAYOUTA, LAYOUTB, LAYOUTC, EXEC_SPACE, MEM_SPACE ) \
extern template struct GEMM< \
     Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<const SCALAR**, LAYOUTB, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<SCALAR**, LAYOUTC, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     false, true>;

#define KOKKOSBLAS3_GEMM_ETI_SPEC_INST_LAYOUTS( SCALAR, LAYOUTA, LAYOUTB, LAYOUTC, EXEC_SPACE, MEM_SPACE ) \
template struct GEMM< \
     Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<const SCALAR**, LAYOUTB, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<SCALAR**, LAYOUTC, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >,  \
     false, true>;

#define KOKKOSBLAS3_GEMM_ETI_SPEC_DECL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
    KOKKOSBLAS3_GEMM_ETI_SPEC_DECL_LAYOUTS(SCALAR, Kokkos::LayoutLeft, Kokkos::LayoutLeft, LAYOUT, EXEC_SPACE, MEM_SPACE) \
    KOKKOSBLAS3_GEMM_ETI_SPEC_DECL_LAYOUTS(SCALAR, Kokkos::LayoutLeft, Kokkos::LayoutRight, LAYOUT, EXEC_SPACE, MEM_SPACE) \
    KOKKOSBLAS3_GEMM_ETI_SPEC_DECL_LAYOUTS(SCALAR, Kokkos::LayoutRight, Kokkos::LayoutLeft, LAYOUT, EXEC_SPACE, MEM_SPACE) \
    KOKKOSBLAS3_GEMM_ETI_SPEC_DECL_LAYOUTS(SCALAR, Kokkos::LayoutRight, Kokkos::LayoutRight, LAYOUT, EXEC_SPACE, MEM_SPACE)

#define KOKKOSBLAS3_GEMM_ETI_SPEC_INST( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
    KOKKOSBLAS3_GEMM_ETI_SPEC_INST_LAYOUTS(SCALAR, Kokkos::LayoutLeft, Kokkos::LayoutLeft, LAYOUT, EXEC_SPACE, MEM_SPACE) \
    KOKKOSBLAS3_GEMM_ETI_SPEC_INST_LAYOUTS(SCALAR, Kokkos::LayoutLeft, Kokkos::LayoutRight, LAYOUT, EXEC_SPACE, MEM_SPACE) \
    KOKKOSBLAS3_GEMM_ETI_SPEC_INST_LAYOUTS(SCALAR, Kokkos::LayoutRight, Kokkos::LayoutLeft, LAYOUT, EXEC_SPACE, MEM_SPACE) \
    KOKKOSBLAS3_GEMM_ETI_SPEC_INST_LAYOUTS(SCALAR, Kokkos::LayoutRight, Kokkos::LayoutRight, LAYOUT, EXEC_SPACE, MEM_SPACE)

#include<KokkosBlas3_gemm_tpl_spec_decl.hpp>
#include<generated_specializations_hpp/KokkosBlas3_gemm_eti_spec_decl.hpp>

#endif // KOKKOSBLAS3_GEMM_SPEC_HPP_
