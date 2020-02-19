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
#ifndef KOKKOSBLAS3_TRSM_SPEC_HPP_
#define KOKKOSBLAS3_TRSM_SPEC_HPP_

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_InnerProductSpaceTraits.hpp"
#include <sstream>

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include<KokkosBlas3_trsm_impl.hpp>
#endif

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template<class AVT, class BVT>
struct trsm_eti_spec_avail {
  enum : bool { value = false };
};
}
}

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::TRSM.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS3_TRSM_ETI_SPEC_AVAIL_LAYOUT( SCALAR, LAYOUTA, LAYOUTB, EXEC_SPACE, MEM_SPACE ) \
    template<> \
    struct trsm_eti_spec_avail< \
         Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         Kokkos::View<SCALAR**, LAYOUTB, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> > \
         > { enum : bool { value = true }; };

#define KOKKOSBLAS3_TRSM_ETI_SPEC_AVAIL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
    KOKKOSBLAS3_TRSM_ETI_SPEC_AVAIL_LAYOUT( SCALAR, LAYOUT, LAYOUT, EXEC_SPACE, MEM_SPACE)

// Include the actual specialization declarations
#include<KokkosBlas3_trsm_tpl_spec_avail.hpp>
#include<generated_specializations_hpp/KokkosBlas3_trsm_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

//
// trsm
//

//Unification layer
template<class AViewType,
         class BViewType,
         bool tpl_spec_avail = trsm_tpl_spec_avail<AViewType, BViewType>::value,
         bool eti_spec_avail = trsm_eti_spec_avail<AViewType, BViewType>::value
        >
struct TRSM{
  static void
  trsm (const char side[],
        const char uplo[],
        const char trans[],
        const char diag[],
        typename BViewType::const_value_type& alpha,
        const AViewType& A,
        const BViewType& B);
};

// Implementation of KokkosBlas::trsm.
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
template<class AViewType,
         class BViewType>
struct TRSM<AViewType, BViewType, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void
  trsm (const char side[],
        const char uplo[],
        const char trans[],
        const char diag[],
        typename BViewType::const_value_type& alpha,
        const AViewType& A,
        const BViewType& B)
  {
    static_assert (Kokkos::Impl::is_view<AViewType>::value,
                   "AViewType must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<BViewType>::value,
                   "BViewType must be a Kokkos::View.");
    static_assert (static_cast<int> (AViewType::rank) == 2,
                   "AViewType must have rank 2.");
    static_assert (static_cast<int> (BViewType::rank) == 2,
                   "BViewType must have rank 2.");

    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY?"KokkosBlas::trsm[ETI]":"KokkosBlas::trsm[noETI]");

    typename AViewType::HostMirror h_A  = Kokkos::create_mirror_view(A);
    typename BViewType::HostMirror h_B  = Kokkos::create_mirror_view(B);

    Kokkos::deep_copy(h_A, A);
    Kokkos::deep_copy(h_B, B);

    SerialTrsm_Invoke<typename AViewType::HostMirror, typename BViewType::HostMirror> (side, uplo, trans, diag, alpha, h_A, h_B);

    Kokkos::deep_copy(B, h_B);

    Kokkos::Profiling::popRegion();
  }
};
#endif //!defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY

} // namespace Impl
} // namespace KokkosBlas


//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::TRSM.  This is NOT for users!!!
// All the declarations of full specializations go in this header
// file. We may spread out definitions (see _DEF macro below) across
// one or more .cpp files.
//

#define KOKKOSBLAS3_TRSM_ETI_SPEC_DECL_LAYOUTS( SCALAR, LAYOUTA, LAYOUTB, EXEC_SPACE, MEM_SPACE ) \
extern template struct TRSM< \
     Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<SCALAR**, LAYOUTB, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     false, true>;

#define KOKKOSBLAS3_TRSM_ETI_SPEC_INST_LAYOUTS( SCALAR, LAYOUTA, LAYOUTB, EXEC_SPACE, MEM_SPACE ) \
template struct TRSM< \
     Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<SCALAR**, LAYOUTB, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     false, true>;

#define KOKKOSBLAS3_TRSM_ETI_SPEC_DECL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
    KOKKOSBLAS3_TRSM_ETI_SPEC_DECL_LAYOUTS(SCALAR, LAYOUT, LAYOUT, EXEC_SPACE, MEM_SPACE)

#define KOKKOSBLAS3_TRSM_ETI_SPEC_INST( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
    KOKKOSBLAS3_TRSM_ETI_SPEC_INST_LAYOUTS(SCALAR, LAYOUT, LAYOUT, EXEC_SPACE, MEM_SPACE)

#include<KokkosBlas3_trsm_tpl_spec_decl.hpp>
#include<generated_specializations_hpp/KokkosBlas3_trsm_eti_spec_decl.hpp>

#endif // KOKKOSBLAS3_TRSM_SPEC_HPP_
