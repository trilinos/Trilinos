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
#ifndef KOKKOSBLAS2_GEMV_SPEC_HPP_
#define KOKKOSBLAS2_GEMV_SPEC_HPP_

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_InnerProductSpaceTraits.hpp"

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include<KokkosBlas2_gemv_impl.hpp>
#endif

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template<class XMV, class YMV, class ZMV>
struct gemv_eti_spec_avail {
  enum : bool { value = false };
};
}
}


//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::GEMV.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS2_GEMV_ETI_SPEC_AVAIL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
    template<> \
    struct gemv_eti_spec_avail< \
         Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         Kokkos::View<const SCALAR*, \
                      typename std::conditional<std::is_same<LAYOUT,Kokkos::LayoutRight>::value, \
                                                             Kokkos::LayoutLeft, LAYOUT>::type, \
                      Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         Kokkos::View<SCALAR*, \
                      typename std::conditional<std::is_same<LAYOUT,Kokkos::LayoutRight>::value, \
                                                             Kokkos::LayoutLeft, LAYOUT>::type, \
                      Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> > \
         > { enum : bool { value = true }; };


// Include the actual specialization declarations
#include<KokkosBlas2_gemv_tpl_spec_avail.hpp>
#include<generated_specializations_hpp/KokkosBlas2_gemv_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

//
// gemv
//

// Implementation of KokkosBlas::gemv.
template<class AViewType,
         class XViewType,
         class YViewType,
         bool tpl_spec_avail = gemv_tpl_spec_avail<AViewType, XViewType, YViewType>::value,
         bool eti_spec_avail = gemv_eti_spec_avail<AViewType, XViewType, YViewType>::value
         >
struct GEMV {
  static void
  gemv (const char trans[],
        typename AViewType::const_value_type& alpha,
        const AViewType& A,
        const XViewType& x,
        typename YViewType::const_value_type& beta,
        const YViewType& y)
  #if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
  {
    static_assert (Kokkos::Impl::is_view<AViewType>::value,
                   "AViewType must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XViewType>::value,
                   "XViewType must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YViewType>::value,
                   "YViewType must be a Kokkos::View.");
    static_assert (static_cast<int> (AViewType::rank) == 2,
                   "AViewType must have rank 2.");
    static_assert (static_cast<int> (XViewType::rank) == 1,
                   "XViewType must have rank 1.");
    static_assert (static_cast<int> (YViewType::rank) == 1,
                   "YViewType must have rank 1.");

    typedef typename AViewType::size_type size_type;
    const size_type numRows = A.extent(0);
    const size_type numCols = A.extent(1);

    // Prefer int as the index type, but use a larger type if needed.
    if (numRows < static_cast<size_type> (INT_MAX) &&
        numCols < static_cast<size_type> (INT_MAX)) {
      singleLevelGemv<AViewType, XViewType, YViewType,int>
         (trans, alpha, A, x, beta, y);
    }
    else {
      singleLevelGemv<AViewType, XViewType, YViewType, int64_t>
         (trans, alpha, A, x, beta, y);
    }
  }
  #else
  ;
  #endif //!defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
};




} // namespace Impl
} // namespace KokkosBlas


//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::GEMV.  This is NOT for users!!!
// All the declarations of full specializations go in this header
// file.  We may spread out definitions (see _DEF macro below) across
// one or more .cpp files.
//

#define KOKKOSBLAS2_GEMV_ETI_SPEC_DECL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
extern template struct GEMV< \
     Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<const SCALAR*, \
                  typename std::conditional<std::is_same<LAYOUT,Kokkos::LayoutRight>::value, \
                                                         Kokkos::LayoutLeft, LAYOUT>::type, \
                  Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<SCALAR*, \
                  typename std::conditional<std::is_same<LAYOUT,Kokkos::LayoutRight>::value, \
                                                         Kokkos::LayoutLeft, LAYOUT>::type, \
                  Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     false, true>;

#define KOKKOSBLAS2_GEMV_ETI_SPEC_INST( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
template struct GEMV< \
     Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<const SCALAR*, \
                  typename std::conditional<std::is_same<LAYOUT,Kokkos::LayoutRight>::value, \
                                                         Kokkos::LayoutLeft, LAYOUT>::type, \
                  Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<SCALAR*, \
                  typename std::conditional<std::is_same<LAYOUT,Kokkos::LayoutRight>::value, \
                                                         Kokkos::LayoutLeft, LAYOUT>::type, \
                  Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >,  \
     false, true>;


#include<KokkosBlas2_gemv_tpl_spec_decl.hpp>
#include<generated_specializations_hpp/KokkosBlas2_gemv_eti_spec_decl.hpp>

#endif // KOKKOSBLAS1_GEMV_SPEC_HPP_
