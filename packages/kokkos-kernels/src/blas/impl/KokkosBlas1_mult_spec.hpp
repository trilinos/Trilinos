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
#ifndef KOKKOSBLAS1_MULT_SPEC_HPP_
#define KOKKOSBLAS1_MULT_SPEC_HPP_

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_InnerProductSpaceTraits.hpp"

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include<KokkosBlas1_mult_impl.hpp>
#endif

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template<class YMV, class AV, class XMV, int rank = XMV::rank>
struct mult_eti_spec_avail {
  enum : bool { value = false };
};
}
}

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Mult for rank == 1.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_MULT_ETI_SPEC_AVAIL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
    template<> \
    struct mult_eti_spec_avail< \
         Kokkos::View<SCALAR*, \
                      LAYOUT, \
                      Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         Kokkos::View<const SCALAR*, \
                      LAYOUT, \
                      Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         Kokkos::View<const SCALAR*, \
                      LAYOUT, \
                      Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         1> { enum : bool { value = true }; };

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Mult for rank == 2.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_MULT_MV_ETI_SPEC_AVAIL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
    template<> \
    struct mult_eti_spec_avail< \
         Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         Kokkos::View<const SCALAR*, \
                      std::conditional<std::is_same<LAYOUT,Kokkos::LayoutRight>::value, \
                               Kokkos::LayoutLeft, LAYOUT>::type, \
                      Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         2> { enum : bool { value = true }; };


// Include the actual specialization declarations
#include<KokkosBlas1_mult_tpl_spec_avail.hpp>
#include<generated_specializations_hpp/KokkosBlas1_mult_eti_spec_avail.hpp>
#include<generated_specializations_hpp/KokkosBlas1_mult_mv_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

//
// mult
//
/// \brief Implementation of entry-wise multiply of multivectors or
///   single vectors (depending on the rank template parameter).
///
/// Compute
///
/// Y(i,j) = alpha*A(i,j)*X(i,j) + gamma*Y(i,j)
///
/// with special cases for alpha, or gamma = 0.
template<class YMV, class AV, class XMV, int rank = XMV::rank,
    bool tpl_spec_avail = mult_tpl_spec_avail<YMV,AV,XMV>::value,
    bool eti_spec_avail = mult_eti_spec_avail<YMV,AV,XMV>::value>
struct Mult {
  static void
    mult (const typename YMV::non_const_value_type& gamma, const YMV& Y,
          const typename XMV::non_const_value_type& alpha, const AV& A, const XMV& X);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
// Partial specialization for YMV, AV, and XMV rank-2 Views.
template<class YMV, class AV, class XMV>
struct Mult<YMV, AV, XMV, 2, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY>
{
  typedef typename YMV::size_type size_type;
  typedef typename YMV::non_const_value_type YMV_scalar;
  typedef typename XMV::non_const_value_type XMV_scalar;

  static void
  mult (const YMV_scalar& gamma, const YMV& Y,
        const XMV_scalar& alpha, const AV& A, const XMV& X)
  {
    static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::Impl::"
                   "Mult<rank 2>::mult: Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<AV>::value, "KokkosBlas::Impl::"
                   "Mult<rank 2>::mult: A is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "Mult<rank 2>::mult: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename YMV::value_type,
                     typename YMV::non_const_value_type>::value,
                   "KokkosBlas::Impl::Mult<rank 2>::mult: Y is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    // Casting to int avoids compiler warnings about comparing
    // different kinds of enum values.
    static_assert ((int) XMV::rank == (int) YMV::rank &&
                   (int) XMV::rank == 2,
                   "KokkosBlas::Impl::Mult<rank 2>::mult: "
                   "X, and Y must have the rank 2.");
    static_assert (AV::rank == 1, "KokkosBlas::Impl::Mult<rank 2>::mult: "
                   "AV must have rank 1.");
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY?"KokkosBlas::mult[ETI]":"KokkosBlas::mult[noETI]");

    #ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::mult<> ETI specialization for < %s , %s , %s >\n",typeid(YMV).name(),typeid(AV).name(),typeid(XMV).name());
    else {
      printf("KokkosBlas1::mult<> non-ETI specialization for < %s , %s , %s >\n",typeid(YMV).name(),typeid(AV).name(),typeid(XMV).name());
    }
    #endif

    const size_type numRows = X.extent(0);
    const size_type numCols = X.extent(1);

    if (numRows < static_cast<int> (INT_MAX) &&
        numRows * numCols < static_cast<int> (INT_MAX)) {
      MV_Mult_Generic<YMV, AV, XMV, int> (gamma, Y, alpha, A, X);
    }
    else {
      MV_Mult_Generic<YMV, AV, XMV, int64_t> (gamma, Y, alpha, A, X);
    }
    Kokkos::Profiling::popRegion();
  }
};

// Partial specialization for YV, AV, and XV rank-1 Views.
template<class YV, class AV, class XV>
struct Mult<YV, AV, XV, 1, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY>
{
  typedef typename YV::size_type size_type;
  typedef typename YV::non_const_value_type YV_scalar;
  typedef typename XV::non_const_value_type XV_scalar;

  static void
  mult (const YV_scalar& gamma, const YV& Y,
        const XV_scalar& alpha, const AV& A, const XV& X)
  {
    // YV, AV, and XV must be Kokkos::View specializations.
    static_assert (Kokkos::Impl::is_view<YV>::value, "KokkosBlas::Impl::"
                   "Mult<rank 1>::mult: Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<AV>::value, "KokkosBlas::Impl::"
                   "Mult<rank 1>::mult: A is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XV>::value, "KokkosBlas::Impl::"
                   "Mult<rank 1>::mult: X is not a Kokkos::View.");
    // XV must be nonconst (else it can't be an output argument).
    static_assert (Kokkos::Impl::is_same<typename YV::value_type,
                     typename YV::non_const_value_type>::value,
                   "KokkosBlas::Impl::Mult<rank 1>::mult: Y is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert ((int) XV::rank == (int) YV::rank && (int) AV::rank == 1,
                   "KokkosBlas::Impl::Mult<rank 1>::mult: "
                   "X, Y, and Z must have rank 1.");
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY?"KokkosBlas::mult[ETI]":"KokkosBlas::mult[noETI]");
    #ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::mult<> ETI specialization for < %s , %s , %s >\n",typeid(YV).name(),typeid(AV).name(),typeid(XV).name());
    else {
      printf("KokkosBlas1::mult<> non-ETI specialization for < %s , %s , %s >\n",typeid(YV).name(),typeid(AV).name(),typeid(XV).name());
    }
    #endif
 
    const size_type numRows = Y.extent(0);
    if (numRows < static_cast<int> (INT_MAX)) {
      V_Mult_Generic<YV, AV, XV, int> (gamma, Y, alpha, A, X);
    }
    else {
      V_Mult_Generic<YV, AV, XV, int64_t> (gamma, Y, alpha, A, X);
    }
    Kokkos::Profiling::popRegion();
  }
};
#endif //!defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY



} // namespace Impl
} // namespace KokkosBlas

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Mult for rank == 1.  This is NOT for users!!!
// All the declarations of full specializations go in this header
// file.  We may spread out definitions (see _INST macro below) across
// one or more .cpp files.
//

#define KOKKOSBLAS1_MULT_ETI_SPEC_DECL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
extern template struct Mult< \
         Kokkos::View<SCALAR*, \
                      LAYOUT, \
                      Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         Kokkos::View<const SCALAR*, \
                      LAYOUT, \
                      Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         Kokkos::View<const SCALAR*, \
                      std::conditional<std::is_same<LAYOUT,Kokkos::LayoutRight>::value, \
                               Kokkos::LayoutLeft, LAYOUT>::type, \
                      Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        1, false, true>;

#define KOKKOSBLAS1_MULT_ETI_SPEC_INST( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
template struct Mult< \
         Kokkos::View<SCALAR*, \
                      LAYOUT, \
                      Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         Kokkos::View<const SCALAR*, \
                      LAYOUT, \
                      Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         Kokkos::View<const SCALAR*, \
                      LAYOUT, \
                      Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        1, false, true>;

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Mult for rank == 2.  This is NOT for users!!!
// All the declarations of full specializations go in this header
// file.  We may spread out definitions (see _DEF macro below) across
// one or more .cpp files.
//

#define KOKKOSBLAS1_MULT_MV_ETI_SPEC_DECL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
extern template struct Mult< \
     Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<const SCALAR*, \
                  std::conditional<std::is_same<LAYOUT,Kokkos::LayoutRight>::value, \
                                            Kokkos::LayoutLeft, LAYOUT>::type, \
                  Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     2, false, true>;

#define KOKKOSBLAS1_MULT_MV_ETI_SPEC_INST( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
template struct Mult< \
     Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<const SCALAR*, \
                  std::conditional<std::is_same<LAYOUT,Kokkos::LayoutRight>::value, \
                                            Kokkos::LayoutLeft, LAYOUT>::type, \
                  Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     2, false, true>;


#include<KokkosBlas1_mult_tpl_spec_decl.hpp>
#include<generated_specializations_hpp/KokkosBlas1_mult_eti_spec_decl.hpp>
#include<generated_specializations_hpp/KokkosBlas1_mult_mv_eti_spec_decl.hpp>

#endif // KOKKOSBLAS1_MULT_SPEC_HPP_
