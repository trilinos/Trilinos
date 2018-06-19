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
#ifndef KOKKOSBLAS1_UPDATE_SPEC_HPP_
#define KOKKOSBLAS1_UPDATE_SPEC_HPP_

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_InnerProductSpaceTraits.hpp"

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include<KokkosBlas1_update_impl.hpp>
#endif

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template<class XMV, class YMV, class ZMV, int rank = ZMV::rank>
struct update_eti_spec_avail {
  enum : bool { value = false };
};
}
}

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Update for rank == 1.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_UPDATE_ETI_SPEC_AVAIL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
    template<> \
    struct update_eti_spec_avail< \
         Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         1> { enum : bool { value = true }; };

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Update for rank == 2.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_UPDATE_MV_ETI_SPEC_AVAIL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
    template<> \
    struct update_eti_spec_avail< \
         Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         2> { enum : bool { value = true }; };


// Include the actual specialization declarations
#include<KokkosBlas1_update_tpl_spec_avail.hpp>
#include<generated_specializations_hpp/KokkosBlas1_update_eti_spec_avail.hpp>
#include<generated_specializations_hpp/KokkosBlas1_update_mv_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

//
// update
//

/// \brief Implementation of KokkosBlas::update for single vectors and
///   multivectors.
///
/// Compute
///
/// Z(i,j) = alpha*X(i,j) + beta*Y(i,j) + gamma*Z(i,j),
///
/// with special cases for alpha, beta, or gamma = 0.
template<class XMV, class YMV, class ZMV, int rank = ZMV::rank,
    bool tpl_spec_avail = update_tpl_spec_avail<XMV,YMV,ZMV>::value,
    bool eti_spec_avail = update_eti_spec_avail<XMV,YMV,ZMV>::value>
struct Update {
  static void
    update (const typename XMV::non_const_value_type& alpha, const XMV& X,
            const typename YMV::non_const_value_type& beta, const YMV& Y,
            const typename ZMV::non_const_value_type& gamma, const ZMV& Z);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
// Partial specialization for XMV, YMV, and ZMV rank-2 Views.
template<class XMV, class YMV, class ZMV>
struct Update<XMV, YMV, ZMV, 2, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY>
{
  typedef typename XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<typename XMV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<typename YMV::non_const_value_type> ATB;
  typedef Kokkos::Details::ArithTraits<typename ZMV::non_const_value_type> ATC;

  static void
  update (const typename XMV::non_const_value_type& alpha, const XMV& X,
          const typename YMV::non_const_value_type& beta, const YMV& Y,
          const typename ZMV::non_const_value_type& gamma, const ZMV& Z)
  {
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "Update<rank 2>::update: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::Impl::"
                   "Update<rank 2>::update: Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<ZMV>::value, "KokkosBlas::Impl::"
                   "Update<rank 2>::update: Z is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename ZMV::value_type,
                     typename ZMV::non_const_value_type>::value,
                   "KokkosBlas::Impl::Update<rank 2>::update: Z is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    // Casting to int avoids compiler warnings about comparing
    // different kinds of enum values.
    static_assert ((int) ZMV::rank == (int) XMV::rank &&
                   (int) ZMV::rank == (int) YMV::rank,
                   "KokkosBlas::Impl::Update<rank 2>::update: "
                   "X, Y, and Z must have the same rank.");
    static_assert (ZMV::rank == 2, "KokkosBlas::Impl::Update<rank 2>::update: "
                   "XMV, YMV, and ZMV must have rank 2.");


    #ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::update<> ETI specialization for < %s , %s , %s >\n",typeid(XMV).name(),typeid(YMV).name(),typeid(ZMV).name());
    else {
      printf("KokkosBlas1::update<> non-ETI specialization for < %s , %s , %s >\n",typeid(XMV).name(),typeid(YMV).name(),typeid(ZMV).name());
    }
    #endif

    const size_type numRows = X.extent(0);
    const size_type numCols = X.extent(1);
    int a = 2, b = 2, c = 2;

    if (alpha == ATA::zero ()) {
      a = 0;
    }
    else {
      a = 2;
    }
    if (beta == ATB::zero ()) {
      b = 0;
    }
    else {
      b = 2;
    }
    if (gamma == ATC::zero ()) {
      c = 0;
    }
    else {
      c = 2;
    }

    if (numCols == static_cast<size_type> (1)) {
      // Special case: ZMV has rank 2, but only 1 column.
      // Dispatch to the rank-1 version for better performance.
      auto X_0 = Kokkos::subview (X, Kokkos::ALL (), 0);
      auto Y_0 = Kokkos::subview (Y, Kokkos::ALL (), 0);
      auto Z_0 = Kokkos::subview (Z, Kokkos::ALL (), 0);

      if (numRows * numCols < static_cast<size_type> (INT_MAX)) {
        typedef int index_type;
        V_Update_Generic<decltype (X_0), decltype (Y_0), decltype (Z_0), index_type> (alpha, X_0, beta, Y_0, gamma, Z_0, a, b, c);
      }
      else {
        typedef typename XMV::size_type index_type;
        V_Update_Generic<decltype (X_0), decltype (Y_0), decltype (Z_0), index_type> (alpha, X_0, beta, Y_0, gamma, Z_0, a, b, c);
      }
    }
    else {
      if (numRows * numCols < static_cast<size_type> (INT_MAX)) {
        typedef int index_type;
        MV_Update_Generic<XMV, YMV, ZMV, index_type> (alpha, X, beta, Y, gamma, Z, a, b, c);
      }
      else {
        typedef typename XMV::size_type index_type;
        MV_Update_Generic<XMV, YMV, ZMV, index_type> (alpha, X, beta, Y, gamma, Z, a, b, c);
      }
    }
  }
};

// Partial specialization for XV, YV, and ZV rank-1 Views.
template<class XV, class YV, class ZV>
struct Update<XV, YV, ZV, 1, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY>
{
  typedef typename XV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<typename XV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<typename YV::non_const_value_type> ATB;
  typedef Kokkos::Details::ArithTraits<typename ZV::non_const_value_type> ATC;

  static void
  update (const typename XV::non_const_value_type& alpha, const XV& X,
          const typename YV::non_const_value_type& beta, const YV& Y,
          const typename ZV::non_const_value_type& gamma, const ZV& Z)
  {
    // XV, YV, and ZV must be Kokkos::View specializations.
    static_assert (Kokkos::Impl::is_view<XV>::value, "KokkosBlas::Impl::"
                   "Update<rank 1>::update: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YV>::value, "KokkosBlas::Impl::"
                   "Update<rank 1>::update: Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<ZV>::value, "KokkosBlas::Impl::"
                   "Update<rank 1>::update: Z is not a Kokkos::View.");
    // ZV must be nonconst (else it can't be an output argument).
    static_assert (Kokkos::Impl::is_same<typename ZV::value_type,
                     typename ZV::non_const_value_type>::value,
                   "KokkosBlas::Impl::Update<rank 1>::update: Z is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert ((int) ZV::rank == (int) XV::rank && (int) ZV::rank == (int) YV::rank,
                   "KokkosBlas::Impl::Update<rank 1>::update: "
                   "X, Y, and Z must have the same rank.");
    static_assert (ZV::rank == 1, "KokkosBlas::Impl::Update<rank 1>::update: "
                   "XV, YV, and ZV must have rank 1.");

    #ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::update<> ETI specialization for < %s , %s , %s >\n",typeid(XV).name(),typeid(YV).name(),typeid(ZV).name());
    else {
      printf("KokkosBlas1::update<> non-ETI specialization for < %s , %s , %s >\n",typeid(XV).name(),typeid(YV).name(),typeid(ZV).name());
    }
    #endif
 
    const size_type numRows = X.extent(0);
    const size_type numCols = X.extent(1);
    int a = 2, b = 2, c = 2;

    if (alpha == ATA::zero ()) {
      a = 0;
    }
    else {
      a = 2;
    }
    if (beta == ATB::zero ()) {
      b = 0;
    }
    else {
      b = 2;
    }
    if (gamma == ATC::zero ()) {
      c = 0;
    }
    else {
      c = 2;
    }

    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      typedef int index_type;
      V_Update_Generic<XV, YV, ZV, index_type> (alpha, X, beta, Y, gamma, Z, a, b, c);
    }
    else {
      typedef typename XV::size_type index_type;
      V_Update_Generic<XV, YV, ZV, index_type> (alpha, X, beta, Y, gamma, Z, a, b, c);
    }
  }
};
#endif //!defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY



} // namespace Impl
} // namespace KokkosBlas

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Update for rank == 1.  This is NOT for users!!!
// All the declarations of full specializations go in this header
// file.  We may spread out definitions (see _INST macro below) across
// one or more .cpp files.
//

#define KOKKOSBLAS1_UPDATE_ETI_SPEC_DECL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
extern template struct Update< \
        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        1, false, true>;

#define KOKKOSBLAS1_UPDATE_ETI_SPEC_INST( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
template struct Update< \
        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        1, false, true>;

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Update for rank == 2.  This is NOT for users!!!
// All the declarations of full specializations go in this header
// file.  We may spread out definitions (see _DEF macro below) across
// one or more .cpp files.
//

#define KOKKOSBLAS1_UPDATE_MV_ETI_SPEC_DECL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
extern template struct Update< \
     Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     2, false, true>;

#define KOKKOSBLAS1_UPDATE_MV_ETI_SPEC_INST( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
template struct Update< \
     Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     2, false, true>;


#include<KokkosBlas1_update_tpl_spec_decl.hpp>
#include<generated_specializations_hpp/KokkosBlas1_update_eti_spec_decl.hpp>
#include<generated_specializations_hpp/KokkosBlas1_update_mv_eti_spec_decl.hpp>

#endif // KOKKOSBLAS1_UPDATE_SPEC_HPP_
