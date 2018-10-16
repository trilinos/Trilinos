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
#ifndef KOKKOS_BLAS1_AXPBY_SPEC_HPP_
#define KOKKOS_BLAS1_AXPBY_SPEC_HPP_

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_InnerProductSpaceTraits.hpp"

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include<KokkosBlas1_axpby_impl.hpp>
#include<KokkosBlas1_axpby_mv_impl.hpp>
#endif

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template<class AV, class XMV, class BV, class YMV, int rank = YMV::Rank>
struct axpby_eti_spec_avail {
  enum : bool { value = false };
};
}
}

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Axpby for rank == 1.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_AXPBY_ETI_SPEC_AVAIL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
    template<> \
    struct axpby_eti_spec_avail< \
         SCALAR, \
         Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         SCALAR, \
         Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         1> { enum : bool { value = true }; };

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Axpby for rank == 2.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_AXPBY_MV_ETI_SPEC_AVAIL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
    template<> \
    struct axpby_eti_spec_avail< \
         SCALAR, \
         Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         SCALAR, \
         Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         2> { enum : bool { value = true }; }; \
    template<> \
    struct axpby_eti_spec_avail< \
         Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,\
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,\
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
         2> { enum : bool { value = true }; };

// Include the actual specialization declarations
#include<KokkosBlas1_axpby_tpl_spec_avail.hpp>
#include<generated_specializations_hpp/KokkosBlas1_axpby_eti_spec_avail.hpp>
#include<generated_specializations_hpp/KokkosBlas1_axpby_mv_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

//
// axpby
//

/// \brief Implementation of KokkosBlas::axpby for (multi)vectors.
///
/// Compute any of the following, depending on the types of the input
/// arguments of axpxy():
///
/// 1. Y(i,j) = av(j)*X(i,j) + bv(j)*Y(i,j) (if R, X, and Y are 2-D,
///    and av and bv are 1-D)
///
/// 2. Y(i,j) = av*X(i,j) + bv*Y(i,j) (if R, X, and Y are 2-D,
///    and av and bv are scalars)
///
/// 3. Y(i) = av()*X(i) + bv()*Y(i) (if R, X, and Y are 1-D, and av
///    and bv are 0-D Views (not scalars))
///
/// 4. Y(i) = av*X(i) + bv*Y(i) (if R, X, and Y are 1-D, and av and bv
///    are scalars)
///
/// Any <i>scalar</i> coefficient of zero has BLAS semantics of
/// ignoring the corresponding (multi)vector entry.  This does NOT
/// apply to coefficients in av and bv vectors, if they are used.
template<class AV, class XMV, class BV, class YMV, int rank = YMV::Rank,
         bool tpl_spec_avail = axpby_tpl_spec_avail<AV,XMV,BV,YMV>::value,
         bool eti_spec_avail = axpby_eti_spec_avail<AV,XMV,BV,YMV>::value>
struct Axpby {
  static void axpby (const AV& av, const XMV& X, const BV& bv, const YMV& Y);
};

template<class AV, class XMV, class BV, class YMV>
struct Axpby<AV,XMV,BV,YMV,0,true,true> {
  static void axpby (const AV& av, const XMV& X, const BV& bv, const YMV& Y) {
    static_assert(YMV::Rank==0,"Oh My God");
  }
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
// Full specialization for XMV and YMV rank-2 Views.
template<class AV, class XMV, class BV, class YMV>
struct Axpby<AV, XMV, BV, YMV, 2, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY>
{
  typedef typename YMV::size_type size_type;

  static void
  axpby (const AV& av, const XMV& X, const BV& bv, const YMV& Y)
  {
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "Axpby<rank-2>::axpby: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::Impl::"
                   "Axpby<rank-2>::axpby: Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename YMV::value_type,
                     typename YMV::non_const_value_type>::value,
                   "KokkosBlas::Impl::Axpby<rank-2>::axpby: Y is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert ((int) YMV::Rank == (int) XMV::Rank,
                   "KokkosBlas::Impl::Axpby<rank-2>::axpby (MV): "
                   "X and Y must have the same rank.");
    static_assert (YMV::Rank == 2, "KokkosBlas::Impl::Axpby<rank-2>::axpby: "
                   "X and Y must have rank 2.");

    #ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::axpby<> ETI specialization for < %s , %s , %s , %s >\n",typeid(AV).name(),typeid(XMV).name(),typeid(BV).name(),typeid(YMV).name());
    else {
      printf("KokkosBlas1::axpby<> non-ETI specialization for < %s , %s , %s , %s >\n",typeid(AV).name(),typeid(XMV).name(),typeid(BV).name(),typeid(YMV).name());
    }
    #endif

    const size_type numRows = X.extent(0);
    const size_type numCols = X.extent(1);
    int a = 2, b = 2;
    if (av.extent(0) == 0) {
      a = 0;
    }
    if (bv.extent(0) == 0) {
      b = 0;
    }

    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      typedef int index_type;
      typedef typename std::conditional<std::is_same<typename XMV::array_layout,Kokkos::LayoutLeft>::value,
        Axpby_MV_Invoke_Right<AV, XMV, BV, YMV, index_type>,
        Axpby_MV_Invoke_Left<AV, XMV, BV, YMV, index_type> >::type Axpby_MV_Invoke_Layout;
      Axpby_MV_Invoke_Layout::run(av, X, bv, Y, a, b);
    }
    else {
      typedef typename XMV::size_type index_type;
      typedef typename std::conditional<std::is_same<typename XMV::array_layout,Kokkos::LayoutLeft>::value,
        Axpby_MV_Invoke_Right<AV, XMV, BV, YMV, index_type>,
        Axpby_MV_Invoke_Left<AV, XMV, BV, YMV, index_type> >::type Axpby_MV_Invoke_Layout;
      Axpby_MV_Invoke_Layout::run(av, X, bv, Y, a, b);
    }
  }
};

// Partial specialization for XMV, and YMV rank-2 Views,
// and AV and BV scalars.
template<class XMV, class YMV>
struct Axpby<typename XMV::non_const_value_type, XMV,
             typename YMV::non_const_value_type, YMV, 2, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY >
{
  typedef typename XMV::non_const_value_type AV;
  typedef typename YMV::non_const_value_type BV;
  typedef typename YMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<typename XMV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<typename YMV::non_const_value_type> ATB;

  static void
  axpby (const AV& alpha, const XMV& X, const BV& beta, const YMV& Y)
  {
    static_assert (Kokkos::Impl::is_view<XMV>::value,
                   "KokkosBlas::Impl::Axpby::axpby (MV): "
                   "X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YMV>::value,
                   "KokkosBlas::Impl::Axpby::axpby (MV): "
                   "Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename YMV::value_type,
                     typename YMV::non_const_value_type>::value,
                   "KokkosBlas::Impl::Axpby::axpby (MV): Y is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert ((int) YMV::Rank == (int) XMV::Rank,
                   "KokkosBlas::Impl::Axpby::axpby (MV): "
                   "X and Y must have the same rank.");
    static_assert (YMV::Rank == 2, "KokkosBlas::Impl::Axpby::axpby (MV): "
                   "X and Y must have rank 2.");


    #ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::axpby<> ETI specialization for < %s , %s , %s , %s >\n",typeid(AV).name(),typeid(XMV).name(),typeid(BV).name(),typeid(YMV).name());
    else {
      printf("KokkosBlas1::axpby<> non-ETI specialization for < %s , %s , %s , %s >\n",typeid(AV).name(),typeid(XMV).name(),typeid(BV).name(),typeid(YMV).name());
    }
    #endif

    const size_type numRows = X.extent(0);
    const size_type numCols = X.extent(1);
    int a, b;
    if (alpha == ATA::zero ()) {
      a = 0;
    }
#if KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2
    else if (alpha == -ATA::one ()) {
      a = -1;
    }
    else if (alpha == ATA::one ()) {
      a = 1;
    }
#endif // KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2
    else {
      a = 2;
    }
    if (beta == ATB::zero ()) {
      b = 0;
    }
#if KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2
    else if (beta == -ATB::one ()) {
      b = -1;
    }
    else if (beta == ATB::one ()) {
      b = 1;
    }
#endif // KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2
    else {
      b = 2;
    }


    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      typedef int index_type;
      typedef typename std::conditional<std::is_same<typename XMV::array_layout,Kokkos::LayoutLeft>::value,
        Axpby_MV_Invoke_Right<AV, XMV, BV, YMV, index_type>,
        Axpby_MV_Invoke_Left<AV, XMV, BV, YMV, index_type> >::type Axpby_MV_Invoke_Layout;
      Axpby_MV_Invoke_Layout::run(alpha, X,
                                                          beta, Y, a, b);
    }
    else {
      typedef typename XMV::size_type index_type;
      typedef typename std::conditional<std::is_same<typename XMV::array_layout,Kokkos::LayoutLeft>::value,
        Axpby_MV_Invoke_Right<AV, XMV, BV, YMV, index_type>,
        Axpby_MV_Invoke_Left<AV, XMV, BV, YMV, index_type> >::type Axpby_MV_Invoke_Layout;
      Axpby_MV_Invoke_Layout::run(alpha, X,
                                                          beta, Y, a, b);
    }
  }
};

// Partial specialization for XV and YV rank-1 Views,
// and AV and BV scalars.
template<class XV, class YV>
struct Axpby<typename XV::non_const_value_type, XV,
             typename YV::non_const_value_type, YV, 1, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY>
{
  typedef typename XV::non_const_value_type AV;
  typedef typename YV::non_const_value_type BV;
  typedef typename YV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<typename XV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<typename YV::non_const_value_type> ATB;

  static void
  axpby (const AV& alpha, const XV& X, const BV& beta, const YV& Y)
  {
    static_assert (Kokkos::Impl::is_view<XV>::value, "KokkosBlas::Impl::"
                   "Axpby<rank-1>::axpby: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YV>::value, "KokkosBlas::Impl::"
                   "Axpby<rank-1>::axpby: Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename YV::value_type,
                     typename YV::non_const_value_type>::value,
                   "KokkosBlas::Impl::Axpby<rank-1>::axpby: Y is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert ((int) YV::Rank == (int) XV::Rank, "KokkosBlas::Impl::"
                   "Axpby<rank-1>::axpby: X and Y must have the same rank.");
    static_assert (YV::Rank == 1, "KokkosBlas::Impl::Axpby<rank-1>::axpby: "
                   "X and Y must have rank 1.");

    #ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::axpby<> ETI specialization for < %s , %s , %s , %s >\n",typeid(AV).name(),typeid(XV).name(),typeid(BV).name(),typeid(YV).name());
    else {
      printf("KokkosBlas1::axpby<> non-ETI specialization for < %s , %s , %s , %s >\n",typeid(AV).name(),typeid(XV).name(),typeid(BV).name(),typeid(YV).name());
    }
    #endif

    const size_type numRows = X.extent(0);
    int a = 2;
    if (alpha == ATA::zero ()) {
      a = 0;
    }
#if KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2
    else if (alpha == -ATA::one ()) {
      a = -1;
    }
    else if (alpha == ATA::one ()) {
      a = 1;
    }
#endif // KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2

    int b = 2;
    if (beta == ATB::zero ()) {
      b = 0;
    }
#if KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2
    else if (beta == -ATB::one ()) {
      b = -1;
    }
    else if (beta == ATB::one ()) {
      b = 1;
    }
#endif // KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2

    if (numRows < static_cast<size_type> (INT_MAX)) {
      typedef int index_type;
      Axpby_Generic<typename XV::non_const_value_type, XV,
        typename YV::non_const_value_type, YV,
        index_type> (alpha, X, beta, Y, 0, a, b);
    }
    else {
      typedef typename XV::size_type index_type;
      Axpby_Generic<typename XV::non_const_value_type, XV,
        typename YV::non_const_value_type, YV,
        index_type> (alpha, X, beta, Y, 0, a, b);
    }
  }
};
#endif //!defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY



} // namespace Impl
} // namespace KokkosBlas

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Axpby for rank == 1.  This is NOT for users!!!
// All the declarations of full specializations go in this header
// file.  We may spread out definitions (see _INST macro below) across
// one or more .cpp files.
//

#define KOKKOSBLAS1_AXPBY_ETI_SPEC_DECL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
extern template struct Axpby< \
        SCALAR, \
        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        SCALAR, \
        Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        1, false, true>;

#define KOKKOSBLAS1_AXPBY_ETI_SPEC_INST( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
template struct Axpby< \
        SCALAR, \
        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        SCALAR, \
        Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        1, false, true>;

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Axpby for rank == 2.  This is NOT for users!!!
// All the declarations of full specializations go in this header
// file.  We may spread out definitions (see _DEF macro below) across
// one or more .cpp files.
//

#define KOKKOSBLAS1_AXPBY_MV_ETI_SPEC_DECL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
extern template struct Axpby< \
     SCALAR, \
     Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     SCALAR, \
     Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     2, false, true>; \
extern template struct Axpby< \
     Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,\
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,\
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     2, false, true>;

#define KOKKOSBLAS1_AXPBY_MV_ETI_SPEC_INST( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
template struct Axpby< \
     SCALAR, \
     Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     SCALAR, \
     Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     2, false, true>; \
template struct Axpby< \
     Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,\
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,\
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     2, false, true>;


#include<KokkosBlas1_axpby_tpl_spec_decl.hpp>
#include<generated_specializations_hpp/KokkosBlas1_axpby_eti_spec_decl.hpp>
#include<generated_specializations_hpp/KokkosBlas1_axpby_mv_eti_spec_decl.hpp>

#endif // KOKKOS_BLAS1_MV_IMPL_AXPBY_HPP_
