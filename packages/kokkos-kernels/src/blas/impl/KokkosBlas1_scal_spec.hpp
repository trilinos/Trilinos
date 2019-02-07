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
#ifndef KOKKOS_BLAS1_IMPL_SCAL_SPEC_HPP_
#define KOKKOS_BLAS1_IMPL_SCAL_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>

// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY 
#include <KokkosBlas1_scal_impl.hpp>
#include <KokkosBlas1_scal_mv_impl.hpp>
#endif

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template<class RV, class AV, class XV, int Xrank = XV::rank>
struct scal_eti_spec_avail {
  enum : bool { value = false };
};
}
}

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Scal for rank == 1.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_SCAL_ETI_SPEC_AVAIL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
    template<> \
    struct scal_eti_spec_avail< \
        Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        SCALAR, \
        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        1> { enum : bool { value = true }; };

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Scal for rank == 2.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_SCAL_MV_ETI_SPEC_AVAIL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
    template<> \
    struct scal_eti_spec_avail< \
        Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        2> { enum : bool { value = true }; }; \
    template<> \
    struct scal_eti_spec_avail< \
        Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        SCALAR, \
        Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        2> { enum : bool { value = true }; };



// Include the actual specialization declarations
#include<KokkosBlas1_scal_tpl_spec_avail.hpp>
#include<generated_specializations_hpp/KokkosBlas1_scal_eti_spec_avail.hpp>
#include<generated_specializations_hpp/KokkosBlas1_scal_mv_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

// Unification layer
template<class RV, class AV, class XV, int XV_Rank = XV::rank,
         bool tpl_spec_avail = scal_tpl_spec_avail<RV,AV,XV>::value,
         bool eti_spec_avail = scal_eti_spec_avail<RV,AV,XV>::value>
struct Scal {
  static void scal (const RV& R, const AV& A, const XV& X);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//! Full specialization of Scal for single vectors (1-D Views).
template<class RV, class XV>
struct Scal<RV,  typename XV::non_const_value_type, XV, 1, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY>
{
  typedef typename XV::non_const_value_type AV;
  typedef typename XV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<typename XV::non_const_value_type> ATA;

  static void
  scal (const RV& R, const AV& alpha, const XV& X)
  {
    static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::Impl::"
                   "Scal<1-D>: RV is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XV>::value, "KokkosBlas::Impl::"
                   "Scal<1-D>: XV is not a Kokkos::View.");
    static_assert (RV::rank == 1, "KokkosBlas::Impl::Scal<1-D>: "
                   "RV is not rank 1.");
    static_assert (XV::rank == 1, "KokkosBlas::Impl::Scal<1-D>: "
                   "XV is not rank 1.");
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY?"KokkosBlas::scal[ETI]":"KokkosBlas::scal[noETI]");

    #ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::scal<1D> ETI specialization for < %s , %s , %s >\n",typeid(RV).name(),typeid(AV).name(),typeid(XV).name());
    else
      printf("KokkosBlas1::scal<1D> non-ETI specialization for < %s , %s , %s >\n",typeid(RV).name(),typeid(AV).name(),typeid(XV).name());
    #endif

    const size_type numRows = X.extent(0);
    int a = 2;
    if (alpha == ATA::zero ()) {
      a = 0;
    }
    else if (alpha == -ATA::one ()) {
      a = -1;
    }
    else if (alpha == ATA::one ()) {
      a = 1;
    }

    if (numRows < static_cast<size_type> (INT_MAX)) {
      typedef int index_type;
      V_Scal_Generic<RV, AV, XV, index_type> (R, alpha, X, a);
    }
    else {
      typedef typename XV::size_type index_type;
      V_Scal_Generic<RV, AV, XV, index_type> (R, alpha, X, a);
    }
    Kokkos::Profiling::popRegion();
  }
};

/// \brief Partial specialization of Scal for 2-D Views and 1-D View AV.
///
/// Compute any of the following:
///
/// 1. R(i,j) = a*X(i,j) for a in -1,0,1
/// 2. R(i,j) = alpha(j)*X(i,j)
template<class RMV, class AV, class XMV>
struct Scal<RMV, AV, XMV, 2, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  typedef typename XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<typename XMV::non_const_value_type> ATA;

  static void
  scal (const RMV& R, const AV& av, const XMV& X)
  {
    static_assert (Kokkos::Impl::is_view<RMV>::value, "KokkosBlas::Impl::"
                   "Scal<2-D>: RMV is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<AV>::value, "KokkosBlas::Impl::"
                   "Scal<2-D>: AV is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "Scal<2-D>: XMV is not a Kokkos::View.");
    static_assert (RMV::rank == 2, "KokkosBlas::Impl::Scal<2-D>: "
                   "RMV is not rank 2.");
    static_assert (AV::rank == 1, "KokkosBlas::Impl::Scal<2-D>: "
                   "AV is not rank 1.");
    static_assert (XMV::rank == 2, "KokkosBlas::Impl::Scal<2-D>: "
                   "XMV is not rank 2.");
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY?"KokkosBlas::scal[ETI]":"KokkosBlas::scal[noETI]");
    #ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::scal<2D> ETI specialization for < %s , %s , %s >\n",typeid(RMV).name(),typeid(AV).name(),typeid(XMV).name());
    else
      printf("KokkosBlas1::scal<2D> non-ETI specialization for < %s , %s , %s >\n",typeid(RMV).name(),typeid(AV).name(),typeid(XMV).name());
    #endif

    const size_type numRows = X.extent(0);
    const size_type numCols = X.extent(1);
    const int a = (av.extent(0) == 0) ? 0 : 2;
    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      typedef int index_type;
      MV_Scal_Invoke_Left<RMV, AV, XMV, index_type> (R, av, X, a);
    }
    else {
      typedef typename XMV::size_type index_type;
      MV_Scal_Invoke_Left<RMV, AV, XMV, index_type> (R, av, X, a);
    }
    Kokkos::Profiling::popRegion();
  }
};

/// \brief Partial specialization of Scal for 2-D Views and scalar AV.
///
/// Compute any of the following:
///
/// 1. R(i,j) = a*X(i,j) for a in -1,0,1
/// 2. R(i,j) = alpha*X(i,j)
template<class RMV, class XMV>
struct Scal<RMV, typename XMV::non_const_value_type, XMV, 2, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  typedef typename XMV::non_const_value_type AV;
  typedef typename XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<typename XMV::non_const_value_type> ATA;

  static void
  scal (const RMV& R, const AV& alpha, const XMV& X)
  {
    static_assert (Kokkos::Impl::is_view<RMV>::value, "KokkosBlas::Impl::"
                   "Scal<2-D, AV=scalar>: RMV is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "Scal<2-D, AV=scalar>: XMV is not a Kokkos::View.");
    static_assert (RMV::rank == 2, "KokkosBlas::Impl::Scal<2-D, AV=scalar>: "
                   "RMV is not rank 2.");
    static_assert (XMV::rank == 2, "KokkosBlas::Impl::Scal<2-D, AV=scalar>: "
                   "XMV is not rank 2.");
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY?"KokkosBlas::scal[ETI]":"KokkosBlas::scal[noETI]");

    #ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::scal<2D> ETI specialization for < %s , %s , %s >\n",typeid(RMV).name(),typeid(AV).name(),typeid(XMV).name());
    else
      printf("KokkosBlas1::scal<2D> non-ETI specialization for < %s , %s , %s >\n",typeid(RMV).name(),typeid(AV).name(),typeid(XMV).name());
    #endif

    const size_type numRows = X.extent(0);
    const size_type numCols = X.extent(1);
    int a = 2;
    if (alpha == ATA::zero ()) {
      a = 0;
    }
    else if (alpha == -ATA::one ()) {
      a = -1;
    }
    else if (alpha == ATA::one ()) {
      a = 1;
    }

    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      typedef int index_type;
      MV_Scal_Invoke_Left<RMV, typename XMV::non_const_value_type, XMV,
        index_type> (R, alpha, X, a);
    }
    else {
      typedef typename XMV::size_type index_type;
      MV_Scal_Invoke_Left<RMV, typename XMV::non_const_value_type, XMV,
        index_type> (R, alpha, X, a);
    }
    Kokkos::Profiling::popRegion();
  }
};
#endif

}
}

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Scal for rank == 2.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_SCAL_ETI_SPEC_DECL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
extern template struct Scal< \
        Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        SCALAR, \
        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        1,false,true>;

#define KOKKOSBLAS1_SCAL_ETI_SPEC_INST( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
template struct Scal< \
        Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        SCALAR, \
        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        1,false,true>;

//
//
// Macro for definition of full specialization of
// KokkosBlas::Impl::Scal for rank == 2.  This is NOT for users!!!  We
// use this macro in one or more .cpp files in this directory.
//
#define KOKKOSBLAS1_SCAL_MV_ETI_SPEC_DECL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
extern template struct Scal< \
        Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        2,false,true>; \
extern template struct Scal< \
        Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        SCALAR, \
        Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        2,false,true>;

#define KOKKOSBLAS1_SCAL_MV_ETI_SPEC_INST( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
template struct Scal< \
        Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, \
                     Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        2,false,true>; \
template struct Scal< \
        Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        SCALAR, \
        Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        2,false,true>;

#include<KokkosBlas1_scal_tpl_spec_decl.hpp>
#include<generated_specializations_hpp/KokkosBlas1_scal_eti_spec_decl.hpp>
#include<generated_specializations_hpp/KokkosBlas1_scal_mv_eti_spec_decl.hpp>

#endif // KOKKOS_BLAS1_MV_IMPL_SCAL_HPP_
