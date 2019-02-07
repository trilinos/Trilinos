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
#ifndef KOKKOS_BLAS1_IMPL_DOT_SPEC_HPP_
#define KOKKOS_BLAS1_IMPL_DOT_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>

// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY 
#include <KokkosBlas1_dot_impl.hpp>
#include <KokkosBlas1_dot_mv_impl.hpp>
#endif

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template<class AV, class XV, class YV, int Xrank = XV::rank, int Yrank = YV::rank>
struct dot_eti_spec_avail {
  enum : bool { value = false };
};
}
}

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Dot for rank == 1.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_DOT_ETI_SPEC_AVAIL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
    template<> \
    struct dot_eti_spec_avail< \
        Kokkos::View<SCALAR, LAYOUT, Kokkos::HostSpace, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        1,1> { enum : bool { value = true }; }; \
    template<> \
    struct dot_eti_spec_avail< \
        Kokkos::View<SCALAR, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        1,1> { enum : bool { value = true }; };

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Dot for rank == 2.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_DOT_MV_ETI_SPEC_AVAIL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
    template<> \
    struct dot_eti_spec_avail< \
        Kokkos::View<SCALAR*, typename std::conditional<std::is_same<LAYOUT,Kokkos::LayoutRight>::value, \
                                                        Kokkos::LayoutLeft, LAYOUT>::type, \
                     Kokkos::Device<Kokkos::DefaultHostExecutionSpace,Kokkos::HostSpace>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        2,2> { enum : bool { value = true }; }; \
    template<> \
    struct dot_eti_spec_avail< \
        Kokkos::View<SCALAR*, typename std::conditional<std::is_same<LAYOUT,Kokkos::LayoutRight>::value, \
                                                        Kokkos::LayoutLeft, LAYOUT>::type, \
                     Kokkos::Device<Kokkos::DefaultHostExecutionSpace,Kokkos::HostSpace>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        2,1> { enum : bool { value = true }; }; \
    template<> \
    struct dot_eti_spec_avail< \
        Kokkos::View<SCALAR*, \
                     typename std::conditional<std::is_same<LAYOUT,Kokkos::LayoutRight>::value, \
                                                            Kokkos::LayoutLeft, LAYOUT>::type, \
                     Kokkos::Device<Kokkos::DefaultHostExecutionSpace,Kokkos::HostSpace>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        1,2> { enum : bool { value = true }; };


// Include the actual specialization declarations
#include<KokkosBlas1_dot_tpl_spec_avail.hpp>
#include<generated_specializations_hpp/KokkosBlas1_dot_eti_spec_avail.hpp>
#include<generated_specializations_hpp/KokkosBlas1_dot_mv_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

// Unification layer
template<class RV, class XV, class YV, int XV_Rank = XV::rank, int YV_Rank = YV::rank,
         bool tpl_spec_avail = dot_tpl_spec_avail<RV,XV,YV>::value,
         bool eti_spec_avail = dot_eti_spec_avail<RV,XV,YV>::value>
struct Dot {
  static void dot (const RV&, const XV& R, const YV& X);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//! Full specialization of Dot for single vectors (1-D Views).
template<class RV, class XV, class YV>
struct Dot<RV, XV, YV, 1, 1, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY>
{
  typedef typename YV::size_type size_type;

  static void dot (const RV& R, const XV& X, const YV& Y)
  {
    static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::Impl::"
                   "Dot<1-D>: RV is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XV>::value, "KokkosBlas::Impl::"
                   "Dot<1-D>: XV is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YV>::value, "KokkosBlas::Impl::"
                   "Dot<1-D>: YV is not a Kokkos::View.");
    static_assert (RV::rank == 0, "KokkosBlas::Impl::Dot<1-D>: "
                   "RV is not rank 0.");
    static_assert (XV::rank == 1, "KokkosBlas::Impl::Dot<1-D>: "
                   "XV is not rank 1.");
    static_assert (YV::rank == 1, "KokkosBlas::Impl::Dot<1-D>: "
                   "YV is not rank 1.");
    static_assert (std::is_same<typename RV::value_type,typename RV::non_const_value_type>::value,
                   "KokkosBlas::Dot<1D>: R is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");

    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY?"KokkosBlas::dot[ETI]":"KokkosBlas::dot[noETI]");
    #ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas::dot<> ETI specialization for < %s , %s >\n",typeid(XV).name(),typeid(YV).name());
    else {
      printf("KokkosBlas::dot<> non-ETI specialization for < %s , %s >\n",typeid(XV).name(),typeid(YV).name());
    }
    #endif
    const size_type numElems = X.extent(0);

    if (numElems < static_cast<size_type> (INT_MAX)) {
      typedef int index_type;
      DotFunctor<RV,XV,YV,index_type> f(X,Y);
      f.run("KokkosBlas::dot<1D>",R);
    }
    else {
      typedef int64_t index_type;
      DotFunctor<RV,XV,YV,index_type> f(X,Y);
      f.run("KokkosBlas::dot<1D>",R);
    }
    Kokkos::Profiling::popRegion();
  }
};

template<class RV,class XV, class YV, int X_Rank, int Y_Rank>
struct Dot<RV, XV, YV, X_Rank, Y_Rank, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  typedef typename YV::size_type size_type;

  static void dot (const RV& R, const XV& X, const YV& Y)
  {
    static_assert (Kokkos::Impl::is_view<XV>::value, "KokkosBlas::Impl::"
                   "Dot<2-D>: XV is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YV>::value, "KokkosBlas::Impl::"
                   "Dot<2-D>: YV is not a Kokkos::View.");
    static_assert (RV::rank == 1, "KokkosBlas::Impl::Dot<2-D>: "
                   "RV is not rank 1.");

    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY?"KokkosBlas::dot[ETI]":"KokkosBlas::dot[noETI]");
    #ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::dot<> ETI specialization for < %s , %s , %s >\n",typeid(RV).name(),typeid(XV).name(),typeid(YV).name());
    else {
      printf("KokkosBlas1::dot<> non-ETI specialization for < %s , %s , %s >\n",typeid(RV).name(),typeid(XV).name(),typeid(YV).name());
    }
    #endif

    const size_type numRows = X.extent(0);
    const size_type numCols = X.extent(1);
    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      typedef int index_type;
      Dot_MV<RV,XV,YV,index_type>::dot(R,X,Y);
    }
    else {
      typedef std::int64_t index_type;
      Dot_MV<RV,XV,YV,index_type>::dot(R,X,Y);
    }
    Kokkos::Profiling::popRegion();
  }
};
#endif

}
}

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Dot for rank == 2.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_DOT_ETI_SPEC_DECL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
extern template struct Dot< \
        Kokkos::View<SCALAR, LAYOUT, Kokkos::HostSpace, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        1,1,false,true>; \
extern template struct Dot< \
        Kokkos::View<SCALAR, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        1,1,false,true>;

#define KOKKOSBLAS1_DOT_ETI_SPEC_INST( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
template struct Dot< \
        Kokkos::View<SCALAR, LAYOUT, Kokkos::HostSpace, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        1,1,false,true>; \
template struct Dot< \
        Kokkos::View<SCALAR, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        1,1,false,true>;

//
//
// Macro for definition of full specialization of
// KokkosBlas::Impl::Dot for rank == 2.  This is NOT for users!!!  We
// use this macro in one or more .cpp files in this directory.
//
#define KOKKOSBLAS1_DOT_MV_ETI_SPEC_DECL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
extern template struct Dot< \
        Kokkos::View<SCALAR*, typename std::conditional<std::is_same<LAYOUT,Kokkos::LayoutRight>::value, \
                                                        Kokkos::LayoutLeft, LAYOUT>::type, \
                     Kokkos::Device<Kokkos::DefaultHostExecutionSpace,Kokkos::HostSpace>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        2,2,false,true>; \
extern template struct Dot< \
        Kokkos::View<SCALAR*, typename std::conditional<std::is_same<LAYOUT,Kokkos::LayoutRight>::value, \
                                                        Kokkos::LayoutLeft, LAYOUT>::type, \
                     Kokkos::Device<Kokkos::DefaultHostExecutionSpace,Kokkos::HostSpace>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        2,1,false,true>; \
extern template struct Dot< \
        Kokkos::View<SCALAR*, typename std::conditional<std::is_same<LAYOUT,Kokkos::LayoutRight>::value, \
                                                        Kokkos::LayoutLeft, LAYOUT>::type, \
                     Kokkos::Device<Kokkos::DefaultHostExecutionSpace,Kokkos::HostSpace>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        1,2,false,true>;

#define KOKKOSBLAS1_DOT_MV_ETI_SPEC_INST( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
template struct Dot< \
        Kokkos::View<SCALAR*, typename std::conditional<std::is_same<LAYOUT,Kokkos::LayoutRight>::value, \
                                                        Kokkos::LayoutLeft, LAYOUT>::type, \
                     Kokkos::Device<Kokkos::DefaultHostExecutionSpace,Kokkos::HostSpace>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        2,2,false,true>; \
template struct Dot< \
        Kokkos::View<SCALAR*, typename std::conditional<std::is_same<LAYOUT,Kokkos::LayoutRight>::value, \
                                                        Kokkos::LayoutLeft, LAYOUT>::type, \
                     Kokkos::Device<Kokkos::DefaultHostExecutionSpace,Kokkos::HostSpace>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        2,1,false,true>; \
template struct Dot< \
        Kokkos::View<SCALAR*, typename std::conditional<std::is_same<LAYOUT,Kokkos::LayoutRight>::value, \
                                                         Kokkos::LayoutLeft, LAYOUT>::type, \
                     Kokkos::Device<Kokkos::DefaultHostExecutionSpace,Kokkos::HostSpace>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
        1,2,false,true>;

#include<KokkosBlas1_dot_tpl_spec_decl.hpp>
#include<generated_specializations_hpp/KokkosBlas1_dot_eti_spec_decl.hpp>
#include<generated_specializations_hpp/KokkosBlas1_dot_mv_eti_spec_decl.hpp>

#endif // KOKKOS_BLAS1_MV_IMPL_DOT_HPP_
