// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef KOKKOS_DYN_RANK_VIEW_SACADO_FAD_CONTIGUOUS_HPP
#define KOKKOS_DYN_RANK_VIEW_SACADO_FAD_CONTIGUOUS_HPP

#include "Sacado_ConfigDefs.h"

// This file is setup to always work even when KokkosContainers (which contains
// Kokkos::DynRankView) isn't enabled.

#if defined(HAVE_SACADO_KOKKOS)

#include "Kokkos_DynRankView.hpp"

#if defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

#include "Kokkos_View_Fad.hpp"

namespace Kokkos {
namespace Impl {

template <>
struct DynRankDimTraits<Kokkos::Impl::ViewSpecializeSacadoFadContiguous> {

  enum : size_t{unspecified = ~size_t(0)};

  // Compute the rank of the view from the nonzero dimension arguments.
  // For views of Fad, the rank is one less than the rank determined by the nonzero dimension args
  KOKKOS_INLINE_FUNCTION
  static size_t computeRank( const size_t N0
                           , const size_t N1
                           , const size_t N2
                           , const size_t N3
                           , const size_t N4
                           , const size_t N5
                           , const size_t N6
                           , const size_t N7 )
  {
    return
      (   (N7 == unspecified && N6 == unspecified && N5 == unspecified && N4 == unspecified && N3 == unspecified && N2 == unspecified && N1 == unspecified && N0 == unspecified) ? 0
      : ( (N7 == unspecified && N6 == unspecified && N5 == unspecified && N4 == unspecified && N3 == unspecified && N2 == unspecified && N1 == unspecified) ? 0
      : ( (N7 == unspecified && N6 == unspecified && N5 == unspecified && N4 == unspecified && N3 == unspecified && N2 == unspecified) ? 1
      : ( (N7 == unspecified && N6 == unspecified && N5 == unspecified && N4 == unspecified && N3 == unspecified) ? 2
      : ( (N7 == unspecified && N6 == unspecified && N5 == unspecified && N4 == unspecified) ? 3
      : ( (N7 == unspecified && N6 == unspecified && N5 == unspecified) ? 4
      : ( (N7 == unspecified && N6 == unspecified) ? 5
      : ( (N7 == unspecified) ? 6
      : 7 ) ) ) ) ) ) ) );
  }

  // Compute the rank of the view from the nonzero layout arguments.
  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static size_t computeRank( const Layout& layout )
  {
    return computeRank( layout.dimension[0]
                      , layout.dimension[1]
                      , layout.dimension[2]
                      , layout.dimension[3]
                      , layout.dimension[4]
                      , layout.dimension[5]
                      , layout.dimension[6]
                      , layout.dimension[7] );
  }

  // Compute the rank of the view from the nonzero layout arguments and possible hidden dim
  template <typename Layout, typename ... P>
  KOKKOS_INLINE_FUNCTION
  static size_t computeRank( const ViewCtorProp<P...>& prop, const Layout& layout )
  {
    size_t rank = computeRank( layout.dimension[0]
                      , layout.dimension[1]
                      , layout.dimension[2]
                      , layout.dimension[3]
                      , layout.dimension[4]
                      , layout.dimension[5]
                      , layout.dimension[6]
                      , layout.dimension[7] );

    // Check if has_common_view_alloc_prop; if so, return rank+1, else rank
    enum { test_traits_check = Kokkos::Impl::check_has_common_view_alloc_prop< P... >::value };
    return (test_traits_check == true) ? rank+1 : rank;
  }

  // Create the layout for the rank-7 view.
  // For Fad we have to move the fad dimension to the last (rank 8 since the DynRankView is rank-7)
  // LayoutLeft or LayoutRight
  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if< (std::is_same<Layout , Kokkos::LayoutRight>::value || std::is_same<Layout , Kokkos::LayoutLeft>::value) , Layout >::type createLayout( const Layout& layout )
  {
    Layout l( layout.dimension[0] != unspecified ? layout.dimension[0] : 1
            , layout.dimension[1] != unspecified ? layout.dimension[1] : 1
            , layout.dimension[2] != unspecified ? layout.dimension[2] : 1
            , layout.dimension[3] != unspecified ? layout.dimension[3] : 1
            , layout.dimension[4] != unspecified ? layout.dimension[4] : 1
            , layout.dimension[5] != unspecified ? layout.dimension[5] : 1
            , layout.dimension[6] != unspecified ? layout.dimension[6] : 1
            , layout.dimension[7] != unspecified ? layout.dimension[7] : 1 );
    const unsigned fad_dim = computeRank(layout);
    const size_t fad_size = layout.dimension[fad_dim];
    l.dimension[fad_dim] = 1;
    l.dimension[7] = fad_size;

    return l;
  }

  //LayoutStride
  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if< (std::is_same<Layout , Kokkos::LayoutStride>::value) , Layout>::type createLayout( const Layout& layout )
  {
    Layout      l( layout.dimension[0] != unspecified ? layout.dimension[0] : 1
                 , layout.stride[0]
                 , layout.dimension[1] != unspecified ? layout.dimension[1] : 1
                 , layout.stride[1]
                 , layout.dimension[2] != unspecified ? layout.dimension[2] : 1
                 , layout.stride[2]
                 , layout.dimension[3] != unspecified ? layout.dimension[3] : 1
                 , layout.stride[3]
                 , layout.dimension[4] != unspecified ? layout.dimension[4] : 1
                 , layout.stride[4]
                 , layout.dimension[5] != unspecified ? layout.dimension[5] : 1
                 , layout.stride[5]
                 , layout.dimension[6] != unspecified ? layout.dimension[6] : 1
                 , layout.stride[6]
                 , layout.dimension[7] != unspecified ? layout.dimension[7] : 1
                 , layout.stride[7]
                 );
    const unsigned fad_dim = computeRank(layout);
    const size_t fad_size = layout.dimension[fad_dim];
    l.dimension[fad_dim] = 1;
    l.dimension[7] = fad_size;

    return l;
  }

  // If fad_dim is stored in ViewCtorProp
  // LayoutLeft or LayoutRight
  template <typename Traits, typename ... P>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if< (std::is_same<typename Traits::array_layout , Kokkos::LayoutRight>::value || std::is_same<typename Traits::array_layout , Kokkos::LayoutLeft>::value) , typename Traits::array_layout >::type createLayout( const ViewCtorProp<P...> & arg_prop, const typename Traits::array_layout& layout )
  {
    using Layout = typename Traits::array_layout;

    Layout l( layout.dimension[0] != unspecified ? layout.dimension[0] : 1
            , layout.dimension[1] != unspecified ? layout.dimension[1] : 1
            , layout.dimension[2] != unspecified ? layout.dimension[2] : 1
            , layout.dimension[3] != unspecified ? layout.dimension[3] : 1
            , layout.dimension[4] != unspecified ? layout.dimension[4] : 1
            , layout.dimension[5] != unspecified ? layout.dimension[5] : 1
            , layout.dimension[6] != unspecified ? layout.dimension[6] : 1
            , layout.dimension[7] != unspecified ? layout.dimension[7] : 1 );

    enum { test_traits_check = Kokkos::Impl::check_has_common_view_alloc_prop< P... >::value };
    if (test_traits_check == true) {
      l.dimension[7] = compute_fad_dim_from_alloc_prop<P...>::eval(arg_prop);
    }
    else {
      const unsigned fad_dim = computeRank(layout);
      const size_t fad_size = layout.dimension[fad_dim];
      l.dimension[fad_dim] = 1;
      l.dimension[7] = fad_size;
    }

    return l;
  }

  // If fad_dim is stored in ViewCtorProp
  //LayoutStride
  template <typename Traits, typename ... P>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if< (std::is_same<typename Traits::array_layout , Kokkos::LayoutStride>::value) , typename Traits::array_layout>::type createLayout( const ViewCtorProp<P...> & arg_prop, const typename Traits::array_layout& layout )
  {
    using Layout = typename Traits::array_layout;

    Layout      l( layout.dimension[0] != unspecified ? layout.dimension[0] : 1
                 , layout.stride[0]
                 , layout.dimension[1] != unspecified ? layout.dimension[1] : 1
                 , layout.stride[1]
                 , layout.dimension[2] != unspecified ? layout.dimension[2] : 1
                 , layout.stride[2]
                 , layout.dimension[3] != unspecified ? layout.dimension[3] : 1
                 , layout.stride[3]
                 , layout.dimension[4] != unspecified ? layout.dimension[4] : 1
                 , layout.stride[4]
                 , layout.dimension[5] != unspecified ? layout.dimension[5] : 1
                 , layout.stride[5]
                 , layout.dimension[6] != unspecified ? layout.dimension[6] : 1
                 , layout.stride[6]
                 , layout.dimension[7] != unspecified ? layout.dimension[7] : 1
                 , layout.stride[7]
                 );

    enum { test_traits_check = Kokkos::Impl::check_has_common_view_alloc_prop< P... >::value };
    if (test_traits_check == true) {
      l.dimension[7] = compute_fad_dim_from_alloc_prop<P...>::eval(arg_prop);
    }
    else {
      const unsigned fad_dim = computeRank(layout);
      const size_t fad_size = layout.dimension[fad_dim];
      l.dimension[fad_dim] = 1;
      l.dimension[7] = fad_size;
    }

    return l;
  }


  // Create a view from the given dimension arguments.
  // This is only necessary because the shmem constructor doesn't take a layout.
  template <typename ViewType, typename ViewArg>
  static ViewType createView( const ViewArg& arg
                            , const size_t N0
                            , const size_t N1
                            , const size_t N2
                            , const size_t N3
                            , const size_t N4
                            , const size_t N5
                            , const size_t N6
                            , const size_t N7 )
  {
    typename ViewType::array_layout l( N0, N1, N2, N3, N4, N5, N6, N7 );
    typename ViewType::array_layout l_fad = createLayout(l);
    return ViewType( arg
                   , l_fad.dimension[0]
                   , l_fad.dimension[1]
                   , l_fad.dimension[2]
                   , l_fad.dimension[3]
                   , l_fad.dimension[4]
                   , l_fad.dimension[5]
                   , l_fad.dimension[6]
                   , l_fad.dimension[7] );
  }

};

}} // end Kokkos::Impl

namespace Kokkos {
namespace Impl {

// Specializations for subdynrankview

template< class SrcTraits , class ... Args >
struct ViewMapping
  < typename std::enable_if<(
      std::is_same< typename SrcTraits::specialize ,
                    Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value
      &&
      (
        std::is_same< typename SrcTraits::array_layout
                    , Kokkos::LayoutRight >::value ||
        std::is_same< typename SrcTraits::array_layout
                    , Kokkos::LayoutStride >::value
      )
    ), Kokkos::Impl::DynRankSubviewTag >::type
  , SrcTraits
  , Args ... >
{
private:

  enum
    { RZ = false
    , R0 = bool(is_integral_extent<0,Args...>::value)
    , R1 = bool(is_integral_extent<1,Args...>::value)
    , R2 = bool(is_integral_extent<2,Args...>::value)
    , R3 = bool(is_integral_extent<3,Args...>::value)
    , R4 = bool(is_integral_extent<4,Args...>::value)
    , R5 = bool(is_integral_extent<5,Args...>::value)
    , R6 = bool(is_integral_extent<6,Args...>::value)
    };

  enum { rank = unsigned(R0) + unsigned(R1) + unsigned(R2) + unsigned(R3)
              + unsigned(R4) + unsigned(R5) + unsigned(R6) };

  typedef Kokkos::LayoutContiguous<Kokkos::LayoutStride,SrcTraits::array_layout::scalar_stride> array_layout ;

  typedef typename SrcTraits::value_type  value_type ;

  typedef value_type******* data_type ;

public:

  typedef Kokkos::ViewTraits
    < data_type
    , array_layout
    , typename SrcTraits::device_type
    , typename SrcTraits::memory_traits > traits_type ;

  typedef Kokkos::View
    < data_type
    , array_layout
    , typename SrcTraits::device_type
    , typename SrcTraits::memory_traits > type ;


  template< class MemoryTraits >
  struct apply {

    static_assert( Kokkos::is_memory_traits< MemoryTraits >::value , "" );

    typedef Kokkos::ViewTraits
      < data_type
      , array_layout
      , typename SrcTraits::device_type
      , MemoryTraits > traits_type ;

    typedef Kokkos::View
      < data_type
      , array_layout
      , typename SrcTraits::device_type
      , MemoryTraits > type ;
  };

  template < class Arg0 = int, class Arg1 = int, class Arg2 = int, class Arg3 = int, class Arg4 = int, class Arg5 = int, class Arg6 = int >
  struct ExtentGenerator {
    template <typename dimension>
    KOKKOS_INLINE_FUNCTION
    static SubviewExtents< 7 , rank > generator ( const dimension & dim , Arg0 arg0 = Arg0(), Arg1 arg1 = Arg1(), Arg2 arg2 = Arg2(), Arg3 arg3 = Arg3(), Arg4 arg4 = Arg4(), Arg5 arg5 = Arg5(), Arg6 arg6 = Arg6() )
    {
      return SubviewExtents< 7 , rank >( dim , arg0 , arg1 , arg2 , arg3 ,
                                         arg4 , arg5 , arg6 );
    }
  };

  template < class Arg0 = int, class Arg1 = int, class Arg2 = int, class Arg3 = int, class Arg4 = int, class Arg5 = int, class Arg6 = int >
  struct ArrayExtentGenerator {
    template <typename dimension>
    KOKKOS_INLINE_FUNCTION
    static SubviewExtents< 8 , rank+1 > generator ( const dimension & dim , Arg0 arg0 = Arg0(), Arg1 arg1 = Arg1(), Arg2 arg2 = Arg2(), Arg3 arg3 = Arg3(), Arg4 arg4 = Arg4(), Arg5 arg5 = Arg5(), Arg6 arg6 = Arg6() )
    {
      return SubviewExtents< 8 , rank+1 >( dim , arg0 , arg1 , arg2 , arg3 ,
                                           arg4 , arg5 , arg6 , Kokkos::ALL() );
    }
  };

  typedef DynRankView< value_type , array_layout , typename SrcTraits::device_type , typename SrcTraits::memory_traits >  ret_type;

  template < typename T , class ... P >
  KOKKOS_INLINE_FUNCTION
  static ret_type subview( const unsigned src_rank , Kokkos::DynRankView< T , P...> const & src , Args ... args )
  {

    typedef ViewMapping< traits_type, typename traits_type::specialize >  DstType ;
    typedef typename std::conditional< (rank==0) , ViewDimension<>
      , typename std::conditional< (rank==1) , ViewDimension<0>
      , typename std::conditional< (rank==2) , ViewDimension<0,0>
      , typename std::conditional< (rank==3) , ViewDimension<0,0,0>
      , typename std::conditional< (rank==4) , ViewDimension<0,0,0,0>
      , typename std::conditional< (rank==5) , ViewDimension<0,0,0,0,0>
      , typename std::conditional< (rank==6) , ViewDimension<0,0,0,0,0,0>
      , ViewDimension<0,0,0,0,0,0,0>
      >::type >::type >::type >::type >::type >::type >::type  DstDimType ;
    typedef typename std::conditional< (rank==0) , ViewDimension<0>
      , typename std::conditional< (rank==1) , ViewDimension<0,0>
      , typename std::conditional< (rank==2) , ViewDimension<0,0,0>
      , typename std::conditional< (rank==3) , ViewDimension<0,0,0,0>
      , typename std::conditional< (rank==4) , ViewDimension<0,0,0,0,0>
      , typename std::conditional< (rank==5) , ViewDimension<0,0,0,0,0,0>
      , typename std::conditional< (rank==6) , ViewDimension<0,0,0,0,0,0,0>
      , ViewDimension<0,0,0,0,0,0,0,0>
      >::type >::type >::type >::type >::type >::type >::type  DstArrayDimType ;

    typedef ViewOffset< DstDimType , Kokkos::LayoutStride > dst_offset_type ;
    typedef ViewOffset< DstArrayDimType , Kokkos::LayoutStride > dst_array_offset_type ;
      typedef typename DstType::handle_type  dst_handle_type ;

      ret_type dst ;

      const SubviewExtents< 7 , rank > extents =
        ExtentGenerator< Args ... >::generator(
          src.m_map.m_impl_offset.m_dim , args... ) ;
      const SubviewExtents< 8 , rank+1 > array_extents =
        ArrayExtentGenerator< Args ... >::generator(
          src.m_map.m_array_offset.m_dim , args... ) ;

      dst_offset_type tempdst( src.m_map.m_impl_offset , extents ) ;
      dst_array_offset_type temparraydst(
        src.m_map.m_array_offset , array_extents ) ;

      dst.m_track = src.m_track ;

      dst.m_map.m_impl_offset.m_dim.N0 = tempdst.m_dim.N0 ;
      dst.m_map.m_impl_offset.m_dim.N1 = tempdst.m_dim.N1 ;
      dst.m_map.m_impl_offset.m_dim.N2 = tempdst.m_dim.N2 ;
      dst.m_map.m_impl_offset.m_dim.N3 = tempdst.m_dim.N3 ;
      dst.m_map.m_impl_offset.m_dim.N4 = tempdst.m_dim.N4 ;
      dst.m_map.m_impl_offset.m_dim.N5 = tempdst.m_dim.N5 ;
      dst.m_map.m_impl_offset.m_dim.N6 = tempdst.m_dim.N6 ;

      dst.m_map.m_impl_offset.m_stride.S0 = tempdst.m_stride.S0;
      dst.m_map.m_impl_offset.m_stride.S1 = tempdst.m_stride.S1;
      dst.m_map.m_impl_offset.m_stride.S2 = tempdst.m_stride.S2;
      dst.m_map.m_impl_offset.m_stride.S3 = tempdst.m_stride.S3;
      dst.m_map.m_impl_offset.m_stride.S4 = tempdst.m_stride.S4;
      dst.m_map.m_impl_offset.m_stride.S5 = tempdst.m_stride.S5;
      dst.m_map.m_impl_offset.m_stride.S6 = tempdst.m_stride.S6;

      // Move last non-unit dim and stride to N7/S7 since subview collapses
      // out all singleton dimensions between the last rank and the fad
      // dimension.  Equivalent to:
      //   dst.m_map.m_array_offset.m_dim.N* = temparraydst.m_dim.N*
      //   dst.m_map.m_array_offset.m_dim.N7 = temparraydst.m_dim.N{rank}
      //   dst.m_map.m_array_offset.m_stride.S* = temparraydst.m_stride.S*
      //   dst.m_map.m_array_offset.m_stride.S7 = temparraydst.m_stride.S{rank}
      AssignFadDimStride<rank,0>::eval( dst.m_map.m_array_offset, temparraydst );

      dst.m_track = src.m_track ;

      dst.m_map.m_impl_handle =
        dst_handle_type(
          src.m_map.m_impl_handle +
          src.m_map.m_array_offset( array_extents.domain_offset(0)
                                  , array_extents.domain_offset(1)
                                  , array_extents.domain_offset(2)
                                  , array_extents.domain_offset(3)
                                  , array_extents.domain_offset(4)
                                  , array_extents.domain_offset(5)
                                  , array_extents.domain_offset(6)
                                  , array_extents.domain_offset(7)
            ) );

      dst.m_map.m_fad_size = src.m_map.m_fad_size;
      dst.m_map.m_original_fad_size = src.m_map.m_original_fad_size;
      dst.m_map.m_fad_stride = src.m_map.m_fad_stride;
      dst.m_map.m_fad_index = src.m_map.m_fad_index;

      dst.m_rank = ( src_rank > 0 ? unsigned(R0) : 0 )
                 + ( src_rank > 1 ? unsigned(R1) : 0 )
                 + ( src_rank > 2 ? unsigned(R2) : 0 )
                 + ( src_rank > 3 ? unsigned(R3) : 0 )
                 + ( src_rank > 4 ? unsigned(R4) : 0 )
                 + ( src_rank > 5 ? unsigned(R5) : 0 )
                 + ( src_rank > 6 ? unsigned(R6) : 0 ) ;

      return dst ;
    }
};

template< class SrcTraits , class ... Args >
struct ViewMapping
  < typename std::enable_if<(
      std::is_same< typename SrcTraits::specialize ,
                    Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value
      &&
      std::is_same< typename SrcTraits::array_layout
                   , Kokkos::LayoutLeft >::value
    ), Kokkos::Impl::DynRankSubviewTag >::type
  , SrcTraits
  , Args ... >
{
private:

  enum
    { RZ = false
    , R0 = bool(is_integral_extent<0,Args...>::value)
    , R1 = bool(is_integral_extent<1,Args...>::value)
    , R2 = bool(is_integral_extent<2,Args...>::value)
    , R3 = bool(is_integral_extent<3,Args...>::value)
    , R4 = bool(is_integral_extent<4,Args...>::value)
    , R5 = bool(is_integral_extent<5,Args...>::value)
    , R6 = bool(is_integral_extent<6,Args...>::value)
    };

  enum { rank = unsigned(R0) + unsigned(R1) + unsigned(R2) + unsigned(R3)
              + unsigned(R4) + unsigned(R5) + unsigned(R6) };

  typedef Kokkos::LayoutContiguous<Kokkos::LayoutStride,SrcTraits::array_layout::scalar_stride> array_layout ;

  typedef typename SrcTraits::value_type  value_type ;

  typedef value_type******* data_type ;

public:

  typedef Kokkos::ViewTraits
    < data_type
    , array_layout
    , typename SrcTraits::device_type
    , typename SrcTraits::memory_traits > traits_type ;

  typedef Kokkos::View
    < data_type
    , array_layout
    , typename SrcTraits::device_type
    , typename SrcTraits::memory_traits > type ;


  template< class MemoryTraits >
  struct apply {

    static_assert( Kokkos::is_memory_traits< MemoryTraits >::value , "" );

    typedef Kokkos::ViewTraits
      < data_type
      , array_layout
      , typename SrcTraits::device_type
      , MemoryTraits > traits_type ;

    typedef Kokkos::View
      < data_type
      , array_layout
      , typename SrcTraits::device_type
      , MemoryTraits > type ;
  };

  template < class Arg0 = int, class Arg1 = int, class Arg2 = int, class Arg3 = int, class Arg4 = int, class Arg5 = int, class Arg6 = int >
  struct ExtentGenerator {
    template <typename dimension>
    KOKKOS_INLINE_FUNCTION
    static SubviewExtents< 7 , rank > generator ( const dimension & dim , Arg0 arg0 = Arg0(), Arg1 arg1 = Arg1(), Arg2 arg2 = Arg2(), Arg3 arg3 = Arg3(), Arg4 arg4 = Arg4(), Arg5 arg5 = Arg5(), Arg6 arg6 = Arg6() )
    {
      return SubviewExtents< 7 , rank >( dim , arg0 , arg1 , arg2 , arg3 ,
                                         arg4 , arg5 , arg6 );
    }
  };

  template < class Arg0 = int, class Arg1 = int, class Arg2 = int, class Arg3 = int, class Arg4 = int, class Arg5 = int, class Arg6 = int >
  struct ArrayExtentGenerator {
    template <typename dimension>
    KOKKOS_INLINE_FUNCTION
    static SubviewExtents< 8 , rank+1 > generator ( const dimension & dim , Arg0 arg0 = Arg0(), Arg1 arg1 = Arg1(), Arg2 arg2 = Arg2(), Arg3 arg3 = Arg3(), Arg4 arg4 = Arg4(), Arg5 arg5 = Arg5(), Arg6 arg6 = Arg6() )
    {
      return SubviewExtents< 8 , rank+1 >( dim , Kokkos::ALL() , arg0 , arg1 ,
                                           arg2 , arg3 , arg4 , arg5 , arg6 );
    }
  };

  typedef DynRankView< value_type , array_layout , typename SrcTraits::device_type , typename SrcTraits::memory_traits >  ret_type;

  template < typename T , class ... P >
  KOKKOS_INLINE_FUNCTION
  static ret_type subview( const unsigned src_rank , Kokkos::DynRankView< T , P...> const & src , Args ... args )
  {

    typedef ViewMapping< traits_type, typename traits_type::specialize >  DstType ;
    typedef typename std::conditional< (rank==0) , ViewDimension<>
      , typename std::conditional< (rank==1) , ViewDimension<0>
      , typename std::conditional< (rank==2) , ViewDimension<0,0>
      , typename std::conditional< (rank==3) , ViewDimension<0,0,0>
      , typename std::conditional< (rank==4) , ViewDimension<0,0,0,0>
      , typename std::conditional< (rank==5) , ViewDimension<0,0,0,0,0>
      , typename std::conditional< (rank==6) , ViewDimension<0,0,0,0,0,0>
      , ViewDimension<0,0,0,0,0,0,0>
      >::type >::type >::type >::type >::type >::type >::type  DstDimType ;
    typedef typename std::conditional< (rank==0) , ViewDimension<0>
      , typename std::conditional< (rank==1) , ViewDimension<0,0>
      , typename std::conditional< (rank==2) , ViewDimension<0,0,0>
      , typename std::conditional< (rank==3) , ViewDimension<0,0,0,0>
      , typename std::conditional< (rank==4) , ViewDimension<0,0,0,0,0>
      , typename std::conditional< (rank==5) , ViewDimension<0,0,0,0,0,0>
      , typename std::conditional< (rank==6) , ViewDimension<0,0,0,0,0,0,0>
      , ViewDimension<0,0,0,0,0,0,0,0>
      >::type >::type >::type >::type >::type >::type >::type  DstArrayDimType ;

    typedef ViewOffset< DstDimType , Kokkos::LayoutStride > dst_offset_type ;
    typedef ViewOffset< DstArrayDimType , Kokkos::LayoutStride > dst_array_offset_type ;
      typedef typename DstType::handle_type  dst_handle_type ;

      ret_type dst ;

      const SubviewExtents< 7 , rank > extents =
        ExtentGenerator< Args ... >::generator(
          src.m_map.m_impl_offset.m_dim , args... ) ;
      const SubviewExtents< 8 , rank+1 > array_extents =
        ArrayExtentGenerator< Args ... >::generator(
          src.m_map.m_array_offset.m_dim , args... ) ;

      dst_offset_type tempdst( src.m_map.m_impl_offset , extents ) ;
      dst_array_offset_type temparraydst(
        src.m_map.m_array_offset , array_extents ) ;

      dst.m_track = src.m_track ;

      dst.m_map.m_impl_offset.m_dim.N0 = tempdst.m_dim.N0 ;
      dst.m_map.m_impl_offset.m_dim.N1 = tempdst.m_dim.N1 ;
      dst.m_map.m_impl_offset.m_dim.N2 = tempdst.m_dim.N2 ;
      dst.m_map.m_impl_offset.m_dim.N3 = tempdst.m_dim.N3 ;
      dst.m_map.m_impl_offset.m_dim.N4 = tempdst.m_dim.N4 ;
      dst.m_map.m_impl_offset.m_dim.N5 = tempdst.m_dim.N5 ;
      dst.m_map.m_impl_offset.m_dim.N6 = tempdst.m_dim.N6 ;

      dst.m_map.m_impl_offset.m_stride.S0 = tempdst.m_stride.S0;
      dst.m_map.m_impl_offset.m_stride.S1 = tempdst.m_stride.S1;
      dst.m_map.m_impl_offset.m_stride.S2 = tempdst.m_stride.S2;
      dst.m_map.m_impl_offset.m_stride.S3 = tempdst.m_stride.S3;
      dst.m_map.m_impl_offset.m_stride.S4 = tempdst.m_stride.S4;
      dst.m_map.m_impl_offset.m_stride.S5 = tempdst.m_stride.S5;
      dst.m_map.m_impl_offset.m_stride.S6 = tempdst.m_stride.S6;

      // dst is always LayoutStride, which uses the last (rank-8) index for
      // the fad dimension, thus we need to move its stride/dimension from the
      // first to the last

      dst.m_map.m_array_offset.m_dim.N0 = temparraydst.m_dim.N1 ;
      dst.m_map.m_array_offset.m_dim.N1 = temparraydst.m_dim.N2 ;
      dst.m_map.m_array_offset.m_dim.N2 = temparraydst.m_dim.N3 ;
      dst.m_map.m_array_offset.m_dim.N3 = temparraydst.m_dim.N4 ;
      dst.m_map.m_array_offset.m_dim.N4 = temparraydst.m_dim.N5 ;
      dst.m_map.m_array_offset.m_dim.N5 = temparraydst.m_dim.N6 ;
      dst.m_map.m_array_offset.m_dim.N6 = temparraydst.m_dim.N7 ;
      dst.m_map.m_array_offset.m_dim.N7 = temparraydst.m_dim.N0 ;

      dst.m_map.m_array_offset.m_stride.S0 = temparraydst.m_stride.S1;
      dst.m_map.m_array_offset.m_stride.S1 = temparraydst.m_stride.S2;
      dst.m_map.m_array_offset.m_stride.S2 = temparraydst.m_stride.S3;
      dst.m_map.m_array_offset.m_stride.S3 = temparraydst.m_stride.S4;
      dst.m_map.m_array_offset.m_stride.S4 = temparraydst.m_stride.S5;
      dst.m_map.m_array_offset.m_stride.S5 = temparraydst.m_stride.S6;
      dst.m_map.m_array_offset.m_stride.S6 = temparraydst.m_stride.S7;
      dst.m_map.m_array_offset.m_stride.S7 = temparraydst.m_stride.S0;

      dst.m_track = src.m_track ;

      dst.m_map.m_impl_handle =
        dst_handle_type(
          src.m_map.m_impl_handle +
          src.m_map.m_array_offset( array_extents.domain_offset(0)
                                  , array_extents.domain_offset(1)
                                  , array_extents.domain_offset(2)
                                  , array_extents.domain_offset(3)
                                  , array_extents.domain_offset(4)
                                  , array_extents.domain_offset(5)
                                  , array_extents.domain_offset(6)
                                  , array_extents.domain_offset(7)
            ) );

      dst.m_map.m_fad_size = src.m_map.m_fad_size;
      dst.m_map.m_original_fad_size = src.m_map.m_original_fad_size;
      dst.m_map.m_fad_stride = src.m_map.m_fad_stride;
      dst.m_map.m_fad_index = src.m_map.m_fad_index;

      dst.m_rank = ( src_rank > 0 ? unsigned(R0) : 0 )
                 + ( src_rank > 1 ? unsigned(R1) : 0 )
                 + ( src_rank > 2 ? unsigned(R2) : 0 )
                 + ( src_rank > 3 ? unsigned(R3) : 0 )
                 + ( src_rank > 4 ? unsigned(R4) : 0 )
                 + ( src_rank > 5 ? unsigned(R5) : 0 )
                 + ( src_rank > 6 ? unsigned(R6) : 0 ) ;

      return dst ;
    }
};

/**\brief  Assign compatible Sacado FAD view mappings.
 *
 *  View<FAD>      = View<FAD>
 */
template< class DstTraits , class SrcTraits >
class ViewMapping< DstTraits , SrcTraits ,
  typename std::enable_if<(
    Kokkos::Impl::MemorySpaceAccess
     < typename DstTraits::memory_space
     , typename SrcTraits::memory_space >::assignable
    &&
    // Destination view has FAD only
    std::is_same< typename DstTraits::specialize
                , Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value
    &&
    // Source view has FAD only
    std::is_same< typename SrcTraits::specialize
                , Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value
  ), Kokkos::Impl::ViewToDynRankViewTag >::type >
{
public:

  enum { is_assignable = true };
  enum { is_assignable_data_type = true };

  typedef Kokkos::Impl::SharedAllocationTracker  TrackType ;
  typedef ViewMapping< DstTraits , typename DstTraits::specialize >  DstType ;
  typedef ViewMapping< SrcTraits , typename SrcTraits::specialize >  SrcFadType ;

  template < typename DT , typename ... DP , typename ST , typename ... SP >
  KOKKOS_INLINE_FUNCTION static
  void assign( Kokkos::DynRankView< DT , DP... > & dst
             , const Kokkos::View< ST , SP... >& src )
    {
      static_assert(
        (
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutLeft >::value ||
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutRight >::value ||
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutStride >::value
        )
        &&
        (
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutLeft >::value ||
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutRight >::value ||
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutStride >::value
        )
        , "View of FAD requires LayoutLeft, LayoutRight, or LayoutStride" );

      static_assert(
        std::is_same< typename DstTraits::value_type
                    , typename SrcTraits::value_type >::value ||
        std::is_same< typename DstTraits::value_type
                    , typename SrcTraits::const_value_type >::value ,
        "View assignment must have same value type or const = non-const" );

      const bool is_left =
        std::is_same<typename DstTraits::array_layout,Kokkos::LayoutLeft>::value;
      typedef typename DstType::offset_type dst_offset_type;
      typedef typename DstType::array_offset_type dst_array_offset_type;
      if (is_left) {
        dst.m_map.m_array_offset =
          dst_array_offset_type(std::integral_constant<unsigned,0>(),
                                src.m_map.m_array_offset.layout() );
      }
      else {
        dst.m_map.m_array_offset =
          dst_array_offset_type(std::integral_constant<unsigned,0>(),
                                permute_fad_layout(src.m_map.m_array_offset.layout(),
                                                   SrcTraits::rank) );
      }
      dst.m_map.m_impl_offset =
        dst_offset_type(std::integral_constant<unsigned,0>(),
                        src.m_map.m_impl_offset.layout() );

      dst.m_map.m_impl_handle = src.m_map.m_impl_handle ;
      dst.m_rank = src.rank ;

      dst.m_map.m_fad_size = src.m_map.m_fad_size ;
      dst.m_map.m_original_fad_size = src.m_map.m_original_fad_size;
      dst.m_map.m_fad_stride = src.m_map.m_fad_stride ;
      dst.m_map.m_fad_index = src.m_map.m_fad_index;
    }
};

/**\brief  Assign compatible Sacado FAD view mappings.
 *
 *  View<ordinary> = View<FAD>
 */
template< class DstTraits , class SrcTraits >
class ViewMapping< DstTraits , SrcTraits ,
  typename std::enable_if<(
    Kokkos::Impl::MemorySpaceAccess
     < typename DstTraits::memory_space
     , typename SrcTraits::memory_space >::assignable
    &&
    // Destination view has ordinary
    std::is_same< typename DstTraits::specialize , void >::value
    &&
    // Source view has FAD only
    std::is_same< typename SrcTraits::specialize
                , Kokkos::Impl::ViewSpecializeSacadoFadContiguous >::value
  ), Kokkos::Impl::ViewToDynRankViewTag >::type >
{
public:

  enum { is_assignable = true };
  enum { is_assignable_data_type = true };

  typedef Kokkos::Impl::SharedAllocationTracker  TrackType ;
  typedef ViewMapping< DstTraits , typename DstTraits::specialize >  DstType ;
  typedef ViewMapping< SrcTraits , typename SrcTraits::specialize >  SrcFadType ;

  template < typename DT , typename ... DP , typename ST , typename ... SP >
  KOKKOS_INLINE_FUNCTION static
  void assign( Kokkos::DynRankView< DT , DP... > & dst
             , const Kokkos::View< ST , SP... >& src )
    {
      static_assert(
        (
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutLeft >::value ||
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutRight >::value ||
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutStride >::value
        )
        &&
        (
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutLeft >::value ||
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutRight >::value ||
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutStride >::value
        )
        , "View of FAD requires LayoutLeft, LayoutRight, or LayoutStride" );

      static_assert(
        std::is_same< typename DstTraits::value_type
                    , typename SrcTraits::value_type >::value ||
        std::is_same< typename DstTraits::value_type
                    , typename SrcTraits::const_value_type >::value ,
        "View assignment must have same value type or const = non-const" );

      dst.m_map.m_impl_offset.m_dim = src.m_map.m_array_offset.m_dim;
      dst.m_map.m_impl_offset.m_stride = src.m_map.m_array_offset.m_stride ;

      dst.m_map.m_impl_handle = src.m_map.m_impl_handle ;
      dst.m_rank = src.rank ;
    }
};

}} //end Kokkos::Impl

#endif //defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

#endif // defined(HAVE_SACADO_KOKKOS)

#endif /* #ifndef KOKKOS_DYN_RANK_VIEW_SACADO_FAD_HPP */
