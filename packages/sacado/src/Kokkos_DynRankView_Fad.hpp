// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef KOKKOS_DYN_RANK_VIEW_SACADO_FAD_HPP
#define KOKKOS_DYN_RANK_VIEW_SACADO_FAD_HPP

#include "Sacado_ConfigDefs.h"

// This file is setup to always work even when KokkosContainers (which contains
// Kokkos::DynRankView) isn't enabled.

#if defined(HAVE_SACADO_KOKKOSCONTAINERS)

#include "Kokkos_DynRankView.hpp"

#if defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

#include "Kokkos_View_Fad.hpp"

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <>
struct DynRankDimTraits<ViewSpecializeSacadoFad> {

  enum : size_t{unspecified = ~size_t(0)};

  // Compute the rank of the view from the nonzero dimension arguments.
  // For views of Fad, the rank is one less than the rank determined by the nonzero dimension args
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

// Utility class that handles calculation of the stride in Fad subview
template <unsigned> struct AssignFadStride {};
template <> struct AssignFadStride<0u> {
  template <typename Src, typename Dst>
  KOKKOS_INLINE_FUNCTION
  static void eval(Dst& dst, const Src& src) {
    dst.m_stride.S0 = 0 ;
    dst.m_stride.S1 = 0 ;
    dst.m_stride.S2 = 0 ;
    dst.m_stride.S3 = 0 ;
    dst.m_stride.S4 = 0 ;
    dst.m_stride.S5 = 0 ;
    dst.m_stride.S6 = 0 ;
    dst.m_stride.S7 = src.m_stride.S0 ;
  }
};
template <> struct AssignFadStride<1u> {
  template <typename Src, typename Dst>
  KOKKOS_INLINE_FUNCTION
  static void eval(Dst& dst, const Src& src) {
    dst.m_stride.S0 = src.m_stride.S0 ;
    dst.m_stride.S1 = 0 ;
    dst.m_stride.S2 = 0 ;
    dst.m_stride.S3 = 0 ;
    dst.m_stride.S4 = 0 ;
    dst.m_stride.S5 = 0 ;
    dst.m_stride.S6 = 0 ;
    dst.m_stride.S7 = src.m_stride.S1 ;
  }
};
template <> struct AssignFadStride<2u> {
  template <typename Src, typename Dst>
  KOKKOS_INLINE_FUNCTION
  static void eval(Dst& dst, const Src& src) {
    dst.m_stride.S0 = src.m_stride.S0 ;
    dst.m_stride.S1 = src.m_stride.S1 ;
    dst.m_stride.S2 = 0 ;
    dst.m_stride.S3 = 0 ;
    dst.m_stride.S4 = 0 ;
    dst.m_stride.S5 = 0 ;
    dst.m_stride.S6 = 0 ;
    dst.m_stride.S7 = src.m_stride.S2 ;
  }
};
template <> struct AssignFadStride<3u> {
  template <typename Src, typename Dst>
  KOKKOS_INLINE_FUNCTION
  static void eval(Dst& dst, const Src& src) {
    dst.m_stride.S0 = src.m_stride.S0 ;
    dst.m_stride.S1 = src.m_stride.S1 ;
    dst.m_stride.S2 = src.m_stride.S2 ;
    dst.m_stride.S3 = 0 ;
    dst.m_stride.S4 = 0 ;
    dst.m_stride.S5 = 0 ;
    dst.m_stride.S6 = 0 ;
    dst.m_stride.S7 = src.m_stride.S3 ;
  }
};
template <> struct AssignFadStride<4u> {
  template <typename Src, typename Dst>
  KOKKOS_INLINE_FUNCTION
  static void eval(Dst& dst, const Src& src) {
    dst.m_stride.S0 = src.m_stride.S0 ;
    dst.m_stride.S1 = src.m_stride.S1 ;
    dst.m_stride.S2 = src.m_stride.S2 ;
    dst.m_stride.S3 = src.m_stride.S3 ;
    dst.m_stride.S4 = 0 ;
    dst.m_stride.S5 = 0 ;
    dst.m_stride.S6 = 0 ;
    dst.m_stride.S7 = src.m_stride.S4 ;
  }
};
template <> struct AssignFadStride<5u> {
  template <typename Src, typename Dst>
  KOKKOS_INLINE_FUNCTION
  static void eval(Dst& dst, const Src& src) {
    dst.m_stride.S0 = src.m_stride.S0 ;
    dst.m_stride.S1 = src.m_stride.S1 ;
    dst.m_stride.S2 = src.m_stride.S2 ;
    dst.m_stride.S3 = src.m_stride.S3 ;
    dst.m_stride.S4 = src.m_stride.S4 ;
    dst.m_stride.S5 = 0 ;
    dst.m_stride.S6 = 0 ;
    dst.m_stride.S7 = src.m_stride.S5 ;
  }
};
template <> struct AssignFadStride<6u> {
  template <typename Src, typename Dst>
  KOKKOS_INLINE_FUNCTION
  static void eval(Dst& dst, const Src& src) {
    dst.m_stride.S0 = src.m_stride.S0 ;
    dst.m_stride.S1 = src.m_stride.S1 ;
    dst.m_stride.S2 = src.m_stride.S2 ;
    dst.m_stride.S3 = src.m_stride.S3 ;
    dst.m_stride.S4 = src.m_stride.S4 ;
    dst.m_stride.S5 = src.m_stride.S5 ;
    dst.m_stride.S6 = 0 ;
    dst.m_stride.S7 = src.m_stride.S6 ;
  }
};
template <> struct AssignFadStride<7u> {
  template <typename Src, typename Dst>
  KOKKOS_INLINE_FUNCTION
  static void eval(Dst& dst, const Src& src) {
    dst.m_stride.S0 = src.m_stride.S0 ;
    dst.m_stride.S1 = src.m_stride.S1 ;
    dst.m_stride.S2 = src.m_stride.S2 ;
    dst.m_stride.S3 = src.m_stride.S3 ;
    dst.m_stride.S4 = src.m_stride.S4 ;
    dst.m_stride.S5 = src.m_stride.S5 ;
    dst.m_stride.S6 = src.m_stride.S6 ;
    dst.m_stride.S7 = src.m_stride.S7 ;
  }
};

template <unsigned> struct AssignDim7 {
  template <typename Src, typename Dst>
  KOKKOS_INLINE_FUNCTION
  static void eval(Dst& dst, const Src& src) {}
};
template <> struct AssignDim7<0u> {
  template <typename Src, typename Dst>
  KOKKOS_INLINE_FUNCTION
  static void eval(Dst& dst, const Src& src) {
    dst.m_dim.N7 = src.m_dim.N7;
  }
};

// Specializations for subdynrankview
template< class SrcTraits , class ... Args >
struct ViewMapping
  < typename std::enable_if<(
      std::is_same< typename SrcTraits::specialize ,
                    ViewSpecializeSacadoFad >::value
      &&
      (
        std::is_same< typename SrcTraits::array_layout
                    , Kokkos::LayoutLeft >::value ||
        std::is_same< typename SrcTraits::array_layout
                    , Kokkos::LayoutRight >::value ||
        std::is_same< typename SrcTraits::array_layout
                    , Kokkos::LayoutStride >::value
      ) 
    ), DynRankSubviewTag >::type
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

  typedef Kokkos::LayoutStride array_layout ;

  typedef typename SrcTraits::value_type  value_type ;

  typedef value_type******* data_type ;

public:

  typedef Kokkos::Experimental::ViewTraits
    < data_type
    , array_layout
    , typename SrcTraits::device_type
    , typename SrcTraits::memory_traits > traits_type ;

  typedef Kokkos::Experimental::View
    < data_type
    , array_layout
    , typename SrcTraits::device_type
    , typename SrcTraits::memory_traits > type ;


  template< class MemoryTraits >
  struct apply {

    static_assert( Kokkos::Impl::is_memory_traits< MemoryTraits >::value , "" );

    typedef Kokkos::Experimental::ViewTraits
      < data_type
      , array_layout
      , typename SrcTraits::device_type
      , MemoryTraits > traits_type ;

    typedef Kokkos::Experimental::View
      < data_type
      , array_layout
      , typename SrcTraits::device_type
      , MemoryTraits > type ;
  };


  //typedef typename SrcTraits::dimension dimension ;

  template < class Arg0 = int, class Arg1 = int, class Arg2 = int, class Arg3 = int, class Arg4 = int, class Arg5 = int, class Arg6 = int >
  struct ExtentGenerator {
    template <typename dimension>
    KOKKOS_INLINE_FUNCTION
    static SubviewExtents< 8 , rank+1 > generator ( const dimension & dim , Arg0 arg0 = Arg0(), Arg1 arg1 = Arg1(), Arg2 arg2 = Arg2(), Arg3 arg3 = Arg3(), Arg4 arg4 = Arg4(), Arg5 arg5 = Arg5(), Arg6 arg6 = Arg6() )
    {
       return SubviewExtents< 8 , rank+1 >( dim , arg0 , arg1 , arg2 , arg3 , arg4 , arg5 , arg6 , Kokkos::ALL() );
    }
  };

  typedef DynRankView< value_type , array_layout , typename SrcTraits::device_type , typename SrcTraits::memory_traits >  ret_type;

  template < typename T , class ... P >
  KOKKOS_INLINE_FUNCTION
  static ret_type subview( const unsigned src_rank , Kokkos::Experimental::View< T******* , P...> const & src
                    , Args ... args )
  {

    typedef ViewMapping< traits_type, void >  DstType ;
    typedef ViewMapping< SrcTraits, void> SrcType;
    enum { FadStaticDim = SrcType::FadStaticDimension };

    typedef typename std::conditional< (rank==0) , ViewDimension<FadStaticDim>
      , typename std::conditional< (rank==1) , ViewDimension<0,FadStaticDim>
      , typename std::conditional< (rank==2) , ViewDimension<0,0,FadStaticDim>
      , typename std::conditional< (rank==3) , ViewDimension<0,0,0,FadStaticDim>
      , typename std::conditional< (rank==4) , ViewDimension<0,0,0,0,FadStaticDim>
      , typename std::conditional< (rank==5) , ViewDimension<0,0,0,0,0,FadStaticDim>
      , typename std::conditional< (rank==6) , ViewDimension<0,0,0,0,0,0,FadStaticDim>
      , ViewDimension<0,0,0,0,0,0,0,FadStaticDim>
      >::type >::type >::type >::type >::type >::type >::type  DstDimType ;

      typedef ViewOffset< DstDimType , Kokkos::LayoutStride > dst_offset_type ;
      typedef typename DstType::handle_type  dst_handle_type ;

      ret_type dst ;

      const SubviewExtents< 8 , rank+1 > extents =
        ExtentGenerator< Args ... >::generator( src.m_map.m_offset.m_dim , args... ) ;

      dst_offset_type tempdst( src.m_map.m_offset , extents ) ;

      dst.m_track = src.m_track ;

      dst.m_map.m_offset.m_dim.N0 = tempdst.m_dim.N0 ;
      dst.m_map.m_offset.m_dim.N1 = tempdst.m_dim.N1 ;
      dst.m_map.m_offset.m_dim.N2 = tempdst.m_dim.N2 ;
      dst.m_map.m_offset.m_dim.N3 = tempdst.m_dim.N3 ;
      dst.m_map.m_offset.m_dim.N4 = tempdst.m_dim.N4 ;
      dst.m_map.m_offset.m_dim.N5 = tempdst.m_dim.N5 ;
      dst.m_map.m_offset.m_dim.N6 = tempdst.m_dim.N6 ;

      // Do this except for when Fad dim is static
      // dst.m_map.m_offset.m_dim.N7 = tempdst.m_dim.N7 ;
      AssignDim7<FadStaticDim>::eval( dst.m_map.m_offset, tempdst );

      // Move last non-unit stride to S7
      // dst.m_map.m_offset.m_stride.S* = tempdst.m_stride.S*;
      // dst.m_map.m_offset.m_stride.S7 = tempdst.m_stride.S{rank}
      AssignFadStride<rank>::eval( dst.m_map.m_offset, tempdst );

      dst.m_track = src.m_track ;

      dst.m_map.m_handle =
        dst_handle_type( src.m_map.m_handle +
                         src.m_map.m_offset( extents.domain_offset(0)
                                           , extents.domain_offset(1)
                                           , extents.domain_offset(2)
                                           , extents.domain_offset(3)
                                           , extents.domain_offset(4)
                                           , extents.domain_offset(5)
                                           , extents.domain_offset(6)
                                           , extents.domain_offset(7)
                         ) );

      dst.m_map.m_fad_size = src.m_map.m_fad_size;

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

}
}
}

#endif //defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

namespace Kokkos {

// Overload of dimension_scalar() for all dynamic-rank views
template <typename T, typename ... P>
KOKKOS_INLINE_FUNCTION
constexpr unsigned
dimension_scalar(const Experimental::DynRankView<T,P...>& view) {
  return dimension_scalar(view.ConstDownCast());
}

}

#endif // defined(HAVE_SACADO_KOKKOSCONTAINERS)

#endif /* #ifndef KOKKOS_DYN_RANK_VIEW_SACADO_FAD_HPP */
