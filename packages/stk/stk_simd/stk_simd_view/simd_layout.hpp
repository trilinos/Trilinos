// Copyright (c) 2013, Sandia Corporation.
 // Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 // the U.S. Government retains certain rights in this software.
 // 
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are
 // met:
 // 
 //     * Redistributions of source code must retain the above copyright
 //       notice, this list of conditions and the following disclaimer.
 // 
 //     * Redistributions in binary form must reproduce the above
 //       copyright notice, this list of conditions and the following
 //       disclaimer in the documentation and/or other materials provided
 //       with the distribution.
 // 
 //     * Neither the name of Sandia Corporation nor the names of its
 //       contributors may be used to endorse or promote products derived
 //       from this software without specific prior written permission.
 // 
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef STK_SIMD_LAYOUT_H
#define STK_SIMD_LAYOUT_H

#include <stk_simd/Simd.hpp>


namespace stk {
namespace simd {
template <typename RealType>
struct LayoutRight;
template <typename RealType>
struct LayoutLeft;
}
}

namespace Kokkos {

namespace Impl {
template < class Dimension , class Layout , typename Enable >
struct ViewOffset;
template < class Dimension , typename RealType >
struct ViewOffset<Dimension, stk::simd::LayoutRight<RealType>, void>;
template < class Dimension , typename RealType >
struct ViewOffset<Dimension, stk::simd::LayoutLeft<RealType>, void>;
} // namespace Impl

} // namespace Kokkos

#include "Kokkos_Core.hpp"

constexpr int simd_double_width = stk::simd::ndoubles;
constexpr int simd_double_width_log_2 = simd_double_width==1 ? 0
  : simd_double_width==2 ? 1
  : simd_double_width==4 ? 2
  : simd_double_width==8 ? 3
  : -1;
constexpr int simd_double_width_m_one = simd_double_width-1;

constexpr int simd_float_width = stk::simd::nfloats;
constexpr int simd_float_width_log_2 = simd_float_width==1 ? 0
  : simd_float_width==2 ? 1
  : simd_float_width==4 ? 2
  : simd_float_width==8 ? 3
  : simd_float_width==16 ? 4
  : -1;
constexpr int simd_float_width_m_one = simd_float_width-1;

template <typename T>
struct SimdSizeTraits {
  static constexpr int simd_width = 1;
  static constexpr int simd_width_log_2 = 0;
  static constexpr int simd_width_m_one = 0;
};

template <>
struct SimdSizeTraits<double> {
  static constexpr int simd_width = simd_double_width;
  static constexpr int simd_width_log_2 = simd_double_width_log_2;
  static constexpr int simd_width_m_one = simd_double_width_m_one;
};

template <>
struct SimdSizeTraits<float> {
  static constexpr int simd_width = simd_float_width;
  static constexpr int simd_width_log_2 = simd_float_width_log_2;
  static constexpr int simd_width_m_one = simd_float_width_m_one;
};


template <typename T>
KOKKOS_INLINE_FUNCTION
constexpr size_t simd_pad(const size_t i) {
  return i + SimdSizeTraits<T>::simd_width_m_one - (i-1)%SimdSizeTraits<T>::simd_width;
}

namespace stk {
namespace simd {

template <typename RealType>
struct LayoutRight {
  //! Tag this class as a kokkos array layout
  typedef LayoutRight array_layout ;

  size_t dimension[ Kokkos::ARRAY_LAYOUT_MAX_RANK ];

  LayoutRight( LayoutRight const & ) = default ;
  LayoutRight( LayoutRight && ) = default ;
  LayoutRight & operator = ( LayoutRight const & ) = default ;
  LayoutRight & operator = ( LayoutRight && ) = default ;

  KOKKOS_INLINE_FUNCTION
  constexpr
  LayoutRight( size_t N0 = 0 , size_t N1 = 0 , size_t N2 = 0 , size_t N3 = 0
                 , size_t N4 = 0 , size_t N5 = 0 , size_t N6 = 0 , size_t N7 = 0 )
    : dimension { N0 , N1 , N2 , N3 , N4 , N5 , N6 , N7 } {}
};


template <typename RealType>
struct LayoutLeft {
  //! Tag this class as a kokkos array layout
  typedef LayoutLeft array_layout ;

  size_t dimension[ Kokkos::ARRAY_LAYOUT_MAX_RANK ];

  LayoutLeft( LayoutLeft const & ) = default ;
  LayoutLeft( LayoutLeft && ) = default ;
  LayoutLeft & operator = ( LayoutLeft const & ) = default ;
  LayoutLeft & operator = ( LayoutLeft && ) = default ;

  KOKKOS_INLINE_FUNCTION
  constexpr
  LayoutLeft( size_t N0 = 0 , size_t N1 = 0 , size_t N2 = 0 , size_t N3 = 0 ,
                  size_t N4 = 0 , size_t N5 = 0 , size_t N6 = 0 , size_t N7 = 0 )
    : dimension { N0 , N1 , N2 , N3 , N4 , N5 , N6 , N7 } {}
};

}}

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
// LayoutRight AND ( 1 >= rank OR 0 == rank_dynamic ) : padding / no striding
template < class Dimension, typename RealType >
struct ViewOffset< Dimension , stk::simd::LayoutRight<RealType> , void>
{
  using is_mapping_plugin = std::true_type ;
  using is_regular        = std::false_type ;

  typedef size_t          size_type ;
  typedef Dimension       dimension_type ;
  typedef stk::simd::LayoutRight<RealType> array_layout ;

  static constexpr int simd_width       = SimdSizeTraits<RealType>::simd_width;
  static constexpr int simd_width_log_2 = SimdSizeTraits<RealType>::simd_width_log_2;
  static constexpr int simd_width_m_one = SimdSizeTraits<RealType>::simd_width_m_one;

  enum { has_padding = true };

  dimension_type m_dim ;

  //----------------------------------------

  // rank 1
  template< typename I0 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0 ) const { return i0 ; }

  // rank 2
  template < typename I0 , typename I1 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0 , I1 const & i1 ) const
    {
      return simd_width * ( i1 + m_dim.N1 * ( i0 >> simd_width_log_2 ) ) + (i0 & simd_width_m_one);
    }

  //rank 3
  template < typename I0, typename I1, typename I2 >
  KOKKOS_INLINE_FUNCTION
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2 ) const
    {
      return simd_width * ( i2 + m_dim.N2 * ( i1 + m_dim.N1 * (i0 >> simd_width_log_2) ) ) + (i0 & simd_width_m_one);
    }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  constexpr array_layout layout() const
    {
      return array_layout( m_dim.N0 , m_dim.N1 , m_dim.N2 , m_dim.N3
                         , m_dim.N4 , m_dim.N5 , m_dim.N6 , m_dim.N7 );
    }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION constexpr
  typename std::enable_if< std::is_integral<iType>::value , size_t >::type
  extent( const iType & r ) const
    { return m_dim.extent(r); }

  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_0() const { return m_dim.N0 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_1() const { return m_dim.N1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_2() const { return m_dim.N2 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_3() const { return m_dim.N3 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_4() const { return m_dim.N4 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_5() const { return m_dim.N5 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_6() const { return m_dim.N6 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_7() const { return m_dim.N7 ; }

  /* Cardinality of the domain index space */
  KOKKOS_INLINE_FUNCTION
  constexpr size_type size() const
    { return m_dim.N0 * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 * m_dim.N5 * m_dim.N6 * m_dim.N7 ; }

  /* Span of the range space */
  KOKKOS_INLINE_FUNCTION
  constexpr size_type span() const
    { return simd_pad<RealType>( m_dim.N0 ) * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 * m_dim.N5 * m_dim.N6 * m_dim.N7 ; }

  KOKKOS_INLINE_FUNCTION constexpr bool span_is_contiguous() const { return true ; }

  /* Strides of dimensions */
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_0() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_1() const { return m_dim.N0 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_2() const { return m_dim.N0 * m_dim.N1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_3() const { return m_dim.N0 * m_dim.N1 * m_dim.N2 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_4() const { return m_dim.N0 * m_dim.N1 * m_dim.N2 * m_dim.N3 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_5() const { return m_dim.N0 * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_6() const { return m_dim.N0 * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 * m_dim.N5 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_7() const { return m_dim.N0 * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 * m_dim.N5 * m_dim.N6 ; }

  // Stride with [ rank ] value is the total length
  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void stride( iType * const s ) const
    {
      Kokkos::abort("Can't stride a Kokkos::View with LayoutRight");
      s[0] = 1 ;
      if ( 0 < dimension_type::rank ) { s[1] = m_dim.N0 ; }
      if ( 1 < dimension_type::rank ) { s[2] = s[1] * m_dim.N1 ; }
      if ( 2 < dimension_type::rank ) { s[3] = s[2] * m_dim.N2 ; }
      if ( 3 < dimension_type::rank ) { s[4] = s[3] * m_dim.N3 ; }
      if ( 4 < dimension_type::rank ) { s[5] = s[4] * m_dim.N4 ; }
      if ( 5 < dimension_type::rank ) { s[6] = s[5] * m_dim.N5 ; }
      if ( 6 < dimension_type::rank ) { s[7] = s[6] * m_dim.N6 ; }
      if ( 7 < dimension_type::rank ) { s[8] = s[7] * m_dim.N7 ; }
    }

  //----------------------------------------

  ViewOffset() = default ;
  ViewOffset( const ViewOffset & ) = default ;
  ViewOffset & operator = ( const ViewOffset & ) = default ;

  template< unsigned TrivialScalarSize >
  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset( std::integral_constant<unsigned,TrivialScalarSize> const & , stk::simd::LayoutRight<RealType> const & arg_layout )
    : m_dim( arg_layout.dimension[0], arg_layout.dimension[1], arg_layout.dimension[2], arg_layout.dimension[3], arg_layout.dimension[4], arg_layout.dimension[5], arg_layout.dimension[6], arg_layout.dimension[7] )
    {}

  template< class DimRHS >
  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset( const ViewOffset< DimRHS , stk::simd::LayoutRight<RealType> , void > & rhs )
    : m_dim( rhs.m_dim.N0 , rhs.m_dim.N1 , rhs.m_dim.N2 , rhs.m_dim.N3
           , rhs.m_dim.N4 , rhs.m_dim.N5 , rhs.m_dim.N6 , rhs.m_dim.N7 )
    {
      static_assert( int(DimRHS::rank) == int(dimension_type::rank) , "ViewOffset assignment requires equal rank" );
      // Also requires equal static dimensions ...
    }

  template< class DimRHS >
  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset( const ViewOffset< DimRHS , Kokkos::LayoutRight , void > & rhs )
    : m_dim( rhs.m_dim.N0, 0, 0, 0, 0, 0, 0, 0 )
    {
      static_assert( DimRHS::rank == 1 && dimension_type::rank == 1 && dimension_type::rank_dynamic == 1
                   , "ViewOffset LayoutRight and LayoutRight are only compatible when rank == 1" );
    }

  template< class DimRHS >
  KOKKOS_INLINE_FUNCTION
  ViewOffset( const ViewOffset< DimRHS , Kokkos::LayoutStride , void > & rhs )
    : m_dim( rhs.m_dim.N0, 0, 0, 0, 0, 0, 0, 0 )
    {
      static_assert( DimRHS::rank == 1 && dimension_type::rank == 1 && dimension_type::rank_dynamic == 1
                   , "ViewOffset LayoutRight and LayoutStride are only compatible when rank == 1" );
      if ( rhs.m_stride.S0 != 1 ) {
        Kokkos::abort("Kokkos::ViewOffset assignment of LayoutRight from LayoutStride  requires stride == 1" );
      }
    }

};

//----------------------------------------------------------------------------
// LayoutLeft AND ( 1 >= rank OR 0 == rank_dynamic ) : padding / no striding
template < class Dimension, typename RealType >
struct ViewOffset< Dimension , stk::simd::LayoutLeft<RealType> , void>
{
  using is_mapping_plugin = std::true_type ;
  using is_regular        = std::false_type ;

  typedef size_t             size_type ;
  typedef Dimension          dimension_type ;
  typedef stk::simd::LayoutLeft<RealType> array_layout ;

  static constexpr int simd_width       = SimdSizeTraits<RealType>::simd_width;
  static constexpr int simd_width_log_2 = SimdSizeTraits<RealType>::simd_width_log_2;
  static constexpr int simd_width_m_one = SimdSizeTraits<RealType>::simd_width_m_one;

  enum { has_padding = true };

  dimension_type m_dim ;

  size_type S0;

  //----------------------------------------

  // rank 1
  template< typename I0 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0 ) const {
    return i0 ;
  }

  // rank 2
  template < typename I0 , typename I1 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0 , I1 const & i1 ) const {
    return i0 + S0 * i1;
  }

  //rank 3
  template < typename I0, typename I1, typename I2 >
  KOKKOS_INLINE_FUNCTION
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2 ) const {
    return i0 + S0 * i1 + S0 * m_dim.N1 * i2;
  }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  constexpr array_layout layout() const
    {
      return array_layout( m_dim.N0 , m_dim.N1 , m_dim.N2 , m_dim.N3
                         , m_dim.N4 , m_dim.N5 , m_dim.N6 , m_dim.N7 );
    }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION constexpr
  typename std::enable_if< std::is_integral<iType>::value , size_t >::type
  extent( const iType & r ) const
    { return m_dim.extent(r); }

  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_0() const { return m_dim.N0 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_1() const { return m_dim.N1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_2() const { return m_dim.N2 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_3() const { return m_dim.N3 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_4() const { return m_dim.N4 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_5() const { return m_dim.N5 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_6() const { return m_dim.N6 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_7() const { return m_dim.N7 ; }

  /* Cardinality of the domain index space */
  KOKKOS_INLINE_FUNCTION
  constexpr size_type size() const
    { return m_dim.N0 * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 * m_dim.N5 * m_dim.N6 * m_dim.N7 ; }

  /* Span of the range space */
  KOKKOS_INLINE_FUNCTION
  constexpr size_type span() const
    { return S0 * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 * m_dim.N5 * m_dim.N6 * m_dim.N7 ; }

  KOKKOS_INLINE_FUNCTION constexpr bool span_is_contiguous() const { return true ; }

  /* Strides of dimensions */
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_0() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_1() const { return S0 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_2() const { return S0 * m_dim.N1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_3() const { return S0 * m_dim.N1 * m_dim.N2 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_4() const { return S0 * m_dim.N1 * m_dim.N2 * m_dim.N3 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_5() const { return S0 * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_6() const { return S0 * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 * m_dim.N5 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_7() const { return S0 * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 * m_dim.N5 * m_dim.N6 ; }

  // Stride with [ rank ] value is the total length
  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void stride( iType * const s ) const
    {
      Kokkos::abort("Can't stride a Kokkos::View with LayoutLeft");
      s[0] = 1 ;
      if ( 0 < dimension_type::rank ) { s[1] = S0 ; }
      if ( 1 < dimension_type::rank ) { s[2] = s[1] * m_dim.N1 ; }
      if ( 2 < dimension_type::rank ) { s[3] = s[2] * m_dim.N2 ; }
      if ( 3 < dimension_type::rank ) { s[4] = s[3] * m_dim.N3 ; }
      if ( 4 < dimension_type::rank ) { s[5] = s[4] * m_dim.N4 ; }
      if ( 5 < dimension_type::rank ) { s[6] = s[5] * m_dim.N5 ; }
      if ( 6 < dimension_type::rank ) { s[7] = s[6] * m_dim.N6 ; }
      if ( 7 < dimension_type::rank ) { s[8] = s[7] * m_dim.N7 ; }
    }

  //----------------------------------------

  ViewOffset() = default ;
  ViewOffset( const ViewOffset & ) = default ;
  ViewOffset & operator = ( const ViewOffset & ) = default ;

  template< unsigned TrivialScalarSize >
  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset
  ( std::integral_constant<unsigned,TrivialScalarSize> const & , stk::simd::LayoutLeft<RealType> const & arg_layout
    )
    : m_dim( arg_layout.dimension[0], arg_layout.dimension[1], arg_layout.dimension[2], arg_layout.dimension[3], arg_layout.dimension[4], arg_layout.dimension[5], arg_layout.dimension[6], arg_layout.dimension[7] )
    , S0( simd_pad<RealType>(arg_layout.dimension[0]) )
    {}

  template< class DimRHS >
  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset( const ViewOffset< DimRHS , stk::simd::LayoutLeft<RealType> , void > & rhs )
    : m_dim( rhs.m_dim.N0 , rhs.m_dim.N1 , rhs.m_dim.N2 , rhs.m_dim.N3 ,
             rhs.m_dim.N4 , rhs.m_dim.N5 , rhs.m_dim.N6 , rhs.m_dim.N7 )
    , S0( simd_pad<RealType>( rhs.m_dim.N0 ) )
    {
      static_assert( int(DimRHS::rank) == int(dimension_type::rank) , "ViewOffset assignment requires equal rank" );
      // Also requires equal static dimensions ...
    }

  template< class DimRHS >
  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset( const ViewOffset< DimRHS , LayoutLeft , void > & rhs )
    : m_dim( rhs.m_dim.N0, 0, 0, 0, 0, 0, 0, 0 )
    , S0( simd_pad<RealType>( rhs.m_dim.N0 ) )
    {
      static_assert( DimRHS::rank == 1 && dimension_type::rank == 1 && dimension_type::rank_dynamic == 1
                   , "ViewOffset LayoutLeft and LayoutLeft are only compatible when rank == 1" );
    }

  template< class DimRHS >
  KOKKOS_INLINE_FUNCTION
  ViewOffset( const ViewOffset< DimRHS , Kokkos::LayoutStride , void > & rhs )
    : m_dim( rhs.m_dim.N0, 0, 0, 0, 0, 0, 0, 0 )
    , S0( simd_pad<RealType>( rhs.m_dim.N0 ) )
    {
      static_assert( DimRHS::rank == 1 && dimension_type::rank == 1 && dimension_type::rank_dynamic == 1
                   , "ViewOffset LayoutLeft and LayoutStride are only compatible when rank == 1" );
      if ( rhs.m_stride.S0 != 1 ) {
        Kokkos::abort("Kokkos::ViewOffset assignment of LayoutLeft from LayoutStride  requires stride == 1" );
      }
    }

};

} // namespace Impl
} // namespace Kokkos

#endif
