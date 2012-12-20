/*
//@HEADER
// ************************************************************************
// 
//    KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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
// Questions? Contact H. Carter Edwards (hcedwar@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSARRAY_CRSPRODUCTTENSORLEGENDRE_HPP
#define KOKKOSARRAY_CRSPRODUCTTENSORLEGENDRE_HPP

#include <cmath>
#include <utility>
#include <vector>
#include <stdexcept>
#include <KokkosArray_Macros.hpp>
#include <KokkosArray_View.hpp>
#include <KokkosArray_BlockCrsMatrix.hpp>
#include <KokkosArray_ProductTensorLegendre.hpp>

#if 0

namespace KokkosArray {

//----------------------------------------------------------------------------

template< typename TensorScalar , class Device >
struct CrsProductTensorLegendre {
  typedef View< unsigned* ,     Device >  array_unsigned_type ;
  typedef View< TensorScalar* , Device >  array_scalar_type ;

  typedef typename array_unsigned_type::HostMirror array_unsigned_host_type ;
  typedef typename array_scalar_type  ::HostMirror array_scalar_host_type ;

  array_unsigned_type  m_entry_offset ;
  array_unsigned_type  m_coordinate ;
  array_scalar_type    m_value ;
  unsigned             m_dimension ;
  unsigned             m_max_row_width ;
  unsigned             m_nonzero_count ;
  unsigned             m_multiply_add_flops ;

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned dimension() const { return m_dimension ; }

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned max_row_width() const { return m_max_row_width ; }

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned multiply_add_flops() const { return m_multiply_add_flops ; }

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned nonzero_count() const { return m_nonzero_count ; }

  CrsProductTensorLegendre()
  : m_entry_offset()
  , m_coordinate()
  , m_value()
  , m_dimension(0)
  , m_max_row_width(0)
  , m_nonzero_count(0)
  , m_multiply_add_flops(0)
  {}

  CrsProductTensorLegendre( const CrsProductTensorLegendre & rhs )
  : m_entry_offset( rhs.m_entry_offset )
  , m_coordinate(   rhs.m_coordinate )
  , m_value(        rhs.m_value )
  , m_dimension(    rhs.m_dimension )
  , m_max_row_width(rhs.m_max_row_width )
  , m_nonzero_count(rhs.m_nonzero_count )
  , m_multiply_add_flops(rhs.m_multiply_add_flops)
  {}

  CrsProductTensorLegendre & operator = ( const CrsProductTensorLegendre & rhs )
  {
    m_entry_offset = rhs.m_entry_offset ;
    m_coordinate   = rhs.m_coordinate ;
    m_value        = rhs.m_value ;
    m_dimension    = rhs.m_dimension ;
    m_max_row_width= rhs.m_max_row_width ;
    m_nonzero_count= rhs.m_nonzero_count ;
    m_multiply_add_flops=rhs.m_multiply_add_flops;
    return *this ;
  }

  CrsProductTensorLegendre( const std::vector<unsigned> & variable_poly_degree ,
                            const unsigned maximum_poly_degree )
  : m_entry_offset()
  , m_coordinate()
  , m_value()
  , m_dimension(0)
  , m_max_row_width(0)
  , m_nonzero_count(0)
  , m_multiply_add_flops(0)
  {
    enum { Align = Impl::is_same<Device,Cuda>::value ? 32 : 1 };

    const KokkosArray::TripleProductTensorLegendreCombinatorialEvaluation
      combinatorial( maximum_poly_degree ,
                     variable_poly_degree );

    m_dimension = combinatorial.bases_count();

    if ( ( 1 << 16 ) < m_dimension ) {
      throw std::runtime_error("CrsProductTensorLegendre tensor dimension too large");
    }

    unsigned entry_count = 0 ;

    // Add to output vector:
    m_multiply_add_flops = m_dimension ;

    for ( unsigned i = 0 ; i < m_dimension ; ++i ) {

      unsigned row_entry_count = 0 ;

      for ( unsigned j = 0 ; j < m_dimension ; ++j ) {
        if ( combinatorial.is_non_zero(i,j,j) ) {
          ++row_entry_count ;
          ++m_nonzero_count ;
          m_multiply_add_flops += 3 ;
        }
      }

      if ( row_entry_count % Align ) { row_entry_count += Align - row_entry_count % Align ; }

      for ( unsigned j = 0 ; j < m_dimension ; ++j ) {
        for ( unsigned k = j+1 ; k < m_dimension ; ++k ) {
          if ( combinatorial.is_non_zero(i,j,k) ) {
            ++row_entry_count ;
            ++m_nonzero_count ;
            m_multiply_add_flops += 5 ;
          }
        }
      }

      if ( row_entry_count % Align ) { row_entry_count += Align - row_entry_count % Align ; }

      m_max_row_width = std::max( m_max_row_width , row_entry_count );

      entry_count += row_entry_count ;
    }

    m_entry_offset = array_unsigned_type( "CrsProductTensorLegendre::entry_offset" , 2 * m_dimension + 1 );
    m_coordinate   = array_unsigned_type( "CrsProductTensorLegendre::coordinate" , entry_count );
    m_value        = array_scalar_type(   "CrsProductTensorLegendre::value" , entry_count );

    array_unsigned_host_type host_entry_offset = create_mirror_view( m_entry_offset );
    array_unsigned_host_type host_coordinate   = create_mirror_view( m_coordinate );
    array_scalar_host_type   host_value        = create_mirror_view( m_value );

    entry_count = 0 ;

    for ( unsigned i = 0 ; i < m_dimension ; ++i ) {

      // Diagonals first:

      host_entry_offset(2*i) = entry_count ;

      for ( unsigned j = 0 ; j < m_dimension ; ++j ) {
        if ( combinatorial.is_non_zero(i,j,j) ) {
          host_coordinate( entry_count ) = j ;
          host_value(      entry_count ) = combinatorial(i,j,j);
          ++entry_count ;
        }
      }

      if ( entry_count % Align ) { entry_count += Align - entry_count % Align ; }

      host_entry_offset(2*i+1) = entry_count ;

      for ( unsigned j = 0 ; j < m_dimension ; ++j ) {
        for ( unsigned k = j+1 ; k < m_dimension ; ++k ) {
          if ( combinatorial.is_non_zero(i,j,k) ) {
            host_coordinate( entry_count ) = ( k << 16 ) | j ;
            host_value(      entry_count ) = combinatorial(i,j,k);
            ++entry_count ;
          }
        }
      }

      if ( entry_count % Align ) { entry_count += Align - entry_count % Align ; }

      host_entry_offset(2*i+2) = entry_count ;
    }

    deep_copy( m_entry_offset , host_entry_offset );
    deep_copy( m_coordinate ,   host_coordinate );
    deep_copy( m_value ,        host_value );
  }

  /** \brief  This data structure's implementation of c += a * b */

  template< typename ScalarTypeA ,
            typename ScalarTypeB ,
            typename ScalarTypeC >
  KOKKOSARRAY_INLINE_FUNCTION
  void multiply_add( const ScalarTypeA a[] ,
                     const ScalarTypeB b[] ,
                           ScalarTypeC c[] ) const
  {
    for ( unsigned ic = 0 ; ic < m_dimension ; ++ic ) {
      const unsigned iBeg        = m_entry_offset(2*ic);
      const unsigned iBegOffDiag = m_entry_offset(2*ic+1);
      const unsigned iEnd        = m_entry_offset(2*ic+2);

      ScalarTypeC tmp = 0 ;

      for ( unsigned i = iBeg ; i < iBegOffDiag ; ++i ) {
        const unsigned j = m_coordinate(i);
        tmp += m_value(i) * a[j] * b[j] ;
      }

      for ( unsigned i = iBegOffDiag ; i < iEnd ; ++i ) {
        const unsigned kj = m_coordinate(i);
        const unsigned j  = kj & 0x0ffff ;
        const unsigned k  = kj >> 16 ;
        tmp += m_value(i) * ( a[j] * b[k] + a[k] * b[j] );
      }

      c[ic] += tmp ;
    }
  }
};

} // namespace KokkosArray

//----------------------------------------------------------------------------

#elif 1

namespace KokkosArray {

template< typename TensorScalar , class Device >
struct CrsProductTensorLegendre {
  typedef View< unsigned* ,     Device >  array_unsigned_type ;
  typedef View< TensorScalar* , Device >  array_scalar_type ;

  typedef typename array_unsigned_type::HostMirror array_unsigned_host_type ;
  typedef typename array_scalar_type  ::HostMirror array_scalar_host_type ;

  array_unsigned_type  m_entry_offset ;
  array_unsigned_type  m_coordinate ;
  array_scalar_type    m_value ;
  unsigned             m_dimension ;
  unsigned             m_max_row_width ;
  unsigned             m_nonzero_count ;
  unsigned             m_multiply_add_flops ;

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned dimension() const { return m_dimension ; }

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned max_row_width() const { return m_max_row_width ; }

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned multiply_add_flops() const { return m_multiply_add_flops ; }

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned nonzero_count() const { return m_nonzero_count ; }

  CrsProductTensorLegendre()
  : m_entry_offset()
  , m_coordinate()
  , m_value()
  , m_dimension(0)
  , m_max_row_width(0)
  , m_nonzero_count(0)
  , m_multiply_add_flops(0)
  {}

  CrsProductTensorLegendre( const CrsProductTensorLegendre & rhs )
  : m_entry_offset( rhs.m_entry_offset )
  , m_coordinate(   rhs.m_coordinate )
  , m_value(        rhs.m_value )
  , m_dimension(    rhs.m_dimension )
  , m_max_row_width(rhs.m_max_row_width )
  , m_nonzero_count(rhs.m_nonzero_count )
  , m_multiply_add_flops(rhs.m_multiply_add_flops)
  {}

  CrsProductTensorLegendre & operator = ( const CrsProductTensorLegendre & rhs )
  {
    m_entry_offset = rhs.m_entry_offset ;
    m_coordinate   = rhs.m_coordinate ;
    m_value        = rhs.m_value ;
    m_dimension    = rhs.m_dimension ;
    m_max_row_width= rhs.m_max_row_width ;
    m_nonzero_count= rhs.m_nonzero_count ;
    m_multiply_add_flops=rhs.m_multiply_add_flops;
    return *this ;
  }

  explicit
  CrsProductTensorLegendre( const std::vector<unsigned> & variable_poly_degree ,
                            const unsigned maximum_poly_degree = 2 )
  : m_entry_offset()
  , m_coordinate()
  , m_value()
  , m_dimension(0)
  , m_max_row_width(0)
  , m_nonzero_count(0)
  , m_multiply_add_flops(0)
  {
    enum { Blocking = 32 };
    enum { Align = Impl::is_same<Device,Cuda>::value ? 32 : 1 };

    const KokkosArray::TripleProductTensorLegendreCombinatorialEvaluation
      combinatorial( variable_poly_degree , maximum_poly_degree );

    m_dimension = combinatorial.bases_count();

    if ( ( 1 << 16 ) < m_dimension ) {
      throw std::runtime_error("CrsProductTensorLegendre tensor dimension too large");
    }

    unsigned entry_count = 0 ;

    // Adding to the output vector.
    m_multiply_add_flops = m_dimension ;

    for ( unsigned i = 0 ; i < m_dimension ; ++i ) {

      unsigned row_entry_count = 0 ;

      for ( unsigned j = 0 ; j < m_dimension ; ++j ) {
        for ( unsigned k = j ; k < m_dimension ; ++k ) {
          if ( combinatorial.is_non_zero(i,j,k) ) { ++row_entry_count ; m_multiply_add_flops += 5 ; }
        }
      }

      m_nonzero_count      += row_entry_count ;
      m_multiply_add_flops += 5 * row_entry_count ; // Two additions and three multiplies:

      if ( row_entry_count % Align ) { row_entry_count += Align - row_entry_count % Align ; }

      m_max_row_width = std::max( m_max_row_width , row_entry_count );

      entry_count += row_entry_count ;
    }

    m_entry_offset = array_unsigned_type( "CrsProductTensorLegendre::entry_offset" , 2 * m_dimension + 1 );
    m_coordinate   = array_unsigned_type( "CrsProductTensorLegendre::coordinate" , entry_count );
    m_value        = array_scalar_type(   "CrsProductTensorLegendre::value" , entry_count );

    array_unsigned_host_type host_entry_offset = create_mirror_view( m_entry_offset );
    array_unsigned_host_type host_coordinate   = create_mirror_view( m_coordinate );
    array_scalar_host_type   host_value        = create_mirror_view( m_value );

    entry_count = 0 ;

    for ( unsigned i = 0 ; i < m_dimension ; ++i ) {

      host_entry_offset(2*i)   = entry_count ;
      host_entry_offset(2*i+1) = entry_count ;

      for ( unsigned jBlock = 0 ; jBlock < m_dimension ; ) {

        const unsigned jEnd = std::min( jBlock + Blocking , m_dimension );

        for ( unsigned kBlock = jBlock ; kBlock < m_dimension ; ) {

          const unsigned kEnd = std::min( kBlock + Blocking , m_dimension );

          for ( unsigned j = jBlock ; j < jEnd ; ++j ) {
          for ( unsigned k = std::max( kBlock , j ) ; k < kEnd ; ++k ) {
            if ( combinatorial.is_non_zero(i,j,k) ) {
              // Diagonal term is treated as an off-diagonal
              //   c[i] += combinatorial(i,j,j) * 0.5 * ( a[j] * b[j] + a[j] * b[j] );
              host_coordinate( entry_count ) = ( k << 16 ) | j ;
              host_value(      entry_count ) = combinatorial(i,j,k) * ( j == k ? 0.5 : 1.0 );
              ++entry_count ;

              if ( host_value.dimension_0() < entry_count ) {
                throw std::runtime_error("CrsProductTensorLegendre Blocking error limit");
              }
            }
          }}

          kBlock = kEnd ;
        }
        jBlock = jEnd ;
      }

      if ( entry_count % Align ) { entry_count += Align - entry_count % Align ; }

      host_entry_offset(2*i+2) = entry_count ;
    }

    if ( host_value.dimension_0() != entry_count ) {
      throw std::runtime_error("CrsProductTensorLegendre Blocking error equal");
    }

    deep_copy( m_entry_offset , host_entry_offset );
    deep_copy( m_coordinate ,   host_coordinate );
    deep_copy( m_value ,        host_value );
  }

  /** \brief  This data structure's implementation of c += a * b */

  template< typename ScalarTypeA ,
            typename ScalarTypeB ,
            typename ScalarTypeC >
  KOKKOSARRAY_INLINE_FUNCTION
  void multiply_add( const ScalarTypeA a[] ,
                     const ScalarTypeB b[] ,
                           ScalarTypeC c[] ) const
  {
    for ( unsigned ic = 0 ; ic < m_dimension ; ++ic ) {
      const unsigned iBeg        = m_entry_offset(2*ic);
      const unsigned iBegOffDiag = m_entry_offset(2*ic+1);
      const unsigned iEnd        = m_entry_offset(2*ic+2);

      ScalarTypeC tmp = 0 ;

      for ( unsigned i = iBeg ; i < iBegOffDiag ; ++i ) {
        const unsigned j = m_coordinate(i);
        tmp += m_value(i) * a[j] * b[j] ;
      }

      for ( unsigned i = iBegOffDiag ; i < iEnd ; ++i ) {
        const unsigned kj = m_coordinate(i);
        const unsigned j  = kj & 0x0ffff ;
        const unsigned k  = kj >> 16 ;
        tmp += m_value(i) * ( a[j] * b[k] + a[k] * b[j] );
      }

      c[ic] += tmp ;
    }
  }
};

} // namespace KokkosArray

#endif

//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< typename MatrixScalar ,
          typename VectorScalar >
class Multiply<
  BlockCrsMatrix< CrsProductTensorLegendre< MatrixScalar , Host > ,
                  MatrixScalar , Host > ,
  View< VectorScalar** , LayoutLeft , Host > ,
  View< VectorScalar** , LayoutLeft , Host > >
{
private:

  typedef BlockCrsMatrix< CrsProductTensorLegendre< MatrixScalar , Host > ,
                          MatrixScalar , Host > matrix_type ;
  typedef View< VectorScalar** , LayoutLeft , Host > vector_type ;

  const matrix_type m_A ;
  const vector_type m_x ;
  const vector_type m_y ;

public:

  typedef Host                             device_type ;
  typedef typename device_type::size_type  size_type ;

  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const size_type iy ) const
  {
    // Compute spatial row 'iy'

    const size_type tensor_dim     = m_A.block.dimension();
    const size_type iBlockEntryEnd = m_A.graph.row_map[ iy + 1 ];
          size_type iBlockEntry    = m_A.graph.row_map[ iy ];

    VectorScalar * const y = & m_y( 0 , iy );

    for ( size_type iyInner = 0 ; iyInner < tensor_dim ; ++iyInner ) {
      y[iyInner] = 0 ;
    }

    for ( ; iBlockEntry < iBlockEntryEnd ; ++iBlockEntry ) {
      const VectorScalar * const x = & m_x( 0 , m_A.graph.entries( iBlockEntry ) );
      const MatrixScalar * const A = & m_A.values( 0 , iBlockEntry );

      m_A.block.multiply_add( A , x , y );
    }
  }

  Multiply( const matrix_type & arg_A ,
            const vector_type & arg_x ,
            const vector_type & arg_y )
  : m_A( arg_A )
  , m_x( arg_x )
  , m_y( arg_y )
  {}

  void run() const
  {
    parallel_for( m_y.dimension_1() , *this );
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

#endif /* #ifndef KOKKOSARRAY_CRSPRODUCTTENSORLEGENDRE_HPP */


