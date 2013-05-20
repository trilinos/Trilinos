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

#ifndef KOKKOSARRAY_SPARSEPRODUCTTENSORLEGENDRE_HPP
#define KOKKOSARRAY_SPARSEPRODUCTTENSORLEGENDRE_HPP

#include <cmath>
#include <utility>
#include <vector>
#include <stdexcept>
#include <KokkosArray_Macros.hpp>
#include <KokkosArray_View.hpp>
#include <KokkosArray_BlockCrsMatrix.hpp>
#include <KokkosArray_ProductTensorLegendre.hpp>

//----------------------------------------------------------------------------

namespace KokkosArray {

enum SparseProductTensorLegendreVariant {
  SparseProductTensorLegendreVariant_Default
};

template< typename TensorScalar , class Device ,
	  SparseProductTensorLegendreVariant Variant = SparseProductTensorLegendreVariant_Default >
struct SparseProductTensorLegendre ;

} // namespace KokkosArray


//----------------------------------------------------------------------------

namespace KokkosArray {

template< typename TensorScalar , class Device , SparseProductTensorLegendreVariant Variant >
struct SparseProductTensorLegendre
{
  enum { unsigned_long_is_8_bytes = Impl::StaticAssert< sizeof(unsigned long) == 8 >::value };

  typedef Device  device_type ;

  typedef View< TensorScalar  * , device_type >  array_scalar_type ;
  typedef View< unsigned long * , device_type >  array_coord_type ;

  array_scalar_type  m_value ;
  array_coord_type   m_coordinate ;
  unsigned           m_dimension ;
  unsigned           m_multiply_add_flops ;


  KOKKOSARRAY_INLINE_FUNCTION
  unsigned entry_count() const { return m_value.dimension_0(); }

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned dimension() const { return m_dimension ; }

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned multiply_add_flops() const { return m_multiply_add_flops ; }

  SparseProductTensorLegendre()
  : m_value()
  , m_coordinate()
  , m_dimension(0)
  , m_multiply_add_flops(0)
  {}

  SparseProductTensorLegendre( const SparseProductTensorLegendre & rhs )
  : m_value(        rhs.m_value )
  , m_coordinate(   rhs.m_coordinate )
  , m_dimension(    rhs.m_dimension )
  , m_multiply_add_flops(rhs.m_multiply_add_flops)
  {}

  SparseProductTensorLegendre & operator = ( const SparseProductTensorLegendre & rhs )
  {
    m_value        = rhs.m_value ;
    m_coordinate   = rhs.m_coordinate ;
    m_dimension    = rhs.m_dimension ;
    m_multiply_add_flops = rhs.m_multiply_add_flops;
    return *this ;
  }

  explicit
  SparseProductTensorLegendre( const std::vector<unsigned> & variable_poly_degree ,
                               const unsigned maximum_poly_degree = 2 )
  : m_value()
  , m_coordinate()
  , m_dimension(0)
  , m_multiply_add_flops(0)
  {
    const KokkosArray::TripleProductTensorLegendreCombinatorialEvaluation
      combinatorial( variable_poly_degree , maximum_poly_degree );

    m_dimension = combinatorial.bases_count();

    unsigned entry_count = 0 ;

    for ( unsigned i = 0 ; i < m_dimension ; ++i ) {
    for ( unsigned j = 0 ; j <= i ; ++j ) {
    for ( unsigned k = 0 ; k <= j ; ++k ) {
      if ( combinatorial.is_non_zero(i,j,k) ) ++entry_count ;
    }}}

    m_coordinate  = array_coord_type(  "SparseProductTensorLegendre::coordinate" , entry_count );
    m_value       = array_scalar_type( "SparseProductTensorLegendre::value" , entry_count );

    typename array_scalar_type::HostMirror host_value = create_mirror_view( m_value );
    typename array_coord_type ::HostMirror host_coord = create_mirror_view( m_coordinate );

    entry_count = 0 ;

    for ( unsigned long i = 0 ; i < m_dimension ; ++i ) {
    for ( unsigned long j = 0 ; j <= i ; ++j ) {
    for ( unsigned long k = 0 ; k <= j ; ++k ) {
      if ( combinatorial.is_non_zero(i,j,k) ) {
        host_value( entry_count ) = combinatorial(i,j,k);
        host_coord( entry_count ) = ( i << 32 ) | ( j << 16 ) | ( k );
        ++entry_count ;

        // Guaranteed by construction that: i >= j >= k
        // When: i != j  then  i >  j >= k  and therefore  i != k
        // When: j != k  then  i >= j >  k  and therefore  i != k

        // Counting 'usefull' flops vs. logic masking flops.
                      m_multiply_add_flops += 5 ;
        if ( i != j ) m_multiply_add_flops += 5 ;
        if ( j != k ) m_multiply_add_flops += 5 ;
      }
    }}}

    deep_copy( m_coordinate , host_coord );
    deep_copy( m_value ,      host_value );
  }

  /** \brief  This data structure's implementation of c += a * b ; */

  template< typename ScalarTypeA ,
            typename ScalarTypeB ,
            typename ScalarTypeC >
  inline
  void multiply_add( const ScalarTypeA a[] ,
                     const ScalarTypeB b[] ,
                           ScalarTypeC c[] ) const
  {
    for ( unsigned ic = 0 ; ic < m_value.dimension_0() ; ++ic ) {

      const TensorScalar  v = m_value(ic);
      const unsigned long coord = m_coordinate(ic);
      const int i = ( coord >> 32 ) & 0x0ffff ;
      const int j = ( coord >> 16 ) & 0x0ffff ;
      const int k = ( coord       ) & 0x0ffff ;
      const int neq_ij = i != j ; // 0 or 1
      const int neq_jk = j != k ; // 0 or 1

      // Guaranteed by construction that: i >= j >= k
      // When: i != j  then  i >  j >= k  and therefore  i != k
      // When: j != k  then  i >= j >  k  and therefore  i != k

                    c[i] += v * ( a[j] * b[k] + a[k] * b[j] * neq_jk );
      if ( neq_ij ) c[j] += v * ( a[k] * b[i] + a[i] * b[k] );
      if ( neq_jk ) c[k] += v * ( a[i] * b[j] + a[j] * b[i] * neq_ij );
    }
  }
};

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< typename MatrixScalar ,
          typename VectorScalar >
class Multiply<
  BlockCrsMatrix<
    SparseProductTensorLegendre< MatrixScalar , Host ,
                                 SparseProductTensorLegendreVariant_Default > ,
    MatrixScalar , Host > ,
  View< VectorScalar** , LayoutLeft , Host > ,
  View< VectorScalar** , LayoutLeft , Host > >
{
private:

  typedef BlockCrsMatrix<
            SparseProductTensorLegendre< MatrixScalar , Host ,
                                         SparseProductTensorLegendreVariant_Default > ,
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

#endif /* #ifndef KOKKOSARRAY_SPARSEPRODUCTTENSORLEGENDRE_HPP */


