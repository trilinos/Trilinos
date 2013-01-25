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

#ifndef KOKKOSARRAY_HOST_SPARSEPRODUCTTENSOR_HPP
#define KOKKOSARRAY_HOST_SPARSEPRODUCTTENSOR_HPP

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------

template< typename ValueType >
class Multiply< SparseProductTensor< 3 , ValueType , Host > , void , void >
{
public:

  typedef Host::size_type size_type ;
  typedef SparseProductTensor< 3 , ValueType , Host > tensor_type ;

  template< typename MatrixValue , typename VectorValue >
  static void apply( const tensor_type & tensor ,
                     const MatrixValue * const a ,
                     const VectorValue * const x ,
                           VectorValue * const y )
  {
    const size_type nEntry = tensor.entry_count();

    for ( size_type iEntry = 0 ; iEntry < nEntry ; ++iEntry ) {
      const size_type i = tensor.coord(iEntry,0);
      const size_type j = tensor.coord(iEntry,1);
      const size_type k = tensor.coord(iEntry,2);
      const ValueType v = tensor.value(iEntry);

      const bool neq_ij = i != j ;
      const bool neq_jk = j != k ;
      const bool neq_ki = k != i ;

      y[k] += neq_ij ? v * ( a[i] * x[j] + x[i] * a[j] )
                     : v * ( a[i] * x[i] );

      if ( neq_jk ) {
        y[j] += neq_ki ? v * ( a[i] * x[k] + x[i] * a[k] )
                       : v * ( a[i] * x[i] );
      }

      if ( neq_ki && neq_ij ) {
        y[i] += neq_jk ? v * ( a[k] * x[j] + x[k] * a[j] )
                       : v * ( a[j] * x[j] );
      }
    }
  }

  static size_type matrix_size( const tensor_type & tensor )
  { return tensor.dimension(); }

  static size_type vector_size( const tensor_type & tensor )
  { return tensor.dimension(); }
};

//----------------------------------------------------------------------------

template< typename ValueType >
class Multiply< CrsProductTensor< 3 , ValueType , Host > , void , void >
{
public:

  typedef Host::size_type size_type ;
  typedef CrsProductTensor< 3 , ValueType , Host > tensor_type ;

  template< typename MatrixValue , typename VectorValue >
  static void apply( const tensor_type & tensor ,
                     const MatrixValue * const a ,
                     const VectorValue * const x ,
                           VectorValue * const y )
  {
    const size_type nDim = tensor.dimension();

    for ( size_type iy = 0 ; iy < nDim ; ++iy ) {

      const size_type iEntryEnd = tensor.entry_end(iy);
            size_type iEntry    = tensor.entry_begin(iy);

      VectorValue ytmp = 0 ;

      for ( ; iEntry < iEntryEnd ; ++iEntry ) {

        const size_type j = tensor.coord(iEntry,0);
        const size_type k = tensor.coord(iEntry,1);

        // ytmp += tensor.value(iEntry) *
	//   ( j == k ? a[j] * x[j] : a[j] * x[k] + a[k] * x[j] );
	ytmp += tensor.value(iEntry) * ( a[j] * x[k] + a[k] * x[j] );
      }

      y[iy] += ytmp ;
    }
  }

  static size_type matrix_size( const tensor_type & tensor )
  { return tensor.dimension(); }

  static size_type vector_size( const tensor_type & tensor )
  { return tensor.dimension(); }
};

//----------------------------------------------------------------------------

template< typename ValueType >
class Multiply< SparseProductTensor< 3 , ValueType , Host > ,
                SymmetricDiagonalSpec< KokkosArray::Host > , void >
{
public:
  typedef Host::size_type size_type ;
  typedef SparseProductTensor< 3 , ValueType , Host > tensor_type ;

  template< typename MatrixValue >
  static void apply( const tensor_type & tensor ,
                     const MatrixValue * const a ,
                           MatrixValue * const M )
  {
    const SymmetricDiagonalSpec< KokkosArray::Host >
      spec( tensor.dimension() );

    const size_type nEntry  = tensor.entry_count();
    const size_type nMatrix = spec.matrix_size();

    for ( size_type iMatrix = 0 ; iMatrix < nMatrix ; ++iMatrix ) {
      M[ iMatrix ] = 0 ;
    }

    for ( size_type iEntry = 0 ; iEntry < nEntry ; ++iEntry ) {
      const size_type i = tensor.coord(iEntry,0);
      const size_type j = tensor.coord(iEntry,1);
      const size_type k = tensor.coord(iEntry,2);
      const ValueType v = tensor.value(iEntry);

      M[ spec.matrix_offset(k,j) ] += v * a[i] ;

      if ( i != j ) {
        M[ spec.matrix_offset(k,i) ] += v * a[j] ;
      }

      if ( i != k && j != k ) {
        M[ spec.matrix_offset(j,i) ] += v * a[k] ;
      }
    }
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

#endif /* #ifndef KOKKOSARRAY_HOST_SPARSEPRODUCTTENSOR_HPP */

