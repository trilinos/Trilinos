// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_CRSMATRIX_HPP
#define STOKHOS_CRSMATRIX_HPP

#include "KokkosArray_View.hpp"
#include "KokkosArray_CrsArray.hpp"

#include "Stokhos_Multiply.hpp"

namespace Stokhos {

/** \brief  CRS matrix.  */

template< typename ValueType , class Device >
class CrsMatrix {
public:
  typedef Device     device_type ;
  typedef ValueType  value_type ;

  KokkosArray::View< value_type[] , device_type >   values ;
  KokkosArray::CrsArray< int , device_type , void , int >  graph ;
};

template< typename MatrixValueType ,
          typename VectorValueType ,
          class Device >
void multiply( const CrsMatrix<MatrixValueType,Device> & A ,
               const KokkosArray::View<VectorValueType[],Device>         & x ,
               const KokkosArray::View<VectorValueType[],Device>         & y )
{
  multiply(A, x, y, DefaultSparseMatOps() );
}

template< typename MatrixValueType ,
          typename VectorValueType ,
          class Device ,
	  class SparseMatOps >
void multiply( const CrsMatrix<MatrixValueType,Device> & A ,
               const KokkosArray::View<VectorValueType[],Device>         & x ,
               const KokkosArray::View<VectorValueType[],Device>         & y ,
	       const SparseMatOps& smo = SparseMatOps() )
{
  typedef CrsMatrix<MatrixValueType,Device>  matrix_type ;
  typedef KokkosArray::View<VectorValueType[],Device>     vector_type ;

  Multiply<matrix_type,vector_type,vector_type,SparseMatOps>::apply( A , x , y );
}

template< typename MatrixValueType ,
          typename VectorValueType ,
	  typename OrdinalType ,
          class Device >
void multiply( const CrsMatrix<MatrixValueType,Device> & A ,
	       const KokkosArray::View<VectorValueType**, KokkosArray::LayoutLeft, Device> & x ,
	       const KokkosArray::View<VectorValueType**, KokkosArray::LayoutLeft, Device> & y ,
	       const std::vector<OrdinalType>& col_indices, 
	       bool use_block_multiply = true )
{
  multiply(A, x, y, col_indices, use_block_multiply, DefaultSparseMatOps() );
}

template< typename MatrixValueType ,
          typename VectorValueType ,
	  typename OrdinalType ,
          class Device ,
	  class SparseMatOps >
void multiply( const CrsMatrix<MatrixValueType,Device> & A ,
	       const KokkosArray::View<VectorValueType**, KokkosArray::LayoutLeft, Device> & x ,
	       const KokkosArray::View<VectorValueType**, KokkosArray::LayoutLeft, Device> & y ,
	       const std::vector<OrdinalType>& col_indices, 
	       bool use_block_multiply = true,
	       const SparseMatOps& smo = SparseMatOps() )
{
  typedef CrsMatrix<MatrixValueType,Device>           matrix_type ;
  typedef KokkosArray::View<VectorValueType[],Device>  vector_type ;
  typedef KokkosArray::View<VectorValueType**, KokkosArray::LayoutLeft, Device> multi_vector_type ;

  if (use_block_multiply)
    MMultiply<matrix_type,multi_vector_type,multi_vector_type,SparseMatOps>::apply( 
      A , x , y , col_indices );
  else {
    for (size_t i=0; i<col_indices.size(); ++i) {
      const vector_type x_view( x , col_indices[i] );
      const vector_type y_view( y , col_indices[i] );
      Multiply<matrix_type,vector_type,vector_type,SparseMatOps>::apply( 
	A , x_view , y_view );
    }
  }
}

template< typename MatrixValueType ,
          typename VectorValueType ,
	  class Device >
void multiply( const CrsMatrix<MatrixValueType,Device> & A ,
	       const std::vector< KokkosArray::View<VectorValueType[], Device> > & x ,
	       const std::vector< KokkosArray::View<VectorValueType[], Device> > & y ,
	       bool use_block_multiply = true )
{
  multiply(A, x, y, use_block_multiply, DefaultSparseMatOps() );
}

template< typename MatrixValueType ,
          typename VectorValueType ,
	  class Device ,
	  class SparseMatOps >
void multiply( const CrsMatrix<MatrixValueType,Device> & A ,
	       const std::vector< KokkosArray::View<VectorValueType[], Device> > & x ,
	       const std::vector< KokkosArray::View<VectorValueType[], Device> > & y ,
	       bool use_block_multiply = true,
	       const SparseMatOps& smo = SparseMatOps() )
{
  typedef CrsMatrix<MatrixValueType,Device>           matrix_type ;
  typedef KokkosArray::View<VectorValueType[],Device> vector_type ;

  if (use_block_multiply)
    MMultiply<matrix_type,vector_type,vector_type,SparseMatOps>::apply( 
      A , x , y  );
  else {
    for (size_t i=0; i<x.size(); ++i) {
      Multiply<matrix_type,vector_type,vector_type,SparseMatOps>::apply( 
	A , x[i] , y[i] );
    }
  }
}

template< typename MatrixValueType ,
          class Device >
void write_matrix_market(const CrsMatrix<MatrixValueType,Device> & A ,
			 const std::string& filename)
{
  MatrixMarketWriter<MatrixValueType,Device>::write(A, filename);
}

template< typename ValueType, 
	  typename VectorValueType ,
          class Device >
void update( const ValueType& alpha, 
	     const KokkosArray::View<VectorValueType[],Device> & x ,
	     const ValueType& beta,
	     const KokkosArray::View<VectorValueType[],Device> & y )
{
  typedef ValueType                          value_type ;
  typedef KokkosArray::View<VectorValueType[],Device>     vector_type ;

  Update<value_type,vector_type>::apply( alpha , x , beta, y );
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Stokhos

#endif /* #ifndef KOKKOSARRAY_CRSMATRIX_HPP */

