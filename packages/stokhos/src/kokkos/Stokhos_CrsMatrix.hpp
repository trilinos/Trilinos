// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
  KokkosArray::CrsArray< int , device_type , device_type , int >  graph ;
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

