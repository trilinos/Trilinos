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

#ifndef STOKHOS_CUDA_CRSMATRIX_HPP
#define STOKHOS_CUDA_CRSMATRIX_HPP

#include <utility>
#include <sstream>
#include <stdexcept>

#include <cuda_runtime.h>
#include <cusparse.h>

#include "KokkosArray_Cuda.hpp"
#include "Cuda/KokkosArray_Cuda_Parallel.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_CrsMatrix.hpp"

namespace Stokhos {

class CudaSparseSingleton {
public:

  cusparseStatus_t   status;
  cusparseHandle_t   handle;
  cusparseMatDescr_t descra;

  static CudaSparseSingleton & singleton();

private:

  CudaSparseSingleton()
  {
    status = cusparseCreate(&handle);
    if(status != CUSPARSE_STATUS_SUCCESS)
    {
      throw std::runtime_error( std::string("ERROR - CUSPARSE Library Initialization failed" ) );
    }

    status = cusparseCreateMatDescr(&descra);
    if(status != CUSPARSE_STATUS_SUCCESS)
    {
      throw std::runtime_error( std::string("ERROR - CUSPARSE Library Matrix descriptor failed" ) );
    }

    cusparseSetMatType(descra , CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descra , CUSPARSE_INDEX_BASE_ZERO);
  }

  CudaSparseSingleton( const CudaSparseSingleton & );
  CudaSparseSingleton & operator = ( const CudaSparseSingleton & );
};

CudaSparseSingleton & CudaSparseSingleton::singleton()
{
  static CudaSparseSingleton s ; return s ;
}


template<>
class Multiply<
  CrsMatrix< float , KokkosArray::Cuda > ,
  KokkosArray::View< float[] , KokkosArray::Cuda > ,
  KokkosArray::View< float[] , KokkosArray::Cuda > >
{
public:
  typedef KokkosArray::Cuda                        device_type ;
  typedef device_type::size_type              size_type ;
  typedef KokkosArray::View< float[] , device_type >  vector_type ;
  typedef CrsMatrix< float , device_type >    matrix_type ;

  //--------------------------------------------------------------------------

  static void apply( const matrix_type & A ,
                     const vector_type & x ,
                     const vector_type & y )
  {
    CudaSparseSingleton & s = CudaSparseSingleton::singleton();
    const float alpha = 1 , beta = 0 ;
    const int n = A.graph.row_map.dimension(0) - 1 ;
    // const int nz = A.graph.entry_count();

    cusparseStatus_t status =
      cusparseScsrmv( s.handle ,
                      CUSPARSE_OPERATION_NON_TRANSPOSE ,
                      n , n ,
                      alpha ,
                      s.descra ,
                      A.values.ptr_on_device() ,
                      A.graph.row_map.ptr_on_device() ,
                      A.graph.entries.ptr_on_device() ,
                      x.ptr_on_device() , 
                      beta ,
                      y.ptr_on_device() );

    if ( CUSPARSE_STATUS_SUCCESS != status ) {
      throw std::runtime_error( std::string("ERROR - cusparseScsrmv " ) );
    }
  }
};

template<>
class Multiply<
  CrsMatrix< double , KokkosArray::Cuda > ,
  KokkosArray::View< double[] , KokkosArray::Cuda > ,
  KokkosArray::View< double[] , KokkosArray::Cuda > >
{
public:
  typedef KokkosArray::Cuda                         device_type ;
  typedef device_type::size_type               size_type ;
  typedef KokkosArray::View< double[] , device_type >  vector_type ;
  typedef CrsMatrix< double , device_type >    matrix_type ;

  //--------------------------------------------------------------------------

  static void apply( const matrix_type & A ,
                     const vector_type & x ,
                     const vector_type & y )
  {
    CudaSparseSingleton & s = CudaSparseSingleton::singleton();
    const double alpha = 1 , beta = 0 ;
    const int n = A.graph.row_map.dimension(0) - 1 ;
    // const int nz = A.graph.entry_count();

    cusparseStatus_t status =
      cusparseDcsrmv( s.handle ,
                      CUSPARSE_OPERATION_NON_TRANSPOSE ,
                      n , n ,
                      alpha ,
                      s.descra ,
                      A.values.ptr_on_device() ,
                      A.graph.row_map.ptr_on_device() ,
                      A.graph.entries.ptr_on_device() ,
                      x.ptr_on_device() , 
                      beta ,
                      y.ptr_on_device() );

    if ( CUSPARSE_STATUS_SUCCESS != status ) {
      throw std::runtime_error( std::string("ERROR - cusparseDcsrmv " ) );
    }
  }
};

template<>
class MMultiply<
  CrsMatrix< float , KokkosArray::Cuda > ,
  KokkosArray::View< float** , KokkosArray::LayoutLeft, KokkosArray::Cuda > ,
  KokkosArray::View< float** , KokkosArray::LayoutLeft, KokkosArray::Cuda > >
{
public:
  typedef KokkosArray::Cuda                           device_type ;
  typedef device_type::size_type                      size_type ;
  typedef KokkosArray::View< float[] , device_type >              vector_type ;
  typedef KokkosArray::View< float** , KokkosArray::LayoutLeft, device_type >  multi_vector_type ;
  typedef CrsMatrix< float , device_type >           matrix_type ;
  typedef int                                         Ordinal ;

  //--------------------------------------------------------------------------

  static void apply( const matrix_type & A ,
                     const multi_vector_type & x ,
                     const multi_vector_type & y ,
		     const std::vector<Ordinal> & col_indices )
  {
    CudaSparseSingleton & s = CudaSparseSingleton::singleton();
    const float alpha = 1 , beta = 0 ;
    const int n = A.graph.row_map.dimension(0) - 1 ;
    // const int nz = A.graph.entry_count();
    const size_t ncol = col_indices.size();

    // Copy columns of x into a contiguous vector
    vector_type xx( "xx" , n * ncol );
    vector_type yy( "yy" , n * ncol );

    for (size_t col=0; col<ncol; col++) {
      const std::pair< size_t , size_t > span( n * col , n * ( col + 1 ) );
      const vector_type xx_view( xx , span );
      const vector_type x_col( x, col_indices[col] );
      KokkosArray::deep_copy(xx_view, x_col);
    }

    // Sparse matrix-times-multivector
    cusparseStatus_t status =
      cusparseScsrmm( s.handle ,
		      CUSPARSE_OPERATION_NON_TRANSPOSE ,
		      n , ncol , n , 
		      alpha ,
		      s.descra ,
		      A.values.ptr_on_device() ,
		      A.graph.row_map.ptr_on_device() ,
		      A.graph.entries.ptr_on_device() ,
		      xx.ptr_on_device() , 
		      n , 
		      beta ,
		      yy.ptr_on_device() ,
		      n );
    
    if ( CUSPARSE_STATUS_SUCCESS != status ) {
      throw std::runtime_error( std::string("ERROR - cusparseDcsrmv " ) );
    }
    
    // Copy columns out of continguous multivector
    for (size_t col=0; col<ncol; col++) {
      const std::pair< size_t , size_t > span( n * col , n * ( col + 1 ) );
      const vector_type yy_view( yy , span );
      const vector_type y_col( y, col_indices[col] );
      KokkosArray::deep_copy(y_col, yy_view );
    }
  }
};

template<>
class MMultiply<
  CrsMatrix< double , KokkosArray::Cuda > ,
  KokkosArray::View< double** , KokkosArray::LayoutLeft, KokkosArray::Cuda > ,
  KokkosArray::View< double** , KokkosArray::LayoutLeft, KokkosArray::Cuda > >
{
public:
  typedef KokkosArray::Cuda                           device_type ;
  typedef device_type::size_type                      size_type ;
  typedef KokkosArray::View< double[] , device_type >              vector_type ;
  typedef KokkosArray::View< double** , KokkosArray::LayoutLeft, device_type >  multi_vector_type ;
  typedef CrsMatrix< double , device_type >           matrix_type ;
  typedef int                                         Ordinal ;

  //--------------------------------------------------------------------------

  static void apply( const matrix_type & A ,
                     const multi_vector_type & x ,
                     const multi_vector_type & y ,
		     const std::vector<Ordinal> & col_indices )
  {
    CudaSparseSingleton & s = CudaSparseSingleton::singleton();
    const double alpha = 1 , beta = 0 ;
    const int n = A.graph.row_map.dimension(0) - 1 ;
    // const int nz = A.graph.entry_count();
    const size_t ncol = col_indices.size();

    // Copy columns of x into a contiguous vector
    vector_type xx( "xx" , n * ncol );
    vector_type yy( "yy" , n * ncol );

    for (size_t col=0; col<ncol; col++) {
      const std::pair< size_t , size_t > span( n * col , n * ( col + 1 ) );
      const vector_type xx_view( xx , span );
      const vector_type x_col( x, col_indices[col] );
      KokkosArray::deep_copy(xx_view, x_col);
    }

    // Sparse matrix-times-multivector
    cusparseStatus_t status =
      cusparseDcsrmm( s.handle ,
		      CUSPARSE_OPERATION_NON_TRANSPOSE ,
		      n , ncol , n , 
		      alpha ,
		      s.descra ,
		      A.values.ptr_on_device() ,
		      A.graph.row_map.ptr_on_device() ,
		      A.graph.entries.ptr_on_device() ,
		      xx.ptr_on_device() , 
		      n , 
		      beta ,
		      yy.ptr_on_device() ,
		      n );
    
    if ( CUSPARSE_STATUS_SUCCESS != status ) {
      throw std::runtime_error( std::string("ERROR - cusparseDcsrmv " ) );
    }
    
    // Copy columns out of continguous multivector
    for (size_t col=0; col<ncol; col++) {
      const std::pair< size_t , size_t > span( n * col , n * ( col + 1 ) );
      const vector_type yy_view( yy , span );
      const vector_type y_col( y, col_indices[col] );
      KokkosArray::deep_copy(y_col, yy_view );
    }
  }
};

template<>
class MMultiply<
  CrsMatrix< float , KokkosArray::Cuda > ,
  KokkosArray::View< float[] , KokkosArray::Cuda > ,
  KokkosArray::View< float[] , KokkosArray::Cuda > >
{
public:
  typedef KokkosArray::Cuda                         device_type ;
  typedef device_type::size_type               size_type ;
  typedef KokkosArray::View< float[] , device_type >  vector_type ;
  typedef CrsMatrix< float , device_type >    matrix_type ;

  //--------------------------------------------------------------------------

  static void apply( const matrix_type & A ,
                     const std::vector<vector_type> & x ,
                     const std::vector<vector_type> & y )
  {
    CudaSparseSingleton & s = CudaSparseSingleton::singleton();
    const float alpha = 1 , beta = 0 ;
    const int n = A.graph.row_map.dimension(0) - 1 ;
    // const int nz = A.graph.entry_count();
    const size_t ncol = x.size();

    // Copy columns of x into a contiguous vector
    vector_type xx( "xx" , n * ncol );
    vector_type yy( "yy" , n * ncol );

    for (size_t col=0; col<ncol; col++) {
      const std::pair< size_t , size_t > span( n * col , n * ( col + 1 ) );
      const vector_type xx_view( xx , span );
      KokkosArray::deep_copy(xx_view, x[col]);
    }

    // Sparse matrix-times-multivector
    cusparseStatus_t status =
      cusparseScsrmm( s.handle ,
		      CUSPARSE_OPERATION_NON_TRANSPOSE ,
		      n , ncol , n , 
		      alpha ,
		      s.descra ,
		      A.values.ptr_on_device() ,
		      A.graph.row_map.ptr_on_device() ,
		      A.graph.entries.ptr_on_device() ,
		      xx.ptr_on_device() , 
		      n , 
		      beta ,
		      yy.ptr_on_device() ,
		      n );
    
    if ( CUSPARSE_STATUS_SUCCESS != status ) {
      throw std::runtime_error( std::string("ERROR - cusparseDcsrmv " ) );
    }
    
    // Copy columns out of continguous multivector
    for (size_t col=0; col<ncol; col++) {
      const std::pair< size_t , size_t > span( n * col , n * ( col + 1 ) );
      const vector_type yy_view( yy , span );
      KokkosArray::deep_copy(y[col], yy_view );
    }
  }
};

template<>
class MMultiply<
  CrsMatrix< double , KokkosArray::Cuda > ,
  KokkosArray::View< double[] , KokkosArray::Cuda > ,
  KokkosArray::View< double[] , KokkosArray::Cuda > >
{
public:
  typedef KokkosArray::Cuda                         device_type ;
  typedef device_type::size_type               size_type ;
  typedef KokkosArray::View< double[] , device_type >  vector_type ;
  typedef CrsMatrix< double , device_type >    matrix_type ;

  //--------------------------------------------------------------------------

  static void apply( const matrix_type & A ,
                     const std::vector<vector_type> & x ,
                     const std::vector<vector_type> & y )
  {
    CudaSparseSingleton & s = CudaSparseSingleton::singleton();
    const double alpha = 1 , beta = 0 ;
    const int n = A.graph.row_map.dimension(0) - 1 ;
    // const int nz = A.graph.entry_count();
    const size_t ncol = x.size();

    // Copy columns of x into a contiguous vector
    vector_type xx( "xx" , n * ncol );
    vector_type yy( "yy" , n * ncol );

    for (size_t col=0; col<ncol; col++) {
      const std::pair< size_t , size_t > span( n * col , n * ( col + 1 ) );
      const vector_type xx_view( xx , span );
      KokkosArray::deep_copy(xx_view, x[col]);
    }

    // Sparse matrix-times-multivector
    cusparseStatus_t status =
      cusparseDcsrmm( s.handle ,
		      CUSPARSE_OPERATION_NON_TRANSPOSE ,
		      n , ncol , n , 
		      alpha ,
		      s.descra ,
		      A.values.ptr_on_device() ,
		      A.graph.row_map.ptr_on_device() ,
		      A.graph.entries.ptr_on_device() ,
		      xx.ptr_on_device() , 
		      n , 
		      beta ,
		      yy.ptr_on_device() ,
		      n );
    
    if ( CUSPARSE_STATUS_SUCCESS != status ) {
      throw std::runtime_error( std::string("ERROR - cusparseDcsrmv " ) );
    }
    
    // Copy columns out of continguous multivector
    for (size_t col=0; col<ncol; col++) {
      const std::pair< size_t , size_t > span( n * col , n * ( col + 1 ) );
      const vector_type yy_view( yy , span );
      KokkosArray::deep_copy(y[col], yy_view );
    }
  }
};

template< typename MatrixValue>
class MatrixMarketWriter<MatrixValue,KokkosArray::Cuda>
{
public:
  typedef KokkosArray::Cuda                         device_type ;
  typedef device_type::size_type                    size_type ;
  typedef CrsMatrix< MatrixValue , device_type >    matrix_type ;

  MatrixMarketWriter() {}
  ~MatrixMarketWriter() {}

  static void write(const matrix_type & A ,
		    const std::string& filename) {}
};

//----------------------------------------------------------------------------

} // namespace Stokhos

#endif /* #ifndef STOKHOS_CUDA_CRSMATRIX_HPP */

