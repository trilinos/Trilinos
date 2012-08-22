/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSARRAY_CUDA_CRSMATRIX_HPP
#define KOKKOSARRAY_CUDA_CRSMATRIX_HPP

#include <utility>
#include <sstream>
#include <stdexcept>

#include <cuda_runtime.h>
#include <cusparse.h>

namespace KokkosArray {
namespace Impl {

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
  typedef View< float[] , device_type >  vector_type ;
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
  typedef View< double[] , device_type >  vector_type ;
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
  CrsMatrix< double , KokkosArray::Cuda > ,
  KokkosArray::View< double[] , KokkosArray::Cuda > ,
  KokkosArray::View< double[] , KokkosArray::Cuda > >
{
public:
  typedef KokkosArray::Cuda                         device_type ;
  typedef device_type::size_type               size_type ;
  typedef View< double[] , device_type >  vector_type ;
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
  typedef Host                                      device_type ;
  typedef device_type::size_type                    size_type ;
  typedef CrsMatrix< MatrixValue , device_type >    matrix_type ;

  MatrixMarketWriter() {}
  ~MatrixMarketWriter() {}

  static void write(const matrix_type & A ,
		    const std::string& filename) {}
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

#endif /* #ifndef KOKKOSARRAY_CUDA_CRSMATRIX_HPP */

