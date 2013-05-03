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

#ifndef KOKKOSARRAY_CRSMATRIX_HPP
#define KOKKOSARRAY_CRSMATRIX_HPP

#if defined( __CUDACC__ )
#include <cusparse.h>
#endif

#include <iostream>
#include <stdexcept>
#include <KokkosArray_CrsArray.hpp>
#include <KokkosArray_Array.hpp>
#include <impl/KokkosArray_ArrayAnalyzeShape.hpp>
#include <impl/KokkosArray_ArrayViewDefault.hpp>

namespace KokkosArray {
namespace Impl {
template< class MatrixType , class InputVector , class OutputVector >
class Multiply ;
}
}

namespace KokkosArray {

/** \brief  CRS matrix.  */

template< typename ValueType , class Device >
class CrsMatrix {
public:
  typedef Device     device_type ;
  typedef ValueType  value_type ;

  View< value_type * , LayoutRight , device_type >   values ;
  CrsArray< int , device_type , device_type , int >  graph ;
};


template< typename ScalarType , unsigned N , class Device >
class CrsMatrix< Array<ScalarType,N> , Device > {
public:
  typedef Device               device_type ;
  typedef Array<ScalarType,N>  value_type ;

  View< value_type * , LayoutRight , device_type >   values ;
  CrsArray< int , device_type , device_type , int >  graph ;
};


template< typename MatrixType ,
          typename InputVectorType ,
          typename OutputVectorType >
void multiply( const MatrixType       & A ,
               const InputVectorType  & x ,
               const OutputVectorType & y )
{
  Impl::Multiply<MatrixType,InputVectorType,OutputVectorType>( A , x , y );
}

}

namespace KokkosArray {
namespace Impl {

template< class MatrixValueType ,
          class InputValueType ,
          class OutputValueType ,
          class Device >
class Multiply< CrsMatrix< MatrixValueType , Device > ,
                View<      InputValueType * , LayoutRight , Device > ,
                View<      OutputValueType* , LayoutRight , Device > >
{
public:

  typedef Device device_type ;

  typedef  CrsMatrix< MatrixValueType , Device > matrix_type ;
  typedef  View<      InputValueType * , LayoutRight , Device > input_type ;
  typedef  View<      OutputValueType* , LayoutRight , Device > output_type ;

  const matrix_type  m_A ;
  const input_type   m_x ;
  const output_type  m_y ;

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const iType iRow ) const
  {
    const int iEntryEnd = m_A.graph.row_map(iRow+1);

    OutputValueType y( 0 );

    for ( int iEntry = m_A.graph.row_map(iRow) ; iEntry < iEntryEnd ; ++iEntry ) {
      y += m_A.values(iEntry) * m_x( m_A.graph.entries(iEntry) );
    }

    m_y(iRow) = y ;
  }

  Multiply( const matrix_type  & A ,
            const input_type   & x ,
            const output_type  & y )
    : m_A( A ), m_x( x ), m_y( y )
  {
    const int nRow = m_A.graph.row_map.dimension_0() - 1 ;
    parallel_for( nRow , *this );
  }
};

template< class MatrixValueType ,
          class InputValueType ,
          class OutputValueType ,
          unsigned N ,
          class Device >
class Multiply< CrsMatrix< Array<MatrixValueType,N> , Device > ,
                View<      Array<InputValueType,N> * , LayoutRight , Device > ,
                View<      Array<OutputValueType,N>* , LayoutRight , Device > >
{
public:

  typedef Device device_type ;

  typedef  CrsMatrix< Array<MatrixValueType,N> , Device > matrix_type ;
  typedef  View<      Array<InputValueType,N> * , LayoutRight , Device > input_type ;
  typedef  View<      Array<OutputValueType,N>* , LayoutRight , Device > output_type ;

  const matrix_type  m_A ;
  const input_type   m_x ;
  const output_type  m_y ;

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const iType iRow ) const
  {
    const int iEntryEnd = m_A.graph.row_map(iRow+1);

    Array<OutputValueType,N> y( 0 );

    for ( int iEntry = m_A.graph.row_map(iRow) ; iEntry < iEntryEnd ; ++iEntry ) {
      y += m_A.values(iEntry) * m_x( m_A.graph.entries(iEntry) );
    }

    m_y(iRow) = y ;
  }

  Multiply( const matrix_type  & A ,
            const input_type   & x ,
            const output_type  & y )
    : m_A( A ), m_x( x ), m_y( y )
  {
    const int nRow = m_A.graph.row_map.dimension_0() - 1 ;
    parallel_for( nRow , *this );
  }
};

//----------------------------------------------------------------------------

#if defined( __CUDACC__ )

template< class MatrixValueType ,
          class InputValueType ,
          class OutputValueType ,
          unsigned N >
class Multiply< CrsMatrix< Array<MatrixValueType,N> , Cuda > ,
                View<  Array<InputValueType,N> * , LayoutRight , Cuda > ,
                View<  Array<OutputValueType,N> *, LayoutRight , Cuda > >
{
public:

  typedef Cuda device_type ;

  typedef CrsMatrix<  Array< MatrixValueType,N> ,   device_type > matrix_type ;
  typedef View<       Array< InputValueType, N> * , LayoutRight , device_type > input_type ;
  typedef View<       Array< OutputValueType,N> * , LayoutRight , device_type > output_type ;

  const matrix_type m_A ;
  const input_type  m_x ;
  const output_type m_y ;

  __device__
  void operator()(void) const
  {
    const int iRow = threadIdx.y + blockDim.y * blockIdx.x ;

    if ( threadIdx.x < N && iRow < m_y.dimension_0() ) {
      const int iEntryEnd = m_A.graph.row_map(iRow+1);

      OutputValueType y = 0 ;

      for ( int iEntry = m_A.graph.row_map(iRow) ; iEntry < iEntryEnd ; ++iEntry ) {
        y += m_A.values(iEntry,threadIdx.x) * m_x( m_A.graph.entries(iEntry) , threadIdx.x );
      }

      m_y(iRow,threadIdx.x) = y ;
    }
  }

  Multiply( const matrix_type & A ,
            const input_type  & x ,
            const output_type & y )
    : m_A( A ), m_x( x ), m_y( y )
  {
    enum { W = CudaTraits::WarpSize };
    enum { NX = ( N + W - 1 ) / W };
    enum { NY = NX < 4 ? 16 / NX : 1 };
    const int nRow = m_A.graph.row_map.dimension_0() - 1 ;

    const dim3 dBlock( W * NX , NY , 1 );
    const dim3 dGrid( ( nRow + NY - 1 ) / NY , 1 , 1 );

    cuda_parallel_launch_local_memory< Multiply ><<< dGrid , dBlock >>>( *this );
  }
};

//----------------------------------------------------------------------------

struct CudaSparseSingleton {
  cusparseHandle_t   handle;
  cusparseMatDescr_t descra;

  CudaSparseSingleton()
  {
    cusparseCreate( & handle );
    cusparseCreateMatDescr( & descra );
    cusparseSetMatType(       descra , CUSPARSE_MATRIX_TYPE_GENERAL );
    cusparseSetMatIndexBase(  descra , CUSPARSE_INDEX_BASE_ZERO );
  }

  static CudaSparseSingleton & singleton();

};

CudaSparseSingleton & CudaSparseSingleton::singleton()
{ static CudaSparseSingleton s ; return s ; }

template<>
class Multiply< CrsMatrix< double , Cuda > ,
                View<      double *, LayoutRight , Cuda > ,
                View<      double *, LayoutRight , Cuda > >
{
public:

  typedef Cuda device_type ;

  typedef CrsMatrix<  double ,   device_type > matrix_type ;
  typedef View<       double * , LayoutRight , device_type > input_type ;
  typedef View<       double * , LayoutRight , device_type > output_type ;

  Multiply( const matrix_type & A ,
            const input_type  & x ,
            const output_type & y )
  {
    const int n = y.dimension_0();

    CudaSparseSingleton & s = CudaSparseSingleton::singleton();

    const double alpha = 1 , beta = 0 ;

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

#endif  /* #if defined( __CUDACC__ ) */

//----------------------------------------------------------------------------

}
}

#endif

