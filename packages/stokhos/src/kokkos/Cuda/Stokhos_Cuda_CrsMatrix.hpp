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

#ifndef STOKHOS_CUDA_CRSMATRIX_HPP
#define STOKHOS_CUDA_CRSMATRIX_HPP

#include "Stokhos_ConfigDefs.h"
#ifdef HAVE_STOKHOS_CUSPARSE

#include <utility>
#include <sstream>
#include <stdexcept>

#include <cuda_runtime.h>
//#include <cusparse.h>
#include <cusparse_v2.h>

#include "Kokkos_Core.hpp"

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
  CrsMatrix< float , Kokkos::Cuda > ,
  Kokkos::View< float* , Kokkos::Cuda > ,
  Kokkos::View< float* , Kokkos::Cuda > ,
  void,
  IntegralRank<1> >
{
public:
  typedef Kokkos::Cuda                        execution_space ;
  typedef execution_space::size_type              size_type ;
  typedef Kokkos::View< float* , execution_space >  vector_type ;
  typedef CrsMatrix< float , execution_space >    matrix_type ;

  //--------------------------------------------------------------------------

  static void apply( const matrix_type & A ,
                     const vector_type & x ,
                     const vector_type & y )
  {
    CudaSparseSingleton & s = CudaSparseSingleton::singleton();
    const float alpha = 1 , beta = 0 ;
    const int n = A.graph.row_map.extent(0) - 1 ;
    const int nz = A.graph.entries.extent(0);

    cusparseStatus_t status =
      cusparseScsrmv( s.handle ,
                      CUSPARSE_OPERATION_NON_TRANSPOSE ,
                      n , n , nz ,
                      &alpha ,
                      s.descra ,
                      A.values.data() ,
                      A.graph.row_map.data() ,
                      A.graph.entries.data() ,
                      x.data() ,
                      &beta ,
                      y.data() );

    if ( CUSPARSE_STATUS_SUCCESS != status ) {
      throw std::runtime_error( std::string("ERROR - cusparseScsrmv " ) );
    }
  }
};

template<>
class Multiply<
  CrsMatrix< double , Kokkos::Cuda > ,
  Kokkos::View< double* , Kokkos::Cuda > ,
  Kokkos::View< double* , Kokkos::Cuda > ,
  void,
  IntegralRank<1> >
{
public:
  typedef Kokkos::Cuda                         execution_space ;
  typedef execution_space::size_type               size_type ;
  typedef Kokkos::View< double* , execution_space >  vector_type ;
  typedef CrsMatrix< double , execution_space >    matrix_type ;

  //--------------------------------------------------------------------------

  static void apply( const matrix_type & A ,
                     const vector_type & x ,
                     const vector_type & y )
  {
    CudaSparseSingleton & s = CudaSparseSingleton::singleton();
    const double alpha = 1 , beta = 0 ;
    const int n = A.graph.row_map.extent(0) - 1 ;
    const int nz = A.graph.entries.extent(0);

    cusparseStatus_t status =
      cusparseDcsrmv( s.handle ,
                      CUSPARSE_OPERATION_NON_TRANSPOSE ,
                      n , n , nz ,
                      &alpha ,
                      s.descra ,
                      A.values.data() ,
                      A.graph.row_map.data() ,
                      A.graph.entries.data() ,
                      x.data() ,
                      &beta ,
                      y.data() );

    if ( CUSPARSE_STATUS_SUCCESS != status ) {
      throw std::runtime_error( std::string("ERROR - cusparseDcsrmv " ) );
    }
  }
};

template <typename Ordinal>
class Multiply<
  CrsMatrix< float , Kokkos::Cuda > ,
  Kokkos::View< float** , Kokkos::LayoutLeft, Kokkos::Cuda > ,
  Kokkos::View< float** , Kokkos::LayoutLeft, Kokkos::Cuda > ,
  std::vector<Ordinal> ,
  IntegralRank<2> >
{
public:
  typedef Kokkos::Cuda execution_space;
  typedef execution_space::size_type size_type;
  typedef Kokkos::View< float*, Kokkos::LayoutLeft, execution_space > vector_type;
  typedef Kokkos::View< float**, Kokkos::LayoutLeft, execution_space > multi_vector_type;
  typedef CrsMatrix< float , execution_space > matrix_type;

  //--------------------------------------------------------------------------

  static void apply( const matrix_type & A ,
                     const multi_vector_type & x ,
                     const multi_vector_type & y ,
                     const std::vector<Ordinal> & col_indices )
  {
    CudaSparseSingleton & s = CudaSparseSingleton::singleton();
    const float alpha = 1 , beta = 0 ;
    const int n = A.graph.row_map.extent(0) - 1 ;
    const int nz = A.graph.entries.extent(0);
    const size_t ncol = col_indices.size();

    // Copy columns of x into a contiguous vector
    vector_type xx( Kokkos::ViewAllocateWithoutInitializing("xx"), n * ncol );
    vector_type yy( Kokkos::ViewAllocateWithoutInitializing("yy"), n * ncol );

    for (size_t col=0; col<ncol; col++) {
      const std::pair< size_t , size_t > span( n * col , n * ( col + 1 ) );
      vector_type xx_view = Kokkos::subview( xx , span );
      vector_type x_col =
        Kokkos::subview( x, Kokkos::ALL(), col_indices[col] );
      Kokkos::deep_copy(xx_view, x_col);
    }

    // Sparse matrix-times-multivector
    cusparseStatus_t status =
      cusparseScsrmm( s.handle ,
                      CUSPARSE_OPERATION_NON_TRANSPOSE ,
                      n , ncol , n , nz ,
                      &alpha ,
                      s.descra ,
                      A.values.data() ,
                      A.graph.row_map.data() ,
                      A.graph.entries.data() ,
                      xx.data() ,
                      n ,
                      &beta ,
                      yy.data() ,
                      n );

    if ( CUSPARSE_STATUS_SUCCESS != status ) {
      throw std::runtime_error( std::string("ERROR - cusparseDcsrmv " ) );
    }

    // Copy columns out of continguous multivector
    for (size_t col=0; col<ncol; col++) {
      const std::pair< size_t , size_t > span( n * col , n * ( col + 1 ) );
      vector_type yy_view = Kokkos::subview( yy , span );
      vector_type y_col =
        Kokkos::subview( y, Kokkos::ALL(), col_indices[col] );
      Kokkos::deep_copy(y_col, yy_view );
    }
  }
};

#define USE_CUSPARSE 1
#if USE_CUSPARSE

template <typename Ordinal>
class Multiply<
  CrsMatrix< double , Kokkos::Cuda > ,
  Kokkos::View< double** , Kokkos::LayoutLeft, Kokkos::Cuda > ,
  Kokkos::View< double** , Kokkos::LayoutLeft, Kokkos::Cuda > ,
  std::vector<Ordinal> ,
  IntegralRank<2> >
{
public:
  typedef Kokkos::Cuda execution_space;
  typedef execution_space::size_type size_type;
  typedef Kokkos::View< double*, Kokkos::LayoutLeft, execution_space > vector_type;
  typedef Kokkos::View< double**, Kokkos::LayoutLeft, execution_space > multi_vector_type;
  typedef CrsMatrix< double , execution_space > matrix_type;

  //--------------------------------------------------------------------------

#define USE_TRANSPOSE 0
#if USE_TRANSPOSE

  // A version that copies the vectors to a transposed 2D view and calls
  // new CUSPARSE function for transpose layout.  Seems to be somewhat
  // slower????

  struct GatherTranspose {
    typedef Kokkos::Cuda execution_space;
    typedef execution_space::size_type size_type;

    multi_vector_type m_xt;
    const multi_vector_type m_x;
    const Kokkos::View<Ordinal*,execution_space> m_col;
    const size_type m_ncol;

    GatherTranspose( multi_vector_type& xt,
                     const multi_vector_type& x,
                     const Kokkos::View<Ordinal*,execution_space>& col ) :
      m_xt(xt), m_x(x), m_col(col), m_ncol(col.extent(0)) {}

    __device__
    inline void operator() (size_type i) const {
      for (size_type j=0; j<m_ncol; ++j)
        m_xt(j,i) = m_x(i,m_col(j));
    }

    static void apply( multi_vector_type& xt,
                       const multi_vector_type& x,
                       const Kokkos::View<Ordinal*,execution_space>& col ) {
      const size_type n = x.extent(0);
      Kokkos::parallel_for( n , GatherTranspose(xt,x,col) );
    }
  };

  static void apply( const matrix_type & A ,
                     const multi_vector_type & x ,
                     const multi_vector_type & y ,
                     const std::vector<Ordinal> & col_indices )
  {
    CudaSparseSingleton & s = CudaSparseSingleton::singleton();
    const double alpha = 1 , beta = 0 ;
    const int n = A.graph.row_map.extent(0) - 1 ;
    const int nz = A.graph.entries.extent(0);
    const size_t ncol = col_indices.size();

    // Copy col_indices to the device
    Kokkos::View<Ordinal*,execution_space> col_indices_dev(
      Kokkos::ViewAllocateWithoutInitializing("col_indices"), ncol);
    typename Kokkos::View<Ordinal*,execution_space>::HostMirror col_indices_host =
      Kokkos::create_mirror_view(col_indices_dev);
    for (size_t i=0; i<ncol; ++i)
      col_indices_host(i) = col_indices[i];
    Kokkos::deep_copy(col_indices_dev, col_indices_host);

    // Copy columns of x into a contiguous multi-vector and transpose
    multi_vector_type xx(
      Kokkos::ViewAllocateWithoutInitializing("xx"), ncol , n );
    GatherTranspose::apply(xx, x, col_indices_dev);

    // Temporary to store result (this is not transposed)
    multi_vector_type yy(
      Kokkos::ViewAllocateWithoutInitializing("yy"), n , ncol );

    // Sparse matrix-times-multivector
    cusparseStatus_t status =
      cusparseDcsrmm2( s.handle ,
                       CUSPARSE_OPERATION_NON_TRANSPOSE ,
                       CUSPARSE_OPERATION_TRANSPOSE ,
                       n , ncol , n , nz ,
                       &alpha ,
                       s.descra ,
                       A.values.data() ,
                       A.graph.row_map.data() ,
                       A.graph.entries.data() ,
                       xx.data() ,
                       ncol ,
                       &beta ,
                       yy.data() ,
                       n );

    if ( CUSPARSE_STATUS_SUCCESS != status ) {
      throw std::runtime_error( std::string("ERROR - cusparseDcsrmv " ) );
    }

    // Copy columns out of continguous multivector
    for (size_t col=0; col<ncol; col++) {
      vector_type yy_view =
        Kokkos::subview( yy ,  Kokkos::ALL(), col );
      vector_type y_col =
        Kokkos::subview( y, Kokkos::ALL(), col_indices[col] );
      Kokkos::deep_copy(y_col, yy_view );
    }
  }
#else
  static void apply( const matrix_type & A ,
                     const multi_vector_type & x ,
                     const multi_vector_type & y ,
                     const std::vector<Ordinal> & col_indices )
  {
    CudaSparseSingleton & s = CudaSparseSingleton::singleton();
    const double alpha = 1 , beta = 0 ;
    const int n = A.graph.row_map.extent(0) - 1 ;
    const int nz = A.graph.entries.extent(0);
    const size_t ncol = col_indices.size();

    // Copy columns of x into a contiguous vector
    vector_type xx( Kokkos::ViewAllocateWithoutInitializing("xx"), n * ncol );
    vector_type yy( Kokkos::ViewAllocateWithoutInitializing("yy"), n * ncol );

    for (size_t col=0; col<ncol; col++) {
      const std::pair< size_t , size_t > span( n * col , n * ( col + 1 ) );
      vector_type xx_view = Kokkos::subview( xx , span );
      vector_type x_col =
        Kokkos::subview( x, Kokkos::ALL(), col_indices[col] );
      Kokkos::deep_copy(xx_view, x_col);
    }

    // Sparse matrix-times-multivector
    cusparseStatus_t status =
      cusparseDcsrmm( s.handle ,
                      CUSPARSE_OPERATION_NON_TRANSPOSE ,
                      n , ncol , n , nz ,
                      &alpha ,
                      s.descra ,
                      A.values.data() ,
                      A.graph.row_map.data() ,
                      A.graph.entries.data() ,
                      xx.data() ,
                      n ,
                      &beta ,
                      yy.data() ,
                      n );

    if ( CUSPARSE_STATUS_SUCCESS != status ) {
      throw std::runtime_error( std::string("ERROR - cusparseDcsrmv " ) );
    }

    // Copy columns out of continguous multivector
    for (size_t col=0; col<ncol; col++) {
      const std::pair< size_t , size_t > span( n * col , n * ( col + 1 ) );
      vector_type yy_view = Kokkos::subview( yy , span );
      vector_type y_col =
        Kokkos::subview( y, Kokkos::ALL(), col_indices[col] );
      Kokkos::deep_copy(y_col, yy_view );
    }
  }
#endif
};

template <>
class Multiply<
  CrsMatrix< float , Kokkos::Cuda > ,
  Kokkos::View< float** , Kokkos::LayoutLeft, Kokkos::Cuda > ,
  Kokkos::View< float** , Kokkos::LayoutLeft, Kokkos::Cuda > ,
  void ,
  IntegralRank<2> >
{
public:
  typedef Kokkos::Cuda execution_space;
  typedef execution_space::size_type size_type;
  typedef Kokkos::View< float**, Kokkos::LayoutLeft, execution_space > multi_vector_type;
  typedef CrsMatrix< float , execution_space > matrix_type;

  //--------------------------------------------------------------------------

  static void apply( const matrix_type & A ,
                     const multi_vector_type & x ,
                     const multi_vector_type & y )
  {
    CudaSparseSingleton & s = CudaSparseSingleton::singleton();
    const float alpha = 1 , beta = 0 ;
    const int n = A.graph.row_map.extent(0) - 1 ;
    const int nz = A.graph.entries.extent(0);
    const size_t ncol = x.extent(1);

    // Sparse matrix-times-multivector
    cusparseStatus_t status =
      cusparseScsrmm( s.handle ,
                      CUSPARSE_OPERATION_NON_TRANSPOSE ,
                      n , ncol , n , nz ,
                      &alpha ,
                      s.descra ,
                      A.values.data() ,
                      A.graph.row_map.data() ,
                      A.graph.entries.data() ,
                      x.data() ,
                      n ,
                      &beta ,
                      y.data() ,
                      n );

    if ( CUSPARSE_STATUS_SUCCESS != status ) {
      throw std::runtime_error( std::string("ERROR - cusparseDcsrmv " ) );
    }
  }
};

template <>
class Multiply<
  CrsMatrix< double , Kokkos::Cuda > ,
  Kokkos::View< double** , Kokkos::LayoutLeft, Kokkos::Cuda > ,
  Kokkos::View< double** , Kokkos::LayoutLeft, Kokkos::Cuda > ,
  void ,
  IntegralRank<2> >
{
public:
  typedef Kokkos::Cuda execution_space;
  typedef execution_space::size_type size_type;
  typedef Kokkos::View< double**, Kokkos::LayoutLeft, execution_space > multi_vector_type;
  typedef CrsMatrix< double , execution_space > matrix_type;

  //--------------------------------------------------------------------------

  static void apply( const matrix_type & A ,
                     const multi_vector_type & x ,
                     const multi_vector_type & y )
  {
    CudaSparseSingleton & s = CudaSparseSingleton::singleton();
    const double alpha = 1 , beta = 0 ;
    const int n = A.graph.row_map.extent(0) - 1 ;
    const int nz = A.graph.entries.extent(0);
    const size_t ncol = x.extent(1);

    // Sparse matrix-times-multivector
    cusparseStatus_t status =
      cusparseDcsrmm( s.handle ,
                      CUSPARSE_OPERATION_NON_TRANSPOSE ,
                      n , ncol , n , nz ,
                      &alpha ,
                      s.descra ,
                      A.values.data() ,
                      A.graph.row_map.data() ,
                      A.graph.entries.data() ,
                      x.data() ,
                      n ,
                      &beta ,
                      y.data() ,
                      n );

    if ( CUSPARSE_STATUS_SUCCESS != status ) {
      throw std::runtime_error( std::string("ERROR - cusparseDcsrmv " ) );
    }
  }
};

#else
// Working on creating a version that doesn't copy vectors to a contiguous
// 2-D view and doesn't call CUSPARSE.  Not done yet.
template <typename Ordinal>
class Multiply<
  CrsMatrix< double , Kokkos::Cuda > ,
  Kokkos::View< double** , Kokkos::LayoutLeft, Kokkos::Cuda > ,
  Kokkos::View< double** , Kokkos::LayoutLeft, Kokkos::Cuda > ,
  std::vector<Ordinal> ,
  IntegralRank<2> >
{
public:
  typedef Kokkos::Cuda execution_space;
  typedef execution_space::size_type size_type;
  typedef Kokkos::View< double*, Kokkos::LayoutLeft, execution_space > vector_type;
  typedef Kokkos::View< double**, Kokkos::LayoutLeft, execution_space > multi_vector_type;
  typedef CrsMatrix< double , execution_space > matrix_type;
  typedef Kokkos::View< size_type*, execution_space > column_indices_type;

  const matrix_type m_A;
  const multi_vector_type m_x;
  multi_vector_type m_y;
  const column_indices_type m_col;
  const size_type m_num_col;

  Multiply( const matrix_type& A,
            const multi_vector_type& x,
            multi_vector_type& y,
            const column_indices_type& col) :
    m_A(A),
    m_x(x),
    m_y(y),
    m_col(col),
    m_num_col(col.extent(0)) {}

  __device__
  inline void operator() ( const size_type iRow ) const {
    const size_type iEntryBegin = m_A.graph.row_map[iRow];
    const size_type iEntryEnd   = m_A.graph.row_map[iRow+1];

    for (size_type j=0; j<m_num_col; j++) {
      size_type iCol = m_col_indices[j];

      scalar_type sum = 0.0;

      for ( size_type iEntry = iEntryBegin ; iEntry < iEntryEnd ; ++iEntry ) {
        sum += m_A.values(iEntry) * m_x(  m_A.graph.entries(iEntry), iCol );
      }

      m_y( iRow, iCol ) = sum;

    }
  }

  //--------------------------------------------------------------------------

  static void apply( const matrix_type & A ,
                     const multi_vector_type & x ,
                     const multi_vector_type & y ,
                     const std::vector<Ordinal> & col_indices )
  {
    // Copy col_indices to the device
    Kokkos::View<Ordinal*,execution_space> col_indices_dev(
      Kokkos::ViewAllocateWithoutInitializing("col_indices"), ncol);
    typename Kokkos::View<Ordinal*,execution_space>::HostMirror col_indices_host =
      Kokkos::create_mirror_view(col_indices_dev);
    for (size_t i=0; i<ncol; ++i)
      col_indices_host(i) = col_indices[i];
    Kokkos::deep_copy(col_indices_dev, col_indices_host);

    const size_t n = A.graph.row_map.extent(0) - 1 ;
    Kokkos::parallel_for( n , Multiply(A,x,y,col_indices_dev) );
  }
};

#endif

template<>
class Multiply<
  CrsMatrix< float , Kokkos::Cuda > ,
  std::vector< Kokkos::View< float* , Kokkos::Cuda > >,
  std::vector< Kokkos::View< float* , Kokkos::Cuda > >,
  void,
  IntegralRank<1> >
{
public:
  typedef Kokkos::Cuda                         execution_space ;
  typedef execution_space::size_type               size_type ;
  typedef Kokkos::View< float* , execution_space >  vector_type ;
  typedef CrsMatrix< float , execution_space >    matrix_type ;

  //--------------------------------------------------------------------------

  static void apply( const matrix_type & A ,
                     const std::vector<vector_type> & x ,
                     const std::vector<vector_type> & y )
  {
    CudaSparseSingleton & s = CudaSparseSingleton::singleton();
    const float alpha = 1 , beta = 0 ;
    const int n = A.graph.row_map.extent(0) - 1 ;
    const int nz = A.graph.entries.extent(0);
    const size_t ncol = x.size();

    // Copy columns of x into a contiguous vector
    vector_type xx( Kokkos::ViewAllocateWithoutInitializing("xx"), n * ncol );
    vector_type yy( Kokkos::ViewAllocateWithoutInitializing("yy"), n * ncol );

    for (size_t col=0; col<ncol; col++) {
      const std::pair< size_t , size_t > span( n * col , n * ( col + 1 ) );
      vector_type xx_view = Kokkos::subview( xx , span );
      Kokkos::deep_copy(xx_view, x[col]);
    }

    // Sparse matrix-times-multivector
    cusparseStatus_t status =
      cusparseScsrmm( s.handle ,
                      CUSPARSE_OPERATION_NON_TRANSPOSE ,
                      n , ncol , n , nz ,
                      &alpha ,
                      s.descra ,
                      A.values.data() ,
                      A.graph.row_map.data() ,
                      A.graph.entries.data() ,
                      xx.data() ,
                      n ,
                      &beta ,
                      yy.data() ,
                      n );

    if ( CUSPARSE_STATUS_SUCCESS != status ) {
      throw std::runtime_error( std::string("ERROR - cusparseDcsrmv " ) );
    }

    // Copy columns out of continguous multivector
    for (size_t col=0; col<ncol; col++) {
      const std::pair< size_t , size_t > span( n * col , n * ( col + 1 ) );
      vector_type yy_view = Kokkos::subview( yy , span );
      Kokkos::deep_copy(y[col], yy_view );
    }
  }
};

template<>
class Multiply<
  CrsMatrix< double , Kokkos::Cuda > ,
  std::vector< Kokkos::View< double* , Kokkos::Cuda > >,
  std::vector< Kokkos::View< double* , Kokkos::Cuda > >,
  void,
  IntegralRank<1> >
{
public:
  typedef Kokkos::Cuda                         execution_space ;
  typedef execution_space::size_type               size_type ;
  typedef Kokkos::View< double* , execution_space >  vector_type ;
  typedef CrsMatrix< double , execution_space >    matrix_type ;

  //--------------------------------------------------------------------------

  static void apply( const matrix_type & A ,
                     const std::vector<vector_type> & x ,
                     const std::vector<vector_type> & y )
  {
    CudaSparseSingleton & s = CudaSparseSingleton::singleton();
    const double alpha = 1 , beta = 0 ;
    const int n = A.graph.row_map.extent(0) - 1 ;
    const int nz = A.graph.entries.extent(0);
    const size_t ncol = x.size();

    // Copy columns of x into a contiguous vector
    vector_type xx( Kokkos::ViewAllocateWithoutInitializing("xx"), n * ncol );
    vector_type yy( Kokkos::ViewAllocateWithoutInitializing("yy"), n * ncol );

    for (size_t col=0; col<ncol; col++) {
      const std::pair< size_t , size_t > span( n * col , n * ( col + 1 ) );
      vector_type xx_view = Kokkos::subview( xx , span );
      Kokkos::deep_copy(xx_view, x[col]);
    }

    // Sparse matrix-times-multivector
    cusparseStatus_t status =
      cusparseDcsrmm( s.handle ,
                      CUSPARSE_OPERATION_NON_TRANSPOSE ,
                      n , ncol , n , nz ,
                      &alpha ,
                      s.descra ,
                      A.values.data() ,
                      A.graph.row_map.data() ,
                      A.graph.entries.data() ,
                      xx.data() ,
                      n ,
                      &beta ,
                      yy.data() ,
                      n );

    if ( CUSPARSE_STATUS_SUCCESS != status ) {
      throw std::runtime_error( std::string("ERROR - cusparseDcsrmv " ) );
    }

    // Copy columns out of continguous multivector
    for (size_t col=0; col<ncol; col++) {
      const std::pair< size_t , size_t > span( n * col , n * ( col + 1 ) );
      vector_type yy_view = Kokkos::subview( yy , span );
      Kokkos::deep_copy(y[col], yy_view );
    }
  }
};

//----------------------------------------------------------------------------

} // namespace Stokhos

#endif /* #ifdef HAVE_STOKHOS_CUSPARSE */

#endif /* #ifndef STOKHOS_CUDA_CRSMATRIX_HPP */
