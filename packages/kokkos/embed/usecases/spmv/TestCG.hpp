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

#ifndef KOKKOSARRAY_CG_HPP
#define KOKKOSARRAY_CG_HPP

#include <cmath>

#include <KokkosArray_View.hpp>
#include <KokkosArray_Array.hpp>
#include <impl/KokkosArray_ArrayAnalyzeShape.hpp>
#include <impl/KokkosArray_ArrayViewDefault.hpp>

#include <impl/KokkosArray_Timer.hpp>

#include <TestBlas1.hpp>
#include <TestCrsMatrix.hpp>
#include <TestGenerateSystem.hpp>

//----------------------------------------------------------------------------

namespace KokkosArray {

template< typename AScalarType , typename VScalarType , class Device >
void cgsolve(
  const CrsMatrix<AScalarType,Device>  & A ,
  const View<VScalarType*,LayoutRight,Device> & b ,
  const View<VScalarType*,LayoutRight,Device> & x ,
  size_t & iteration ,
  double & normr ,
  double & iter_time ,
  const size_t maximum_iteration = 200 ,
  const double tolerance = std::numeric_limits<VScalarType>::epsilon() )
{
  typedef View<VScalarType*,LayoutRight,Device> vector_type ;

  const size_t count = b.dimension_0();

  vector_type p ( "cg::p" , count );
  vector_type r ( "cg::r" , count );
  vector_type Ap( "cg::Ap", count );

  /* r = b - A * x ; */

  /* p  = x      */ deep_copy( p , x );
  /* Ap = A * p  */ multiply( A , p , Ap );
  /* r  = b - Ap */ waxpby( count , 1.0 , b , -1.0 , Ap , r );
  /* p  = r      */ deep_copy( p , r );

  double old_rdot = dot( count , r );

  normr     = std::sqrt( old_rdot );
  iteration = 0 ;

  KokkosArray::Impl::Timer wall_clock ;

  while ( tolerance < normr && iteration < maximum_iteration ) {

    /* pAp_dot = dot( p , Ap = A * p ) */

    /* Ap = A * p  */ multiply( A , p , Ap );

    const double pAp_dot = dot( count , p , Ap );
    const double alpha   = old_rdot / pAp_dot ;

    /* x += alpha * p ;  */ axpy( count,  alpha, p , x );
    /* r -= alpha * Ap ; */ axpy( count, -alpha, Ap, r );

    const double r_dot = dot( count , r );
    const double beta  = r_dot / old_rdot ;

    /* p = r + beta * p ; */ xpby( count , r , beta , p );

    normr = std::sqrt( old_rdot = r_dot );
    ++iteration ;
  }

  iter_time = wall_clock.seconds();
}

}

//----------------------------------------------------------------------------

namespace Test {

struct PerfCGSolve {
  double seconds_per_iter ;
  size_t row_count ;
  size_t entry_count ;
};

template< typename Scalar , class Device >
PerfCGSolve test_cgsolve_scalar( const int nGrid ,
                                const int iterMax ,
                                const char * const /* verify_label */ )
{
  typedef Scalar value_type ;

  typedef KokkosArray::CrsArray<int,Device,Device,int>      crsarray_type ;
  typedef KokkosArray::CrsMatrix<value_type,Device>         matrix_type ;
  typedef KokkosArray::View<value_type*,KokkosArray::LayoutRight,Device> vector_type ;

  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const unsigned fem_length = nGrid * nGrid * nGrid ;

  Test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------

  vector_type x = vector_type( "x" , fem_length );
  vector_type y = vector_type( "y" , fem_length );

  typename vector_type::HostMirror hx        = KokkosArray::create_mirror( x );
  typename vector_type::HostMirror hy_result = KokkosArray::create_mirror( y );

  for ( unsigned i = 0 ; i < fem_length ; ++i ) {
    hx(i) = Test::generate_vector_coefficient( fem_length , 1 , i , 0 );
  }

  KokkosArray::deep_copy( x , hx );

  //------------------------------

  matrix_type matrix ;

  matrix.graph = KokkosArray::create_crsarray<crsarray_type>( std::string("testing") , fem_graph );

  const unsigned fem_graph_length = matrix.graph.entries.dimension_0();

  matrix.values = vector_type( "matrix" , fem_graph_length );

  {
    typename vector_type::HostMirror hM =
      KokkosArray::create_mirror( matrix.values );

    for ( size_t iRow = 0 , iEntry = 0 ; iRow < fem_length ; ++iRow ) {

      hy_result(iRow) = 0 ;

      for ( size_t iRowEntry = 0 ; iRowEntry < fem_graph[ iRow ].size() ; ++iRowEntry , ++iEntry ) {

        const size_t iCol = fem_graph[ iRow ][ iRowEntry ];

        hM(iEntry) = Test::generate_matrix_coefficient( fem_length , 1 , iRow, iCol, 0 );
        hy_result(iRow) += hM(iEntry) * hx(iCol);
      }
    }

    KokkosArray::deep_copy( matrix.values , hM );
  }

  size_t iter_count = 0 ;
  double iter_time  = 0 ;
  double norm_resid = 0 ;

  cgsolve( matrix , x , y , iter_count , norm_resid , iter_time , iterMax , 1e-14 );

  PerfCGSolve perf ;

  perf.seconds_per_iter = iter_time ;
  perf.row_count   = fem_length ;
  perf.entry_count = fem_graph_length ;

  return perf ;
}

//----------------------------------------------------------------------------

template< typename Scalar , unsigned N , class Device >
PerfCGSolve test_cgsolve_array( const int nGrid ,
                                const int iterMax ,
                                const char * const /* verify_label */ )
{
  typedef KokkosArray::Array<Scalar,N> value_type ;

  typedef KokkosArray::CrsArray<int,Device,Device,int>      crsarray_type ;
  typedef KokkosArray::CrsMatrix<value_type,Device>         matrix_type ;
  typedef KokkosArray::View<value_type*,KokkosArray::LayoutRight,Device> vector_type ;

  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const unsigned fem_length = nGrid * nGrid * nGrid ;

  Test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------

  vector_type x = vector_type( "x" , fem_length );
  vector_type y = vector_type( "y" , fem_length );

  typename vector_type::HostMirror hx        = KokkosArray::create_mirror( x );
  typename vector_type::HostMirror hy_result = KokkosArray::create_mirror( y );

  for ( unsigned i = 0 ; i < fem_length ; ++i ) {
    for ( unsigned j = 0 ; j < N ; ++j ) {
      hx(i)[j] = Test::generate_vector_coefficient( fem_length , N , i , j );
    }
  }

  KokkosArray::deep_copy( x , hx );

  //------------------------------

  matrix_type matrix ;

  matrix.graph = KokkosArray::create_crsarray<crsarray_type>( std::string("testing") , fem_graph );

  const unsigned fem_graph_length = matrix.graph.entries.dimension_0();

  matrix.values = vector_type( "matrix" , fem_graph_length );

  {
    typename vector_type::HostMirror hM =
      KokkosArray::create_mirror( matrix.values );

    for ( size_t iRow = 0 , iEntry = 0 ; iRow < fem_length ; ++iRow ) {

      for ( unsigned k = 0 ; k < N ; ++k ) { hy_result(iRow)[k] = 0 ; }

      for ( size_t iRowEntry = 0 ; iRowEntry < fem_graph[ iRow ].size() ; ++iRowEntry , ++iEntry ) {

        const size_t iCol = fem_graph[ iRow ][ iRowEntry ];

        for ( unsigned k = 0 ; k < N ; ++k ) {
          hM(iEntry)[k] = Test::generate_matrix_coefficient( fem_length , N , iRow, iCol, k );
          hy_result(iRow)[k] += hM(iEntry)[k] * hx(iCol)[k];
        }
      }
    }

    KokkosArray::deep_copy( matrix.values , hM );
  }

  size_t iter_count = 0 ;
  double iter_time  = 0 ;
  double norm_resid = 0 ;

  cgsolve( matrix , x , y , iter_count , norm_resid , iter_time , iterMax , 1e-14 );

  PerfCGSolve perf ;

  perf.seconds_per_iter = iter_time ;
  perf.row_count   = fem_length ;
  perf.entry_count = fem_graph_length ;

  return perf ;
}


//----------------------------------------------------------------------------

template< class Device >
void test_cgsolve_driver( const char * const label )
{
  PerfCGSolve perf_array ;

  std::cout << std::endl ;
  std::cout << "\"CGSolve " << label << "\" Samples: scalar\"" << std::endl ;
  std::cout << "\"FEM-ROWS\" , \"FEM-ENTRIES\" , \"TIME/ITER\" , \"TIME/ITER/ROW\"" << std::endl ;

  for ( unsigned j = 4 ; j <= 128 ; j *= 2 ) {
    perf_array = test_cgsolve_scalar<double,Device>( j , 100 , 0 );

    std::cout << perf_array.row_count << " , "
              << perf_array.entry_count << " , "
              << perf_array.seconds_per_iter << " , "
              << perf_array.seconds_per_iter / perf_array.row_count
              << std::endl ;
  }

  std::cout << std::endl ;
  std::cout << "\"CGSolve " << label << "\" Samples: Array<32>\"" << std::endl ;
  std::cout << "\"FEM-ROWS\" , \"FEM-ENTRIES\" , \"TOTAL-ROWS\" , \"TOTAL-ENTRIES\" , \"TIME/ITER\" , \"TIME/ITER/ROW\"" << std::endl ;

  for ( unsigned j = 4 ; j <= 128 ; j *= 2 ) {
    perf_array = test_cgsolve_array<double,32,Device>( j , 100 , 0 );

    std::cout << perf_array.row_count << " , "
              << perf_array.entry_count << " , "
              << perf_array.row_count * 32 << " , "
              << perf_array.entry_count * 32 << " , "
              << perf_array.seconds_per_iter << " , "
              << perf_array.seconds_per_iter / perf_array.row_count
              << std::endl ;
  }

  std::cout << std::endl ;

}

}


#endif

