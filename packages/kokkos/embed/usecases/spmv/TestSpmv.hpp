
#include <iostream>
#include <typeinfo>
#include <utility>

#include <KokkosArray_View.hpp>

#include <KokkosArray_Array.hpp>
#include <impl/KokkosArray_ArrayAnalyzeShape.hpp>
#include <impl/KokkosArray_ArrayViewDefault.hpp>

#include <impl/KokkosArray_Timer.hpp>


#include <TestCrsMatrix.hpp>
#include <TestGenerateSystem.hpp>


namespace Test {

struct PerfSpmv {
  double seconds_per_iter ;
  double giga_flops ;
  double giga_bytes ;
  size_t row_count ;
  size_t entry_count ;
};

template< typename ScalarType , class Device >
PerfSpmv test_spmv_scalar( const int nGrid ,
                           const int iterCount ,
                           const char * const verify_label )
{
  typedef ScalarType value_type ;

  typedef KokkosArray::View< value_type* , KokkosArray::LayoutRight , Device > vector_type ;

  //------------------------------

  typedef KokkosArray::CrsMatrix<value_type,Device> matrix_type ;
  typedef KokkosArray::CrsArray<int,Device,Device,int> crsarray_type ;

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

  //------------------------------

  PerfSpmv result = { 0 , 0 , 0 , 0 , 0 };

  for ( int iter = 0 ; iter < iterCount ; ++iter ) {

    KokkosArray::Impl::Timer clock ;
    for ( int i = 0 ; i < 100 ; ++i ) {
      KokkosArray::multiply( matrix , x , y );
    }
    Device::fence();

    const double dt = clock.seconds() / ((double) 100 );

    result.seconds_per_iter =
      0 < result.seconds_per_iter && result.seconds_per_iter < dt ? result.seconds_per_iter : dt ;
  }

  result.row_count   = fem_length ;
  result.entry_count = fem_graph_length ;
  result.giga_flops  = 1.0e-9 * 2.0 * fem_graph_length ;
  result.giga_bytes  = 1.0e-9 *
    ( hy_result.dimension_0() * 2 * sizeof(int)             // read: row entry range
    + matrix.values.dimension_0() * sizeof(int)             // read: column index
    + matrix.values.dimension_0() * 2 * sizeof(ScalarType)  // read: coefficient A and X value
    + hy_result.dimension_0() * sizeof(ScalarType)          // write: Y value
    );

  //------------------------------

  if ( 0 != verify_label ) {

    std::cout << "test_spmv verify " << verify_label << " : "
              << fem_length << " x 1" << std::endl ;

    KokkosArray::deep_copy( hx , y );

    for ( unsigned i = 0 ; i < fem_length ; ++i ) {
      const double tol  = KokkosArray::Impl::is_same<ScalarType,double>::value ? 1.0e-14 : 1.0e-5 ;
      const double diff = fabs( hx(i) - hy_result(i) );
      const double mag  = 0 < fabs( hy_result(i) ) ? fabs( hy_result(i) ) : 1 ;
      if ( diff > tol && ( diff / mag ) > tol ) {
        std::cout << "error at(" << i << ") = " << hx(i) << " != " << hy_result(i) << std::endl ;
      }
    }
  }

  //------------------------------

  return result ;
}

template< typename ScalarType , unsigned N , class Device >
PerfSpmv test_spmv( const int nGrid ,
                    const int iterCount ,
                    const char * const verify_label )
{
  typedef KokkosArray::Array<ScalarType,N> value_type ;

  typedef KokkosArray::View< value_type* , KokkosArray::LayoutRight , Device > vector_type ;

  //------------------------------

  typedef KokkosArray::CrsMatrix<value_type,Device> matrix_type ;
  typedef KokkosArray::CrsArray<int,Device,Device,int> crsarray_type ;

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

  //------------------------------

  PerfSpmv result = { 0 , 0 , 0 , 0 , 0 };

  for ( int iter = 0 ; iter < iterCount ; ++iter ) {

    KokkosArray::Impl::Timer clock ;
    for ( int i = 0 ; i < 100 ; ++i ) {
      KokkosArray::multiply( matrix , x , y );
    }
    Device::fence();

    const double dt = clock.seconds() / ((double) 100 );

    result.seconds_per_iter =
      0 < result.seconds_per_iter && result.seconds_per_iter < dt ? result.seconds_per_iter : dt ;
  }

  result.row_count   = fem_length ;
  result.entry_count = fem_graph_length ;
  result.giga_flops  = 1.0e-9 * 2.0 * fem_graph_length * N ;
  result.giga_bytes  = 1.0e-9 *
    ( hy_result.dimension_0() * 2 * sizeof(int)                 // read: row entry range
    + matrix.values.dimension_0() * sizeof(int)                 // read: column index
    + matrix.values.dimension_0() * N * 2 * sizeof(ScalarType)  // read: coefficient A and X value
    + hy_result.dimension_0() * N * sizeof(ScalarType)          // write: Y value
    );

  //------------------------------

  if ( 0 != verify_label ) {

    std::cout << "test_spmv verify " << verify_label << " : "
              << fem_length << " x " << N << std::endl ;

    KokkosArray::deep_copy( hx , y );

    for ( unsigned i = 0 ; i < fem_length ; ++i ) {
      for ( unsigned j = 0 ; j < N ; ++j ) {
        const double tol  = KokkosArray::Impl::is_same<ScalarType,double>::value ? 1.0e-14 : 1.0e-5 ;
        const double diff = fabs( hx(i)[j] - hy_result(i)[j] );
        const double mag  = 0 < fabs( hy_result(i)[j] ) ? fabs( hy_result(i)[j] ) : 1 ;
        if ( diff > tol && ( diff / mag ) > tol ) {
          std::cout << "error at(" << i << ")[" << j << "] = " << hx(i)[j] << " != " << hy_result(i)[j] << std::endl ;
        }
      }
    }
  }

  //------------------------------

  return result ;
}

//----------------------------------------------------------------------------

template< class Device >
void test_spmv_driver( const char * const label )
{
  PerfSpmv perf_spmv , perf_spmv_scalar ;

  perf_spmv_scalar = test_spmv_scalar<double,Device>( 5 , 1 , label );
  perf_spmv        = test_spmv<double,8,Device>( 5 , 1 , label );

  std::cout << std::endl ;
  std::cout << "\"SPMV " << label << " Sample: scalar\"" << std::endl ;
  std::cout << "\"ROWS\" , \"ENTRIES\" , \"TIME\" , \"GBYTE\" , \"GBYTE/SEC\" , \"GFLOP\" , \"GFLOP/SEC\"" << std::endl ;


  for ( unsigned j = 8 ; j <= 128 ; j *= 2 ) {
    perf_spmv_scalar = test_spmv_scalar<double,Device>( j , 10 , 0 );

    std::cout << perf_spmv_scalar.row_count << " , "
              << perf_spmv_scalar.entry_count << " , "
              << perf_spmv_scalar.seconds_per_iter << " , "
              << perf_spmv_scalar.giga_bytes << " , "
              << perf_spmv_scalar.giga_bytes / perf_spmv_scalar.seconds_per_iter << " , "
              << perf_spmv_scalar.giga_flops << " , "
              << perf_spmv_scalar.giga_flops / perf_spmv_scalar.seconds_per_iter
              << std::endl ;
  }

  std::cout << std::endl ;
  std::cout << "\"SPMV " << label << " Samples: Array<32>\"" << std::endl ;
  std::cout << "\"ROWS\" , \"ENTRIES\" , \"TIME\" , \"GBYTE\" , \"GBYTE/SEC\" , \"GFLOP\" , \"GFLOP/SEC\"" << std::endl ;

  for ( unsigned j = 4 ; j <= 64 ; j *= 2 ) {
    perf_spmv = test_spmv<double,32,Device>( j , 10 , 0 );

    std::cout << perf_spmv.row_count * 32 << " , "
              << perf_spmv.entry_count * 32 << " , "
              << perf_spmv.seconds_per_iter << " , "
              << perf_spmv.giga_bytes << " , "
              << perf_spmv.giga_bytes / perf_spmv.seconds_per_iter << " , "
              << perf_spmv.giga_flops << " , "
              << perf_spmv.giga_flops / perf_spmv.seconds_per_iter
              << std::endl ;
  }

  std::cout << std::endl ;
}

}

