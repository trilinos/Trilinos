// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

// Kokkos CrsMatrix
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_spmv.hpp"


// Utilities
#include "Kokkos_Timer.hpp"

template< typename IntType >
inline
IntType map_fem_graph_coord( const IntType & N ,
                             const IntType & i ,
                             const IntType & j ,
                             const IntType & k )
{
  return k + N * ( j + N * i );
}

inline
size_t generate_fem_graph( size_t N ,
                           std::vector< std::vector<size_t> > & graph )
{
  graph.resize( N * N * N , std::vector<size_t>() );

  size_t total = 0 ;

  for ( int i = 0 ; i < (int) N ; ++i ) {
  for ( int j = 0 ; j < (int) N ; ++j ) {
  for ( int k = 0 ; k < (int) N ; ++k ) {

    const size_t row = map_fem_graph_coord((int)N,i,j,k);

    graph[row].reserve(27);

    for ( int ii = -1 ; ii < 2 ; ++ii ) {
    for ( int jj = -1 ; jj < 2 ; ++jj ) {
    for ( int kk = -1 ; kk < 2 ; ++kk ) {
      if ( 0 <= i + ii && i + ii < (int) N &&
           0 <= j + jj && j + jj < (int) N &&
           0 <= k + kk && k + kk < (int) N ) {
        size_t col = map_fem_graph_coord((int)N,i+ii,j+jj,k+kk);

        graph[row].push_back(col);
      }
    }}}
    total += graph[row].size();
  }}}

  return total ;
}

template <typename ScalarType, typename OrdinalType, typename Device>
void
test_spmm(const OrdinalType ensemble_length,
          const OrdinalType nGrid,
          const OrdinalType iterCount,
          std::vector<double>& scalar_perf,
          std::vector<double>& block_left_perf,
          std::vector<double>& block_right_perf)
{
  typedef ScalarType value_type;
  typedef OrdinalType ordinal_type;
  typedef Device execution_space;
  typedef Kokkos::View< value_type*, execution_space > vector_type;
  typedef Kokkos::View< value_type**, Kokkos::LayoutLeft, execution_space > left_multivec_type;
  //typedef Kokkos::View< value_type**, Kokkos::LayoutRight, execution_space > right_multivec_type;
  typedef KokkosSparse::CrsMatrix< value_type, ordinal_type, execution_space > matrix_type;
  typedef typename matrix_type::StaticCrsGraphType matrix_graph_type;
  typedef typename matrix_type::values_type matrix_values_type;

  //------------------------------
  // Generate graph for "FEM" box structure:

  std::vector< std::vector<size_t> > fem_graph;
  const size_t fem_length = nGrid * nGrid * nGrid;
  const size_t graph_length = generate_fem_graph( nGrid , fem_graph );

  //------------------------------
  // Generate input vectors:

  std::vector<vector_type> x(ensemble_length);
  std::vector<vector_type> y(ensemble_length);
  for (ordinal_type e=0; e<ensemble_length; ++e) {
    x[e] = vector_type(Kokkos::ViewAllocateWithoutInitializing("x"), fem_length);
    y[e] = vector_type(Kokkos::ViewAllocateWithoutInitializing("y"), fem_length);

    Kokkos::deep_copy( x[e] , value_type(1.0) );
    Kokkos::deep_copy( y[e] , value_type(0.0) );
  }
  left_multivec_type xl(Kokkos::ViewAllocateWithoutInitializing("xl"), fem_length, ensemble_length);
  left_multivec_type yl(Kokkos::ViewAllocateWithoutInitializing("yl"), fem_length, ensemble_length);
  // right_multivec_type xr(Kokkos::ViewAllocateWithoutInitializing("xr"), fem_length, ensemble_length);
  // right_multivec_type yr(Kokkos::ViewAllocateWithoutInitializing("yr"), fem_length, ensemble_length);
  Kokkos::deep_copy(xl, value_type(1.0));
  //Kokkos::deep_copy(xr, value_type(1.0));
  Kokkos::deep_copy(yl, value_type(0.0));
  //Kokkos::deep_copy(yr, value_type(0.0));

  //------------------------------
  // Generate matrix

  matrix_graph_type matrix_graph =
    Kokkos::create_staticcrsgraph<matrix_graph_type>(
      std::string("test crs graph"), fem_graph);
  matrix_values_type matrix_values =
    matrix_values_type(Kokkos::ViewAllocateWithoutInitializing("matrix"), graph_length);
  matrix_type matrix("matrix", fem_length, matrix_values, matrix_graph);
  Kokkos::deep_copy( matrix_values , value_type(1.0) );

  //------------------------------
  // Scalar multiply

  {
    // warm up
    for (ordinal_type iter = 0; iter < iterCount; ++iter) {
      for (ordinal_type e=0; e<ensemble_length; ++e) {
        KokkosSparse::spmv( "N", value_type(1.0), matrix, x[e] , value_type(0.0) , y[e]);
      }
    }

    execution_space().fence();
    Kokkos::Timer clock ;
    for (ordinal_type iter = 0; iter < iterCount; ++iter) {
      for (ordinal_type e=0; e<ensemble_length; ++e) {
        KokkosSparse::spmv( "N", value_type(1.0), matrix, x[e] , value_type(0.0) , y[e]);
      }
    }
    execution_space().fence();

    const double seconds_per_iter = clock.seconds() / ((double) iterCount );
    const double flops = 1.0e-9 * 2.0 * graph_length * ensemble_length;

    scalar_perf.resize(5);
    scalar_perf[0] = fem_length;
    scalar_perf[1] = ensemble_length;
    scalar_perf[2] = graph_length;
    scalar_perf[3] = seconds_per_iter;
    scalar_perf[4] = flops / seconds_per_iter;
  }

  //------------------------------
  // Block-left multiply

  {
    // warm up
    for (ordinal_type iter = 0; iter < iterCount; ++iter) {
      KokkosSparse::spmv( "N", value_type(1.0), matrix, xl , value_type(0.0) , yl);
    }

    execution_space().fence();
    Kokkos::Timer clock ;
    for (ordinal_type iter = 0; iter < iterCount; ++iter) {
      KokkosSparse::spmv( "N", value_type(1.0), matrix, xl , value_type(0.0) , yl);
    }
    execution_space().fence();

    const double seconds_per_iter = clock.seconds() / ((double) iterCount );
    const double flops = 1.0e-9 * 2.0 * graph_length * ensemble_length;

    block_left_perf.resize(5);
    block_left_perf[0] = fem_length;
    block_left_perf[1] = ensemble_length;
    block_left_perf[2] = graph_length;
    block_left_perf[3] = seconds_per_iter;
    block_left_perf[4] = flops / seconds_per_iter;
  }

#if 0
  //------------------------------
  // Block-right multiply

  {
    // warm up
    for (ordinal_type iter = 0; iter < iterCount; ++iter) {
      KokkosSparse::spmv( "N", value_type(1.0), matrix, xr , value_type(0.0) , yr);
    }

    execution_space().fence();
    Kokkos::Timer clock ;
    for (ordinal_type iter = 0; iter < iterCount; ++iter) {
      KokkosSparse::spmv( "N", value_type(1.0), matrix, xr , value_type(0.0) , yr);
    }
    execution_space().fence();

    const double seconds_per_iter = clock.seconds() / ((double) iterCount );
    const double flops = 1.0e-9 * 2.0 * graph_length * ensemble_length;

    block_right_perf.resize(5);
    block_right_perf[0] = fem_length;
    block_right_perf[1] = ensemble_length;
    block_right_perf[2] = graph_length;
    block_right_perf[3] = seconds_per_iter;
    block_right_perf[4] = flops / seconds_per_iter;
  }
#endif

}

template <typename Scalar, typename Ordinal, typename Device>
void performance_test_driver( const Ordinal nGrid,
                              const Ordinal nIter,
                              const Ordinal ensemble_min,
                              const Ordinal ensemble_max,
                              const Ordinal ensemble_step )
{
  std::cout.precision(8);
  std::cout << std::endl
            << "\"Grid Size\" , "
            << "\"FEM Size\" , "
            << "\"FEM Graph Size\" , "
            << "\"Ensemble Size\" , "
            << "\"Scalar SpMM Time\" , "
            << "\"Scalar SpMM Speedup\" , "
            << "\"Scalar SpMM GFLOPS\" , "
            << "\"Block-Left SpMM Speedup\" , "
            << "\"Block-Left SpMM GFLOPS\" , "
    //<< "\"Block_Right SpMM Speedup\" , "
    //<< "\"Block_Right SpMM GFLOPS\" , "
            << std::endl;

  std::vector<double> perf_scalar, perf_block_left, perf_block_right;
  for (Ordinal e=ensemble_min; e<=ensemble_max; e+=ensemble_step) {

    test_spmm<Scalar,Ordinal,Device>(
      e, nGrid, nIter, perf_scalar, perf_block_left, perf_block_right );

    std::cout << nGrid << " , "
              << perf_scalar[0] << " , "
              << perf_scalar[2] << " , "
              << perf_scalar[1] << " , "
              << perf_scalar[3] << " , "
              << perf_scalar[4] / perf_scalar[4] << " , "
              << perf_scalar[4] << " , "
              << perf_block_left[4]/ perf_scalar[4] << " , "
              << perf_block_left[4] << " , "
      //<< perf_block_right[4]/ perf_scalar[4] << " , "
      //<< perf_block_right[4] << " , "
              << std::endl;

  }
}
