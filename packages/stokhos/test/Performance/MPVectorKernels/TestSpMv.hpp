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
#include <iostream>

// MP::Vector and Matrix
#include "Stokhos_Sacado_Kokkos.hpp"
#include "Kokkos_CrsMatrix.hpp"
#include "Kokkos_CrsMatrix_MP_Vector.hpp"

// Compile-time loops
#include "Sacado_mpl_range_c.hpp"
#include "Sacado_mpl_for_each.hpp"
#include "Sacado_mpl_integral_c.hpp"

// Utilities
#include "impl/Kokkos_Timer.hpp"

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

template <typename StorageType, typename MultiplyTag>
std::vector<double>
test_mpvector_spmv(const int ensemble_length,
                   const int nGrid,
                   const int iterCount,
                   Kokkos::DeviceConfig dev_config,
                   MultiplyTag tag)
{
  typedef StorageType storage_type;
  typedef typename storage_type::value_type value_type;
  typedef typename storage_type::ordinal_type ordinal_type;
  typedef typename storage_type::execution_space execution_space;
  typedef Sacado::MP::Vector<StorageType> VectorType;
  typedef Kokkos::LayoutRight Layout;
  typedef Kokkos::View< VectorType*, Layout, execution_space > vector_type;
  typedef Kokkos::CrsMatrix< VectorType, ordinal_type, execution_space > matrix_type;
  typedef typename matrix_type::StaticCrsGraphType matrix_graph_type;
  typedef typename matrix_type::values_type matrix_values_type;

  //------------------------------
  // Generate graph for "FEM" box structure:

  std::vector< std::vector<size_t> > fem_graph;
  const size_t fem_length = nGrid * nGrid * nGrid;
  const size_t graph_length = generate_fem_graph( nGrid , fem_graph );

  //------------------------------
  // Generate input multivector:

  vector_type x =
    vector_type(Kokkos::ViewAllocateWithoutInitializing("x"), fem_length, ensemble_length);
  vector_type y =
    vector_type(Kokkos::ViewAllocateWithoutInitializing("y"), fem_length, ensemble_length);

  //------------------------------

  matrix_graph_type matrix_graph =
    Kokkos::create_staticcrsgraph<matrix_graph_type>(
      std::string("test crs graph"), fem_graph);
  matrix_values_type matrix_values =
    matrix_values_type(Kokkos::ViewAllocateWithoutInitializing("matrix"), graph_length, ensemble_length);
  matrix_type matrix("block_matrix", fem_length, matrix_values, matrix_graph);
  matrix.dev_config = dev_config;

  //------------------------------
  // Fill:

  {
    // The VectorType may be dynamic (with allocated memory)
    // so cannot pass a VectorType value to the device.
    // Get an array-of-intrinsic View and fill that view.
    typename vector_type::array_type xx( x );
    typename vector_type::array_type yy( y );
    typename matrix_values_type::array_type mm( matrix_values );

    Kokkos::deep_copy( xx , value_type(1.0) );
    Kokkos::deep_copy( yy , value_type(1.0) );
    Kokkos::deep_copy( mm , value_type(1.0) );
  }

  //------------------------------

  // One iteration to warm up
  Stokhos::multiply( matrix, x, y, tag );

  execution_space::fence();
  Kokkos::Impl::Timer clock ;
  for (int iter = 0; iter < iterCount; ++iter) {
    Stokhos::multiply( matrix, x, y, tag );
  }
  execution_space::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );
  const double flops = 1.0e-9 * 2.0 * graph_length * ensemble_length;

  std::vector<double> perf(5);
  perf[0] = fem_length;
  perf[1] = ensemble_length;
  perf[2] = graph_length;
  perf[3] = seconds_per_iter;
  perf[4] = flops / seconds_per_iter;
  return perf;
}

template <typename ScalarType, typename OrdinalType, typename Device>
std::vector<double>
test_scalar_spmv(const int ensemble_length,
                 const int nGrid,
                 const int iterCount,
                 Kokkos::DeviceConfig dev_config)
{
  typedef ScalarType value_type;
  typedef OrdinalType ordinal_type;
  typedef Device execution_space;
  typedef Kokkos::View< value_type*, execution_space > vector_type;
  typedef Kokkos::CrsMatrix< value_type, ordinal_type, execution_space > matrix_type;
  typedef typename matrix_type::StaticCrsGraphType matrix_graph_type;
  typedef typename matrix_type::values_type matrix_values_type;

  //------------------------------
  // Generate graph for "FEM" box structure:

  std::vector< std::vector<size_t> > fem_graph;
  const size_t fem_length = nGrid * nGrid * nGrid;
  const size_t graph_length = generate_fem_graph( nGrid , fem_graph );

  //------------------------------
  // Generate input multivector:

  std::vector<vector_type> x(ensemble_length);
  std::vector<vector_type> y(ensemble_length);
  for (int e=0; e<ensemble_length; ++e) {
    x[e] = vector_type(Kokkos::ViewAllocateWithoutInitializing("x"), fem_length);
    y[e] = vector_type(Kokkos::ViewAllocateWithoutInitializing("y"), fem_length);

    Kokkos::deep_copy( x[e] , value_type(1.0) );
    Kokkos::deep_copy( y[e] , value_type(0.0) );
  }

  //------------------------------

  std::vector<matrix_type> matrix(ensemble_length);
  for (int e=0; e<ensemble_length; ++e) {
    matrix_graph_type matrix_graph =
      Kokkos::create_staticcrsgraph<matrix_graph_type>(
        std::string("test crs graph"), fem_graph);
    matrix_values_type matrix_values =
      matrix_values_type(Kokkos::ViewAllocateWithoutInitializing("matrix"), graph_length);
    matrix[e] = matrix_type("matrix", fem_length, matrix_values, matrix_graph);

    Kokkos::deep_copy( matrix[e].values , value_type(1.0) );
  }

  //------------------------------

  // One iteration to warm up
  for (int iter = 0; iter < iterCount; ++iter) {
    for (int e=0; e<ensemble_length; ++e) {
      Kokkos::MV_Multiply( y[e], matrix[e], x[e] );
    }
  }

  execution_space::fence();
  Kokkos::Impl::Timer clock ;
  for (int iter = 0; iter < iterCount; ++iter) {
    for (int e=0; e<ensemble_length; ++e) {
      Kokkos::MV_Multiply( y[e], matrix[e], x[e] );
    }
  }
  execution_space::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );
  const double flops = 1.0e-9 * 2.0 * graph_length * ensemble_length;

  std::vector<double> perf(5);
  perf[0] = fem_length;
  perf[1] = ensemble_length;
  perf[2] = graph_length;
  perf[3] = seconds_per_iter;
  perf[4] = flops / seconds_per_iter;
  return perf;
}

template <class Storage>
struct PerformanceDriverOp {
  typedef typename Storage::value_type Scalar;
  typedef typename Storage::ordinal_type Ordinal;
  typedef typename Storage::execution_space Device;
  const int nGrid, nIter;
  Kokkos::DeviceConfig dev_config;

  PerformanceDriverOp(const int nGrid_, const int nIter_,
                      Kokkos::DeviceConfig dev_config_) :
    nGrid(nGrid_), nIter(nIter_), dev_config(dev_config_) {}

  template <typename ArgT>
  void operator() (ArgT arg) const {
    const int ensemble = ArgT::value;
    typedef typename Storage::template apply_N<ensemble> NewStorageApply;
    typedef typename NewStorageApply::type storage_type;

    const std::vector<double> perf_scalar =
      test_scalar_spmv<Scalar,Ordinal,Device>(
        ensemble, nGrid, nIter, dev_config );

    const std::vector<double> perf_mpvector =
      test_mpvector_spmv<storage_type>(
        ensemble, nGrid, nIter, dev_config, Stokhos::DefaultMultiply() );

    std::cout << nGrid << " , "
              << perf_scalar[0] << " , "
              << perf_scalar[2] << " , "
              << perf_scalar[1] << " , "
              << perf_scalar[3] << " , "
              << perf_scalar[4] / perf_scalar[4] << " , "
              << perf_scalar[4] << " , "
              << perf_mpvector[4]/ perf_scalar[4] << " , "
              << perf_mpvector[4] << " , "
              << std::endl;
  }
};

template <class Storage, int entry_min, int entry_max, int entry_step>
void performance_test_driver( const int nGrid,
                              const int nIter,
                              Kokkos::DeviceConfig dev_config)
{
  std::cout.precision(8);
  std::cout << std::endl
            << "\"Grid Size\" , "
            << "\"FEM Size\" , "
            << "\"FEM Graph Size\" , "
            << "\"Ensemble Size\" , "
            << "\"Scalar SpMv Time\" , "
            << "\"Scalar SpMv Speedup\" , "
            << "\"Scalar SpMv GFLOPS\" , "
            << "\"MPVector SpMv Speedup\" , "
            << "\"MPVector SpMv GFLOPS\" , "
            << std::endl;

  // Loop over [entry_min, entry_max] vector entries per thread
  typedef Sacado::mpl::range_c< int, entry_min, entry_max+1, entry_step > Range;
  PerformanceDriverOp<Storage> op(nGrid, nIter, dev_config);
  Sacado::mpl::for_each<Range> f(op);
}
