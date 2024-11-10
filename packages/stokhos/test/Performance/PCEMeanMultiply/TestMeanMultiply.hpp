// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

// PCE scalar type
#include "Stokhos_Sacado_Kokkos_UQ_PCE.hpp"

// Kokkos CrsMatrix
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_spmv.hpp"
#include "Kokkos_CrsMatrix_UQ_PCE.hpp"
#include "Kokkos_CrsMatrix_UQ_PCE_Cuda.hpp"

// Stokhos
#include "Stokhos_LegendreBasis.hpp"
#include "Stokhos_TotalOrderBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"

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
test_mean_multiply(const OrdinalType order,
                   const OrdinalType dim,
                   const OrdinalType nGrid,
                   const OrdinalType iterCount,
                   std::vector<double>& scalar_perf,
                   std::vector<double>& block_left_perf,
                   std::vector<double>& block_right_perf,
                   std::vector<double>& pce_perf,
                   std::vector<double>& block_pce_perf)
{
  typedef ScalarType value_type;
  typedef OrdinalType ordinal_type;
  typedef Device execution_space;

  typedef Stokhos::DynamicStorage<ordinal_type,value_type,execution_space> storage_type;
  typedef Sacado::UQ::PCE<storage_type> pce_type;

  typedef Kokkos::View< value_type*, Kokkos::LayoutLeft, execution_space > scalar_vector_type;
  typedef Kokkos::View< value_type**, Kokkos::LayoutLeft, execution_space > scalar_left_multi_vector_type;
   typedef Kokkos::View< value_type**, Kokkos::LayoutRight, execution_space > scalar_right_multi_vector_type;
  typedef Kokkos::View< pce_type*, Kokkos::LayoutLeft, execution_space > pce_vector_type;
  typedef Kokkos::View< pce_type**, Kokkos::LayoutLeft, execution_space > pce_multi_vector_type;

  typedef KokkosSparse::CrsMatrix< value_type, ordinal_type, execution_space > scalar_matrix_type;
  typedef KokkosSparse::CrsMatrix< pce_type, ordinal_type, execution_space > pce_matrix_type;
  typedef typename scalar_matrix_type::StaticCrsGraphType matrix_graph_type;
  typedef typename scalar_matrix_type::values_type scalar_matrix_values_type;
  typedef typename pce_matrix_type::values_type pce_matrix_values_type;

  typedef Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> abstract_basis_type;
  typedef Stokhos::LegendreBasis<ordinal_type,value_type> basis_type;
  typedef Stokhos::LexographicLess<Stokhos::MultiIndex<ordinal_type> > order_type;
  typedef Stokhos::TotalOrderBasis<ordinal_type,value_type,order_type> product_basis_type;
  typedef Stokhos::Sparse3Tensor<ordinal_type,value_type> cijk_type;
  typedef typename pce_type::cijk_type kokkos_cijk_type;

  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::Array;

  // Number of columns for PCE multi-vector apply
  const ordinal_type num_pce_col = 5;

  // Create Stochastic Galerkin basis and expansion
  Array< RCP<const abstract_basis_type> > bases(dim);
  for (ordinal_type i=0; i<dim; ++i) {
    bases[i] = Teuchos::rcp(new basis_type(order, true));
  }
  RCP<const product_basis_type> basis = rcp(new product_basis_type(bases));
  RCP<cijk_type> cijk = basis->computeTripleProductTensor();
  kokkos_cijk_type kokkos_cijk =
    Stokhos::create_product_tensor<execution_space>(*basis, *cijk);
  Kokkos::setGlobalCijkTensor(kokkos_cijk);

  //------------------------------
  // Generate graph for "FEM" box structure:

  std::vector< std::vector<size_t> > fem_graph;
  const size_t fem_length = nGrid * nGrid * nGrid;
  const size_t graph_length = generate_fem_graph( nGrid , fem_graph );

  //------------------------------
  // Generate input vectors:

  ordinal_type pce_size = basis->size();
  scalar_left_multi_vector_type xl(Kokkos::ViewAllocateWithoutInitializing("scalar left x"), fem_length, pce_size);
  scalar_left_multi_vector_type yl(Kokkos::ViewAllocateWithoutInitializing("scalar right y"), fem_length, pce_size);
  scalar_right_multi_vector_type xr(Kokkos::ViewAllocateWithoutInitializing("scalar right x"), fem_length, pce_size);
  scalar_right_multi_vector_type yr(Kokkos::ViewAllocateWithoutInitializing("scalar right y"), fem_length, pce_size);
  std::vector<scalar_vector_type> x_col(pce_size), y_col(pce_size);
  for (ordinal_type i=0; i<pce_size; ++i) {
    x_col[i] = scalar_vector_type (Kokkos::ViewAllocateWithoutInitializing("scalar x col"), fem_length);
    y_col[i] = scalar_vector_type(Kokkos::ViewAllocateWithoutInitializing("scalar y col"), fem_length);
    Kokkos::deep_copy( x_col[i] , value_type(1.0) );
    Kokkos::deep_copy( y_col[i] , value_type(0.0) );
  }
  pce_vector_type x_pce =
    Kokkos::make_view<pce_vector_type>(Kokkos::ViewAllocateWithoutInitializing("pce x"),
                                       kokkos_cijk, fem_length, pce_size);
  pce_vector_type y_pce =
    Kokkos::make_view<pce_vector_type>(Kokkos::ViewAllocateWithoutInitializing("pce y"),
                                       kokkos_cijk, fem_length, pce_size);
  pce_multi_vector_type x_multi_pce =
    Kokkos::make_view<pce_multi_vector_type>(Kokkos::ViewAllocateWithoutInitializing("pce multi x"),
                                             kokkos_cijk, fem_length,
                                             num_pce_col, pce_size);
  pce_multi_vector_type y_multi_pce =
    Kokkos::make_view<pce_multi_vector_type>(Kokkos::ViewAllocateWithoutInitializing("pce multi y"),
                                             kokkos_cijk, fem_length,
                                             num_pce_col, pce_size);

  Kokkos::deep_copy( xl , value_type(1.0) );
  Kokkos::deep_copy( yl , value_type(0.0) );
  Kokkos::deep_copy( xr , value_type(1.0) );
  Kokkos::deep_copy( yr , value_type(0.0) );
  Kokkos::deep_copy( x_pce , value_type(1.0) );
  Kokkos::deep_copy( y_pce , value_type(0.0) );
  Kokkos::deep_copy( x_multi_pce , value_type(1.0) );
  Kokkos::deep_copy( y_multi_pce , value_type(0.0) );

  //------------------------------
  // Generate matrix

  matrix_graph_type matrix_graph =
    Kokkos::create_staticcrsgraph<matrix_graph_type>(
      std::string("test crs graph"), fem_graph);
  scalar_matrix_values_type scalar_matrix_values =
    scalar_matrix_values_type(Kokkos::ViewAllocateWithoutInitializing("scalar matrix"), graph_length);
  pce_matrix_values_type pce_matrix_values =
    Kokkos::make_view<pce_matrix_values_type>(Kokkos::ViewAllocateWithoutInitializing("pce matrix"), kokkos_cijk, graph_length, 1);
  scalar_matrix_type scalar_matrix("scalar matrix", fem_length,
                                   scalar_matrix_values, matrix_graph);
  pce_matrix_type pce_matrix("pce matrix", fem_length,
                             pce_matrix_values, matrix_graph);

  Kokkos::deep_copy( scalar_matrix_values , value_type(1.0) );
  Kokkos::deep_copy( pce_matrix_values , value_type(1.0) );

  //------------------------------
  // Scalar multiply

  {
    // warm up
    for (ordinal_type iter = 0; iter < iterCount; ++iter) {
      for (ordinal_type col=0; col<pce_size; ++col) {
        // scalar_vector_type xc =
        //   Kokkos::subview(x, Kokkos::ALL(), col);
        // scalar_vector_type yc =
        //   Kokkos::subview(y, Kokkos::ALL(), col);
        // Kokkos::MV_Multiply( yc, scalar_matrix, xc );
        KokkosSparse::spmv(  "N" , value_type(1.0) , scalar_matrix, x_col[col] , value_type(0.0) ,y_col[col]);
      }
    }

    execution_space().fence();
    Kokkos::Timer clock ;
    for (ordinal_type iter = 0; iter < iterCount; ++iter) {
      for (ordinal_type col=0; col<pce_size; ++col) {
        // scalar_vector_type xc =
        //   Kokkos::subview(x, Kokkos::ALL(), col);
        // scalar_vector_type yc =
        //   Kokkos::subview(y, Kokkos::ALL(), col);
        // Kokkos::MV_Multiply( yc, scalar_matrix, xc );
        KokkosSparse::spmv(  "N" , value_type(1.0) , scalar_matrix, x_col[col] , value_type(0.0) ,y_col[col]);
      }
    }
    execution_space().fence();

    const double seconds_per_iter = clock.seconds() / ((double) iterCount );
    const double flops = 1.0e-9 * 2.0 * graph_length * pce_size;

    scalar_perf.resize(5);
    scalar_perf[0] = fem_length;
    scalar_perf[1] = pce_size;
    scalar_perf[2] = graph_length;
    scalar_perf[3] = seconds_per_iter;
    scalar_perf[4] = flops / seconds_per_iter;
  }

  //------------------------------
  // Block-left multiply

  {
    // warm up
    for (ordinal_type iter = 0; iter < iterCount; ++iter) {
      KokkosSparse::spmv(  "N" , value_type(1.0) , scalar_matrix, xl , value_type(0.0) ,yl);
    }

    execution_space().fence();
    Kokkos::Timer clock ;
    for (ordinal_type iter = 0; iter < iterCount; ++iter) {
      KokkosSparse::spmv(  "N" , value_type(1.0) , scalar_matrix, xl , value_type(0.0) ,yl);
    }
    execution_space().fence();

    const double seconds_per_iter = clock.seconds() / ((double) iterCount );
    const double flops = 1.0e-9 * 2.0 * graph_length * pce_size;

    block_left_perf.resize(5);
    block_left_perf[0] = fem_length;
    block_left_perf[1] = pce_size;
    block_left_perf[2] = graph_length;
    block_left_perf[3] = seconds_per_iter;
    block_left_perf[4] = flops / seconds_per_iter;
  }

  //------------------------------
  // Block-right multiply

  {
    // warm up
    for (ordinal_type iter = 0; iter < iterCount; ++iter) {
      KokkosSparse::spmv(  "N" , value_type(1.0) , scalar_matrix, xr , value_type(0.0) ,yr);
    }

    execution_space().fence();
    Kokkos::Timer clock ;
    for (ordinal_type iter = 0; iter < iterCount; ++iter) {
      KokkosSparse::spmv(  "N" , value_type(1.0) , scalar_matrix, xr , value_type(0.0) ,yr);
    }
    execution_space().fence();

    const double seconds_per_iter = clock.seconds() / ((double) iterCount );
    const double flops = 1.0e-9 * 2.0 * graph_length * pce_size;

    block_right_perf.resize(5);
    block_right_perf[0] = fem_length;
    block_right_perf[1] = pce_size;
    block_right_perf[2] = graph_length;
    block_right_perf[3] = seconds_per_iter;
    block_right_perf[4] = flops / seconds_per_iter;
  }

  //------------------------------
  // PCE multiply

  {
    // warm up
    for (ordinal_type iter = 0; iter < iterCount; ++iter) {
      KokkosSparse::spmv(  "N" , value_type(1.0) , pce_matrix, x_pce , value_type(0.0) ,y_pce);
    }

    execution_space().fence();
    Kokkos::Timer clock ;
    for (ordinal_type iter = 0; iter < iterCount; ++iter) {
      KokkosSparse::spmv(  "N" , value_type(1.0) , pce_matrix, x_pce , value_type(0.0) ,y_pce);
    }
    execution_space().fence();

    const double seconds_per_iter = clock.seconds() / ((double) iterCount );
    const double flops = 1.0e-9 * 2.0 * graph_length * pce_size;

    pce_perf.resize(5);
    pce_perf[0] = fem_length;
    pce_perf[1] = pce_size;
    pce_perf[2] = graph_length;
    pce_perf[3] = seconds_per_iter;
    pce_perf[4] = flops / seconds_per_iter;
  }

  //------------------------------
  // PCE multi-vector multiply

  {
    // warm up
    for (ordinal_type iter = 0; iter < iterCount; ++iter) {
      KokkosSparse::spmv(  "N" , value_type(1.0) , pce_matrix, x_multi_pce , value_type(0.0) ,y_multi_pce);
    }

    execution_space().fence();
    Kokkos::Timer clock ;
    for (ordinal_type iter = 0; iter < iterCount; ++iter) {
      KokkosSparse::spmv(  "N" , value_type(1.0) , pce_matrix, x_multi_pce , value_type(0.0) ,y_multi_pce);
    }
    execution_space().fence();

    const double seconds_per_iter = clock.seconds() / ((double) iterCount );
    const double flops = 1.0e-9 * 2.0 * graph_length * pce_size * num_pce_col;

    block_pce_perf.resize(5);
    block_pce_perf[0] = fem_length;
    block_pce_perf[1] = pce_size;
    block_pce_perf[2] = graph_length;
    block_pce_perf[3] = seconds_per_iter;
    block_pce_perf[4] = flops / seconds_per_iter;
  }

}

template <typename Scalar, typename Ordinal, typename Device>
void performance_test_driver( const Ordinal nGrid,
                              const Ordinal nIter,
                              const Ordinal order,
                              const Ordinal min_var,
                              const Ordinal max_var )
{
  std::cout.precision(8);
  std::cout << std::endl
            << "\"Grid Size\" , "
            << "\"FEM Size\" , "
            << "\"FEM Graph Size\" , "
            << "\"Dimension\" , "
            << "\"Order\" , "
            << "\"PCE Size\" , "
            << "\"Scalar SpMM Time\" , "
            << "\"Scalar SpMM Speedup\" , "
            << "\"Scalar SpMM GFLOPS\" , "
            << "\"Block-Left SpMM Speedup\" , "
            << "\"Block-Left SpMM GFLOPS\" , "
            << "\"Block-Right SpMM Speedup\" , "
            << "\"Block-Right SpMM GFLOPS\" , "
            << "\"PCE SpMM Speedup\" , "
            << "\"PCE SpMM GFLOPS\" , "
            << "\"Block PCE SpMM Speedup\" , "
            << "\"Block PCE SpMM GFLOPS\" , "
            << std::endl;

  std::vector<double> perf_scalar, perf_block_left, perf_block_right,
    perf_pce, perf_block_pce;
  for (Ordinal dim=min_var; dim<=max_var; ++dim) {

    test_mean_multiply<Scalar,Ordinal,Device>(
      order, dim, nGrid, nIter, perf_scalar, perf_block_left, perf_block_right,
      perf_pce, perf_block_pce );

    std::cout << nGrid << " , "
              << perf_scalar[0] << " , "
              << perf_scalar[2] << " , "
              << dim << " , "
              << order << " , "
              << perf_scalar[1] << " , "
              << perf_scalar[3] << " , "
              << perf_scalar[4] / perf_scalar[4] << " , "
              << perf_scalar[4] << " , "
              << perf_block_left[4]/ perf_scalar[4] << " , "
              << perf_block_left[4] << " , "
              << perf_block_right[4]/ perf_scalar[4] << " , "
              << perf_block_right[4] << " , "
              << perf_pce[4]/ perf_scalar[4] << " , "
              << perf_pce[4] << " , "
              << perf_block_pce[4]/ perf_scalar[4] << " , "
              << perf_block_pce[4] << " , "
              << std::endl;

  }
}

#define INST_PERF_DRIVER(SCALAR, ORDINAL, DEVICE)                       \
  template void performance_test_driver< SCALAR, ORDINAL, DEVICE >(     \
    const ORDINAL nGrid, const ORDINAL nIter,  const ORDINAL order,     \
    const ORDINAL min_var, const ORDINAL max_var);
