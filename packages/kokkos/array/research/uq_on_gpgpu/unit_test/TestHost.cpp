
// #include <Kokkos_ProductTensor.hpp>

#include <Kokkos_SymmetricDiagonalSpec.hpp>
#include <Kokkos_BlockCrsMatrix.hpp>

//

#include <Kokkos_Host.hpp>

#include <Kokkos_Host_macros.hpp>
// #include <impl/Kokkos_ProductTensor_macros.hpp>
#include <impl/Kokkos_SymmetricDiagonalSpec_macros.hpp>
#include <impl/Kokkos_BlockCrsMatrix_macros.hpp>
#include <Kokkos_Clear_macros.hpp>

//

// #include <TestSparseProductTensor.hpp>
#include <TestBlockCrsMatrix.hpp>

int mainHost()
{
  Kokkos::Host::initialize( Kokkos::Host::SetThreadCount(4) );

//  unit_test::test_dense<Kokkos::Host>();
//  unit_test::test_diagonal<Kokkos::Host>();
//  unit_test::test_other<Kokkos::Host>();

  unit_test::test_block_crs_matrix<Kokkos::Host>( 1 , 2 );
  unit_test::test_block_crs_matrix<Kokkos::Host>( 1 , 5 );
  unit_test::test_block_crs_matrix<Kokkos::Host>( 2 , 1 );
  unit_test::test_block_crs_matrix<Kokkos::Host>( 3 , 1 );

  unit_test::test_block_crs_matrix<Kokkos::Host>( 10 , 8 );
  unit_test::test_block_crs_matrix<Kokkos::Host>( 11 , 8 );
  unit_test::test_block_crs_matrix<Kokkos::Host>( 12 , 10 );
  unit_test::test_block_crs_matrix<Kokkos::Host>( 13 , 10 );

  Kokkos::Host::finalize();

  return 0 ;
}

