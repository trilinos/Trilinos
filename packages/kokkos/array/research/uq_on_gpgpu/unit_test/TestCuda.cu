
// #include <Kokkos_ProductTensor.hpp>

#include <Kokkos_SymmetricDiagonalSpec.hpp>
#include <Kokkos_BlockCrsMatrix.hpp>

//

#include <Kokkos_Host.hpp>
#include <Kokkos_Cuda.hpp>
#include <Cuda/Kokkos_Cuda_SymmetricDiagonalSpec.hpp>
#include <Cuda/Kokkos_Cuda_BlockCrsMatrix.hpp>

//

// #include <TestSparseProductTensor.hpp>
#include <TestBlockCrsMatrix.hpp>

int mainCuda()
{
  Kokkos::Cuda::initialize();

//  unit_test::test_dense<Kokkos::Cuda>();
//  unit_test::test_diagonal<Kokkos::Cuda>();
//  unit_test::test_other<Kokkos::Cuda>();

  unit_test::test_block_crs_matrix<Kokkos::Cuda>( 1 , 2 );
  unit_test::test_block_crs_matrix<Kokkos::Cuda>( 1 , 5 );
  unit_test::test_block_crs_matrix<Kokkos::Cuda>( 2 , 1 );
  unit_test::test_block_crs_matrix<Kokkos::Cuda>( 3 , 1 );

  unit_test::test_block_crs_matrix<Kokkos::Cuda>( 10 , 8 );
  unit_test::test_block_crs_matrix<Kokkos::Cuda>( 11 , 8 );
  unit_test::test_block_crs_matrix<Kokkos::Cuda>( 12 , 10 );
  unit_test::test_block_crs_matrix<Kokkos::Cuda>( 13 , 10 );

  Kokkos::Cuda::finalize();

  return 0 ;
}

