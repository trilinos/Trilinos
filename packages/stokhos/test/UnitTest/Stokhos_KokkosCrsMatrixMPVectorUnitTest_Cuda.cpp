// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Stokhos_KokkosCrsMatrixMPVectorUnitTest.hpp"

// Instantiate test for Cuda device
using Kokkos::Cuda;
CRSMATRIX_MP_VECTOR_TESTS_DEVICE( Cuda )

template <typename Storage, typename Ordinal, typename MultiplyOp,
          Ordinal NumPerThread, Ordinal ThreadsPerVector>
bool test_cuda_embedded_vector(Ordinal num_blocks,
                               Ordinal num_vec_threads,
                               Ordinal num_row_threads,
                               Teuchos::FancyOStream& out) {
  typedef Kokkos::Cuda Device;

  const Ordinal VectorSize = NumPerThread * ThreadsPerVector;
  typedef typename Storage::template apply_N<VectorSize>::type storage_type;
  typedef Sacado::MP::Vector<storage_type> Vector;

  const Ordinal nGrid = 5;
  KokkosSparse::DeviceConfig dev_config(num_blocks, num_vec_threads, num_row_threads);

  bool success = test_embedded_vector<Vector>(
    nGrid, VectorSize, dev_config, MultiplyOp(), out);

  return success;
}

// Test default configuration
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(
  Kokkos_CrsMatrix_MP, Multiply_Default, Storage, MultiplyOp )
{
  typedef typename Storage::ordinal_type Ordinal;
  const Ordinal NumPerThread = 1;
  const Ordinal ThreadsPerVector = 16;

  const Ordinal num_blocks = 0;
  const Ordinal num_vec_threads = 0;
  const Ordinal num_row_threads = 0;

  success =
    test_cuda_embedded_vector<Storage,Ordinal,MultiplyOp,NumPerThread,ThreadsPerVector>(num_blocks, num_vec_threads, num_row_threads, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(
  Kokkos_CrsMatrix_MP, Multiply_1, Storage, MultiplyOp )
{
  typedef typename Storage::ordinal_type Ordinal;
  const Ordinal NumPerThread = 1;
  const Ordinal ThreadsPerVector = 16;

  const Ordinal num_blocks = 10;
  const Ordinal num_vec_threads = ThreadsPerVector;
  const Ordinal num_row_threads = 4;

  success =
    test_cuda_embedded_vector<Storage,Ordinal,MultiplyOp,NumPerThread,ThreadsPerVector>(num_blocks, num_vec_threads, num_row_threads, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(
  Kokkos_CrsMatrix_MP, Multiply_2, Storage, MultiplyOp )
{
  typedef typename Storage::ordinal_type Ordinal;
  const Ordinal NumPerThread = 2;
  const Ordinal ThreadsPerVector = 16;

  const Ordinal num_blocks = 10;
  const Ordinal num_vec_threads = ThreadsPerVector;
  const Ordinal num_row_threads = 4;

  success =
    test_cuda_embedded_vector<Storage,Ordinal,MultiplyOp,NumPerThread,ThreadsPerVector>(num_blocks, num_vec_threads, num_row_threads, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(
  Kokkos_CrsMatrix_MP, Multiply_3, Storage, MultiplyOp )
{
  typedef typename Storage::ordinal_type Ordinal;
  const Ordinal NumPerThread = 3;
  const Ordinal ThreadsPerVector = 16;

  const Ordinal num_blocks = 10;
  const Ordinal num_vec_threads = ThreadsPerVector;
  const Ordinal num_row_threads = 4;

  success =
    test_cuda_embedded_vector<Storage,Ordinal,MultiplyOp,NumPerThread,ThreadsPerVector>(num_blocks, num_vec_threads, num_row_threads, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(
  Kokkos_CrsMatrix_MP, Multiply_4, Storage, MultiplyOp )
{
  typedef typename Storage::ordinal_type Ordinal;
  const Ordinal NumPerThread = 4;
  const Ordinal ThreadsPerVector = 16;

  const Ordinal num_blocks = 10;
  const Ordinal num_vec_threads = ThreadsPerVector;
  const Ordinal num_row_threads = 4;

  success =
    test_cuda_embedded_vector<Storage,Ordinal,MultiplyOp,NumPerThread,ThreadsPerVector>(num_blocks, num_vec_threads, num_row_threads, out);
}

#define CRS_MATRIX_MP_VECTOR_MULTIPLY_TESTS_STORAGE_OP( STORAGE, OP ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                               \
    Kokkos_CrsMatrix_MP, Multiply_Default, STORAGE, OP )                \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_CrsMatrix_MP, Multiply_1, STORAGE, OP )                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_CrsMatrix_MP, Multiply_2, STORAGE, OP )                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_CrsMatrix_MP, Multiply_3, STORAGE, OP )                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_CrsMatrix_MP, Multiply_4, STORAGE, OP )

// Notes:  SFS, DS are defined in main test header (we are also being lazy
// and not putting ordinal/scalar/device in the names, assuming we will only
// do one combination).
#define CRS_MATRIX_MP_VECTOR_MULTIPLY_TESTS_ORDINAL_SCALAR_DEVICE( ORDINAL, SCALAR, DEVICE ) \
  CRS_MATRIX_MP_VECTOR_MULTIPLY_TESTS_STORAGE_OP( SFS, DefaultMultiply ) \
  CRS_MATRIX_MP_VECTOR_MULTIPLY_TESTS_STORAGE_OP( SFS, KokkosMultiply ) \
  CRS_MATRIX_MP_VECTOR_MULTIPLY_TESTS_STORAGE_OP( DS, DefaultMultiply ) \
  CRS_MATRIX_MP_VECTOR_MULTIPLY_TESTS_STORAGE_OP( DS, KokkosMultiply )

CRS_MATRIX_MP_VECTOR_MULTIPLY_TESTS_ORDINAL_SCALAR_DEVICE(int, double, Cuda)

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Initialize Cuda
  Kokkos::InitializationSettings init_args;
  init_args.set_device_id(0);
  Kokkos::initialize( init_args );
  Kokkos::print_configuration(std::cout);

  // Run tests
  int ret = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finish up
  Kokkos::finalize();

  return ret;
}
