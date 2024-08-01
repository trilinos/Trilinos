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

#include "Kokkos_Core.hpp"

// Instantiate test for OpenMP device
using Kokkos::OpenMP;
CRSMATRIX_MP_VECTOR_TESTS_DEVICE( OpenMP )

template <typename Storage, typename Ordinal, typename MultiplyOp,
          Ordinal NumPerThread, Ordinal ThreadsPerVector>
bool test_host_embedded_vector(Ordinal num_hyper_threads,
                               Ordinal num_cores,
                               Teuchos::FancyOStream& out) {
  const Ordinal VectorSize = NumPerThread * ThreadsPerVector;
  typedef typename Storage::template apply_N<VectorSize>::type storage_type;
  typedef Sacado::MP::Vector<storage_type> Vector;

  const Ordinal nGrid = 5;

  bool success = true;
  if (num_hyper_threads >= ThreadsPerVector) {
    int row_threads = num_hyper_threads / ThreadsPerVector;
    KokkosSparse::DeviceConfig dev_config(num_cores, ThreadsPerVector, row_threads);

    success = test_embedded_vector<Vector>(
      nGrid, VectorSize, dev_config, MultiplyOp(), out);
  }
  return success;
}

size_t num_cores, num_hyper_threads;

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(
  Kokkos_CrsMatrix_MP, Multiply_1, Storage, MultiplyOp )
{
  typedef typename Storage::ordinal_type Ordinal;
  const Ordinal NumPerThread = 16;
  const Ordinal ThreadsPerVector = 1;
  success =
    test_host_embedded_vector<Storage,Ordinal,MultiplyOp,NumPerThread,ThreadsPerVector>(num_hyper_threads, num_cores, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(
  Kokkos_CrsMatrix_MP, Multiply_2, Storage, MultiplyOp )
{
  typedef typename Storage::ordinal_type Ordinal;
  const Ordinal NumPerThread = 8;
  const Ordinal ThreadsPerVector = 2;
  success =
    test_host_embedded_vector<Storage,Ordinal,MultiplyOp,NumPerThread,ThreadsPerVector>(num_hyper_threads, num_cores, out);
}

#define CRS_MATRIX_MP_VECTOR_MULTIPLY_TESTS_STORAGE_OP( STORAGE, OP )   \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_CrsMatrix_MP, Multiply_1, STORAGE, OP )                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_CrsMatrix_MP, Multiply_2, STORAGE, OP )

// Notes:  SFS, DS are defined in main test header (we are also being lazy
// and not putting ordinal/scalar/device in the names, assuming we will only
// do one combination).  We can't do DefaultMultiply for DS because it
// uses partitioning
#define CRS_MATRIX_MP_VECTOR_MULTIPLY_TESTS_ORDINAL_SCALAR_DEVICE( ORDINAL, SCALAR, DEVICE ) \
  CRS_MATRIX_MP_VECTOR_MULTIPLY_TESTS_STORAGE_OP( SFS, DefaultMultiply ) \
  CRS_MATRIX_MP_VECTOR_MULTIPLY_TESTS_STORAGE_OP( SFS, KokkosMultiply ) \
  CRS_MATRIX_MP_VECTOR_MULTIPLY_TESTS_STORAGE_OP( DS, DefaultMultiply ) \
  CRS_MATRIX_MP_VECTOR_MULTIPLY_TESTS_STORAGE_OP( DS, KokkosMultiply )

CRS_MATRIX_MP_VECTOR_MULTIPLY_TESTS_ORDINAL_SCALAR_DEVICE(int, double, OpenMP)

int main( int argc, char* argv[] ) {
  // Setup the MPI session
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Initialize Kokkos
  Kokkos::initialize(argc, argv);

  // Run tests
  int ret = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finish up
  Kokkos::finalize();

  return ret;
}
