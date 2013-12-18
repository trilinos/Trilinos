// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Stokhos_KokkosMPVectorKernelsUnitTest.hpp"
#include "Kokkos_CrsMatrix_MP_Vector_Cuda.hpp"

#include "Kokkos_Cuda.hpp"

using namespace KokkosMPVectorKernelsUnitTest;

template <typename Ordinal, typename Scalar, typename MultiplyTag,
          Ordinal NumPerThread, Ordinal ThreadsPerVector>
bool test_cuda_static_fixed_embedded_vector(Ordinal num_blocks,
                                            Ordinal num_vec_threads,
                                            Ordinal num_row_threads,
                                            Teuchos::FancyOStream& out) {
  typedef Kokkos::Cuda Device;

  const Ordinal VectorSize = NumPerThread * ThreadsPerVector;
  typedef Stokhos::StaticFixedStorage<Ordinal,Scalar,VectorSize,Device> Storage;
  typedef Sacado::MP::Vector<Storage> Vector;

  const Ordinal nGrid = 5;
  Kokkos::DeviceConfig dev_config(num_blocks, num_vec_threads, num_row_threads);

  bool success = test_embedded_vector<Vector>(
    nGrid, VectorSize, dev_config, MultiplyTag(), out);

  return success;
}

// Test default configuration
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_CrsMatrix_MP, Multiply_Default, Ordinal, Scalar, MultiplyTag )
{
  const Ordinal NumPerThread = 1;
  const Ordinal ThreadsPerVector = 16;

  const Ordinal num_blocks = 0;
  const Ordinal num_vec_threads = 0;
  const Ordinal num_row_threads = 0;

  success =
    test_cuda_static_fixed_embedded_vector<Ordinal,Scalar,MultiplyTag,NumPerThread,ThreadsPerVector>(num_blocks, num_vec_threads, num_row_threads, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_CrsMatrix_MP, Multiply_1, Ordinal, Scalar, MultiplyTag )
{
  const Ordinal NumPerThread = 1;
  const Ordinal ThreadsPerVector = 16;

  const Ordinal num_blocks = 10;
  const Ordinal num_vec_threads = ThreadsPerVector;
  const Ordinal num_row_threads = 4;

  success =
    test_cuda_static_fixed_embedded_vector<Ordinal,Scalar,MultiplyTag,NumPerThread,ThreadsPerVector>(num_blocks, num_vec_threads, num_row_threads, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_CrsMatrix_MP, Multiply_2, Ordinal, Scalar, MultiplyTag )
{
  const Ordinal NumPerThread = 2;
  const Ordinal ThreadsPerVector = 16;

  const Ordinal num_blocks = 10;
  const Ordinal num_vec_threads = ThreadsPerVector;
  const Ordinal num_row_threads = 4;

  success =
    test_cuda_static_fixed_embedded_vector<Ordinal,Scalar,MultiplyTag,NumPerThread,ThreadsPerVector>(num_blocks, num_vec_threads, num_row_threads, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_CrsMatrix_MP, Multiply_3, Ordinal, Scalar, MultiplyTag )
{
  const Ordinal NumPerThread = 3;
  const Ordinal ThreadsPerVector = 16;

  const Ordinal num_blocks = 10;
  const Ordinal num_vec_threads = ThreadsPerVector;
  const Ordinal num_row_threads = 4;

  success =
    test_cuda_static_fixed_embedded_vector<Ordinal,Scalar,MultiplyTag,NumPerThread,ThreadsPerVector>(num_blocks, num_vec_threads, num_row_threads, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_CrsMatrix_MP, Multiply_4, Ordinal, Scalar, MultiplyTag )
{
  const Ordinal NumPerThread = 4;
  const Ordinal ThreadsPerVector = 16;

  const Ordinal num_blocks = 10;
  const Ordinal num_vec_threads = ThreadsPerVector;
  const Ordinal num_row_threads = 4;

  success =
    test_cuda_static_fixed_embedded_vector<Ordinal,Scalar,MultiplyTag,NumPerThread,ThreadsPerVector>(num_blocks, num_vec_threads, num_row_threads, out);
}

#define CRS_MATRIX_MP_VECTOR_TESTS_ORDINAL_SCALAR_TAG( ORDINAL, SCALAR, TAG ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(                                 \
    Kokkos_CrsMatrix_MP, Multiply_Default,  ORDINAL, SCALAR, TAG )      \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(                                 \
    Kokkos_CrsMatrix_MP, Multiply_1,  ORDINAL, SCALAR, TAG )            \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(                                 \
    Kokkos_CrsMatrix_MP, Multiply_2,  ORDINAL, SCALAR, TAG )            \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(                                 \
    Kokkos_CrsMatrix_MP, Multiply_3,  ORDINAL, SCALAR, TAG )            \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(                                 \
    Kokkos_CrsMatrix_MP, Multiply_4,  ORDINAL, SCALAR, TAG )

using Stokhos::DefaultMultiply;
using Stokhos::EnsembleMultiply;
CRS_MATRIX_MP_VECTOR_TESTS_ORDINAL_SCALAR_TAG(int, double, DefaultMultiply)
CRS_MATRIX_MP_VECTOR_TESTS_ORDINAL_SCALAR_TAG(int, double, EnsembleMultiply)

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Initialize Cuda
  Kokkos::Cuda::host_mirror_device_type::initialize();
  Kokkos::Cuda::initialize(Kokkos::Cuda::SelectDevice(0));
  Kokkos::Cuda::print_configuration(std::cout);

  // Run tests
  int ret = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finish up
  Kokkos::Cuda::host_mirror_device_type::finalize();
  Kokkos::Cuda::finalize();

  return ret;
}
