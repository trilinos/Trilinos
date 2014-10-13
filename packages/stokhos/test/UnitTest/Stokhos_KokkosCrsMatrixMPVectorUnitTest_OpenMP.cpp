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
    Kokkos::DeviceConfig dev_config(num_cores, ThreadsPerVector, row_threads);

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
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Initialize threads
  num_cores =
    Kokkos::hwloc::get_available_numa_count() *
    Kokkos::hwloc::get_available_cores_per_numa();
  num_hyper_threads =
    Kokkos::hwloc::get_available_threads_per_core();
  Kokkos::OpenMP::initialize(num_cores * num_hyper_threads);
  //Kokkos::OpenMP::print_configuration(std::cout);

  // Run tests
  int ret = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finish up
  Kokkos::OpenMP::finalize();

  return ret;
}
