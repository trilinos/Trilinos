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

// Utilities
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

// Device
#include "Kokkos_Core.hpp"

// Kernels
#include "Stokhos_ConfigDefs.h"
#include "Stokhos_OpenMP_CrsProductTensor.hpp"
#ifdef HAVE_STOKHOS_MKL
#include "Stokhos_OpenMP_MKL_CrsMatrix.hpp"
#endif

// Tests
#include "Stokhos_KokkosArrayKernelsUnitTest.hpp"

using namespace KokkosKernelsUnitTest;

UnitTestSetup<Kokkos::OpenMP> setup;

// Test declarations
#include "Stokhos_KokkosArrayKernelsUnitTestDecl.hpp"
#include "Stokhos_KokkosArrayKernelsUnitTest_Host.hpp"

// Tests using OpenMP device
using Kokkos::OpenMP;
UNIT_TEST_GROUP_SCALAR_DEVICE( double, OpenMP )
UNIT_TEST_GROUP_SCALAR_HOST_DEVICE( double, OpenMP )

#ifdef HAVE_STOKHOS_MKL
TEUCHOS_UNIT_TEST( Kokkos_SG_SpMv, double_OpenMP_CrsMatrixFree_MKL ) {
  typedef double Scalar;
  typedef Kokkos::OpenMP Device;
  typedef Stokhos::MKLMultiply SparseMatOps;
  success = test_crs_matrix_free<Scalar,Device,SparseMatOps>(
    setup, out);
}
#endif

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  const size_t team_count =
  Kokkos::hwloc::get_available_numa_count() *
    Kokkos::hwloc::get_available_cores_per_numa();
  const size_t threads_per_team =
    Kokkos::hwloc::get_available_threads_per_core();

  // Initialize openmp
  Kokkos::OpenMP::initialize( team_count * threads_per_team );
  //Kokkos::OpenMP::print_configuration( std::cout );

  // Setup (has to happen after initialization)
  setup.setup();

  // Run tests
  int ret = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finish up
  Kokkos::OpenMP::finalize();

  return ret;
}
