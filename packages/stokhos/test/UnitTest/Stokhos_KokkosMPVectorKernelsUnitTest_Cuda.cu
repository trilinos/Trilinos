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

typedef UnitTestSetup<int,double,Kokkos::Cuda> SetupType;

TEUCHOS_UNIT_TEST( Stokhos_KokkosMPVectorKernels,
                   EmbeddedVector_MPKernel_Cuda ) {
  typedef int Ordinal;
  typedef double Scalar;
  typedef Kokkos::Cuda Device;
  const Ordinal VectorSize = SetupType::local_vec_size;
  typedef Stokhos::StaticFixedStorage<Ordinal,Scalar,VectorSize,Device> Storage;
  typedef Sacado::MP::Vector<Storage> Vector;

  int vector_threads = 16;
  int row_threads = 4;
  int num_cores = 10;
  Kokkos::DeviceConfig dev_config(num_cores, vector_threads, row_threads);

  SetupType setup;
  setup.setup(vector_threads);
  success = test_embedded_vector<Vector>(
    setup, dev_config, Stokhos::DefaultMultiply(), out);
  setup.cleanup();
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosMPVectorKernels,
                   EmbeddedVector_Ensemble_Cuda ) {
  typedef int Ordinal;
  typedef double Scalar;
  typedef Kokkos::Cuda Device;
  const Ordinal VectorSize = SetupType::local_vec_size;
  typedef Stokhos::StaticFixedStorage<Ordinal,Scalar,VectorSize,Device> Storage;
  typedef Sacado::MP::Vector<Storage> Vector;

  int vector_threads = 16;
  int row_threads = 4;
  int num_cores = 10;
  Kokkos::DeviceConfig dev_config(num_cores, vector_threads, row_threads);

  SetupType setup;
  setup.setup(vector_threads);
  success = test_embedded_vector<Vector>(
    setup, dev_config, Stokhos::EnsembleMultiply(), out);
  setup.cleanup();
}

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
