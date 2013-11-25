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

#include "Kokkos_Cuda.hpp"

using namespace KokkosMPVectorKernelsUnitTest;

typedef UnitTestSetup<int,double,Kokkos::Cuda> SetupType;
SetupType setup;

TEUCHOS_UNIT_TEST( Stokhos_KokkosMPVectorKernels, EmbeddedVector_Left_Cuda ) {
  typedef int Ordinal;
  typedef double Scalar;
  typedef Kokkos::Cuda Device;
  typedef Kokkos::LayoutLeft Layout;
  const Ordinal VectorSize = SetupType::local_vec_size;
  typedef Stokhos::StaticFixedStorage<Ordinal,Scalar,VectorSize,Device> Storage;
  typedef Sacado::MP::Vector<Storage> Vector;

  success = test_embedded_vector<Vector,Layout>(setup, out);
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosMPVectorKernels, EmbeddedVector_Right_Cuda ) {
  typedef int Ordinal;
  typedef double Scalar;
  typedef Kokkos::Cuda Device;
  typedef Kokkos::LayoutRight Layout;
  const Ordinal VectorSize = SetupType::local_vec_size;
  typedef Stokhos::StaticFixedStorage<Ordinal,Scalar,VectorSize,Device> Storage;
  typedef Sacado::MP::Vector<Storage> Vector;

  success = test_embedded_vector<Vector,Layout>(setup, out);
}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Initialize Cuda
  Kokkos::Cuda::host_mirror_device_type::initialize();
  Kokkos::Cuda::initialize(Kokkos::Cuda::SelectDevice(0));
  Kokkos::Cuda::print_configuration(std::cout);

  // Setup
  const size_t threads_per_team = 256;
  setup.setup(threads_per_team);

  // Run tests
  int ret = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Clean up
  setup.cleanup();

  // Finish up
  Kokkos::Cuda::host_mirror_device_type::finalize();
  Kokkos::Cuda::finalize();

  return ret;
}
