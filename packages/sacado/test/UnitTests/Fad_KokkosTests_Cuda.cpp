// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Fad_KokkosTests.hpp"

#include "Kokkos_Core.hpp"

// Instantiate tests for Cuda device.  We can only test DFad is UVM is enabled.
#if defined(KOKKOS_ENABLE_CUDA_UVM)
#define SACADO_TEST_DFAD 1
#else
#define SACADO_TEST_DFAD 0
#endif
using Kokkos::Cuda;
VIEW_FAD_TESTS_D( Cuda )

// Tests special size alignment for SFad on Cuda is correct
TEUCHOS_UNIT_TEST(Kokkos_View_Fad, SFadCudaAligned)
{
  const int StaticDim = 64;
  const int Stride = 32;
  const int LocalDim = 2;
  typedef Sacado::Fad::SFad<double,StaticDim> FadType;
  typedef Kokkos::LayoutContiguous<Kokkos::LayoutLeft,Stride> Layout;
  typedef Kokkos::Cuda Device;
  typedef Kokkos::View<FadType*,Layout,Device> ViewType;

  typedef typename ViewType::traits TraitsType;
  typedef Kokkos::Impl::ViewMapping< TraitsType , void > MappingType;
  const int view_static_dim = MappingType::FadStaticDimension;
  TEUCHOS_TEST_EQUALITY(view_static_dim, StaticDim, out, success);

  typedef typename Kokkos::ThreadLocalScalarType<ViewType>::type local_fad_type;
  const bool issfd = is_sfad<local_fad_type>::value;
  const int static_dim = Sacado::StaticSize<local_fad_type>::value;
  TEUCHOS_TEST_EQUALITY(issfd, true, out, success);
  TEUCHOS_TEST_EQUALITY(static_dim, LocalDim, out, success);

  const size_t num_rows = 11;
  const size_t fad_size = StaticDim;

  ViewType v("v", num_rows, fad_size+1);
  const size_t span = v.span();
  TEUCHOS_TEST_EQUALITY(span, num_rows*(StaticDim+1), out, success);
}

TEUCHOS_UNIT_TEST(Kokkos_View_Fad, SFadCudaNotAligned)
{
  const int StaticDim = 50;
  const int Stride = 32;
  const int LocalDim = 0;
  typedef Sacado::Fad::SFad<double,StaticDim> FadType;
  typedef Kokkos::LayoutContiguous<Kokkos::LayoutLeft,Stride> Layout;
  typedef Kokkos::Cuda Device;
  typedef Kokkos::View<FadType*,Layout,Device> ViewType;

  typedef typename ViewType::traits TraitsType;
  typedef Kokkos::Impl::ViewMapping< TraitsType , void > MappingType;
  const int view_static_dim = MappingType::FadStaticDimension;
  TEUCHOS_TEST_EQUALITY(view_static_dim, StaticDim, out, success);

  typedef typename Kokkos::ThreadLocalScalarType<ViewType>::type local_fad_type;
  const bool issfd = is_sfad<local_fad_type>::value;
  const int static_dim = Sacado::StaticSize<local_fad_type>::value;
  TEUCHOS_TEST_EQUALITY(issfd, false, out, success);
  TEUCHOS_TEST_EQUALITY(static_dim, LocalDim, out, success);

  const size_t num_rows = 11;
  const size_t fad_size = StaticDim;

  ViewType v("v", num_rows, fad_size+1);
  const size_t span = v.span();
  TEUCHOS_TEST_EQUALITY(span, num_rows*(StaticDim+1), out, success);
}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Initialize Cuda
  Kokkos::InitArguments init_args;
  init_args.device_id = 0;
  Kokkos::initialize( init_args );
  Kokkos::print_configuration(std::cout);

  int res = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finalize Cuda
  Kokkos::finalize();

  return res;
}
