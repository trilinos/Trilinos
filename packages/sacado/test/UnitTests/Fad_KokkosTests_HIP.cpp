// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Kokkos_Macros.hpp"

// Temporarily disable DFad testing on HIP. HIP does not support "new"
// on device so temporary allocations don't work.
#ifdef KOKKOS_ENABLE_HIP
#define SACADO_TEST_DFAD 0
#else
#define SACADO_TEST_DFAD 1
#endif

#include "Fad_KokkosTests.hpp"

// Instantiate tests for HIP device.  DFAD is disabled since HIP doesn't support UVM.
using Kokkos::HIP;
VIEW_FAD_TESTS_D( HIP )

// Tests special size alignment for SFad on HIP is correct
TEUCHOS_UNIT_TEST(Kokkos_View_Fad, SFadHipAligned)
{
  const int StaticDim = 64;
  const int Stride = 32;
  const int LocalDim = 2;
  typedef Sacado::Fad::SFad<double,StaticDim> FadType;
  typedef Kokkos::LayoutContiguous<Kokkos::LayoutLeft,Stride> Layout;
  typedef Kokkos::HIP Device;
  typedef Kokkos::View<FadType*,Layout,Device> ViewType;

  typedef typename ViewType::traits TraitsType;
  typedef Kokkos::Impl::ViewMapping< TraitsType , typename TraitsType::specialize > MappingType;
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

TEUCHOS_UNIT_TEST(Kokkos_View_Fad, SFadHipNotAligned)
{
  const int StaticDim = 50;
  const int Stride = 32;
  const int LocalDim = 0;
  typedef Sacado::Fad::SFad<double,StaticDim> FadType;
  typedef Kokkos::LayoutContiguous<Kokkos::LayoutLeft,Stride> Layout;
  typedef Kokkos::HIP Device;
  typedef Kokkos::View<FadType*,Layout,Device> ViewType;

  typedef typename ViewType::traits TraitsType;
  typedef Kokkos::Impl::ViewMapping< TraitsType , typename TraitsType::specialize > MappingType;
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

  // Initialize HIP
  Kokkos::InitializationSettings init_args;
  init_args.set_device_id(0);
  Kokkos::initialize( init_args );
  Kokkos::print_configuration(std::cout);

  int res = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finalize HIP
  Kokkos::finalize();

  return res;
}
