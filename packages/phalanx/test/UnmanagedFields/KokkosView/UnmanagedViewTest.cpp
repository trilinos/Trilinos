// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <limits>

#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_Field_UnmanagedAllocator.hpp"
#include "Phalanx_DataLayout_DynamicLayout.hpp"
#include "Phalanx_FieldTag_Tag.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Phalanx_Print.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "MyTraits.hpp"

#include "View_TestEvaluators.hpp"

TEUCHOS_UNIT_TEST(unmanaged_fields, basic)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

  // Put the field manager outside all evalaution code so that we can
  // test that memory management is handled correctly.
  FieldManager<MyTraits> fm;

  RCP<DataLayout> dl = rcp(new Layout("H-Grad<CELL,BASIS>",100,4));

  Tag<double> tag_a("a",dl);
  Tag<double> tag_b("b",dl);
  Tag<double> tag_c("c",dl);
  Tag<double> tag_d("d",dl);

  const auto tol = std::numeric_limits<double>::epsilon() * 100.0;

  {
    // Register evaluators
    {
      auto plist = Teuchos::parameterList("EvalUnmanaged");
      plist->set("dl",dl);
      using EU = EvalUnmanaged<MyTraits::Residual,MyTraits>;
      RCP<EU> e = rcp(new EU(*plist));
      e->setName("EvalUnmanaged");
      fm.registerEvaluator<MyTraits::Residual>(e);
    }
    {
      using ED = EvalDummy<MyTraits::Residual,MyTraits>;
      auto plist = Teuchos::parameterList("EvalDummy");
      plist->set("dl",dl);
      RCP<ED> e = rcp(new ED(*plist));
      e->setName("EvalDummy");
      fm.registerEvaluator<MyTraits::Residual>(e);
    }

    // Require fields
    {
      fm.requireField<MyTraits::Residual>(tag_a);

      fm.requireField<MyTraits::Residual>(tag_c);
    }

    Kokkos::View<double**,typename PHX::DevLayout<double>::type,PHX::Device> unmanaged_a("a",dl->extent(0),dl->extent(1));
    Kokkos::View<double**,typename PHX::DevLayout<double>::type,PHX::Device> unmanaged_b("b",dl->extent(0),dl->extent(1));
    Kokkos::View<double**,typename PHX::DevLayout<double>::type,PHX::Device> unmanaged_c("c",dl->extent(0),dl->extent(1));
    Kokkos::View<double**,typename PHX::DevLayout<double>::type,PHX::Device> unmanaged_d("d",dl->extent(0),dl->extent(1));
    Kokkos::deep_copy(unmanaged_a,0.0);
    Kokkos::deep_copy(unmanaged_b,5.0);
    Kokkos::deep_copy(unmanaged_c,0.0);
    Kokkos::deep_copy(unmanaged_d,6.0);

    auto host_a = Kokkos::create_mirror_view(unmanaged_a);
    auto host_b = Kokkos::create_mirror_view(unmanaged_b);
    auto host_c = Kokkos::create_mirror_view(unmanaged_c);
    auto host_d = Kokkos::create_mirror_view(unmanaged_d);
    Kokkos::deep_copy(host_a,unmanaged_a);
    Kokkos::deep_copy(host_b,unmanaged_b);
    Kokkos::deep_copy(host_c,unmanaged_c);
    Kokkos::deep_copy(host_d,unmanaged_d);

    for (int cell = 0; cell < unmanaged_a.extent_int(0); ++cell) {
      for (int basis = 0; basis < unmanaged_a.extent_int(1); ++basis) {
        TEST_FLOATING_EQUALITY(host_a(cell,basis),0.0,tol);
        TEST_FLOATING_EQUALITY(host_b(cell,basis),5.0,tol);
        TEST_FLOATING_EQUALITY(host_c(cell,basis),0.0,tol);
        TEST_FLOATING_EQUALITY(host_d(cell,basis),6.0,tol);
      }
    }

    // Register some unmanaged fields before postRegistrationSetup()
    // is called and some unmanaged fields after. There are two
    // branches through the code base depending on whether
    // postRegistrationSetup() has been called and we want to cover
    // both branches in testing.
    fm.setUnmanagedField<MyTraits::Residual>(tag_a,unmanaged_a);
    fm.setUnmanagedField<MyTraits::Residual>(tag_b,unmanaged_b);
    fm.postRegistrationSetup(0);
    fm.setUnmanagedField<MyTraits::Residual>(tag_c,unmanaged_c);
    fm.setUnmanagedField<MyTraits::Residual>(tag_d,unmanaged_d);
    fm.evaluateFields<MyTraits::Residual>(0);

    Kokkos::deep_copy(host_a,unmanaged_a);
    Kokkos::deep_copy(host_b,unmanaged_b);
    Kokkos::deep_copy(host_c,unmanaged_c);
    Kokkos::deep_copy(host_d,unmanaged_d);

    for (int cell = 0; cell < unmanaged_a.extent_int(0); ++cell) {
      for (int basis = 0; basis < unmanaged_a.extent_int(1); ++basis) {
        TEST_FLOATING_EQUALITY(host_a(cell,basis),5.0,tol);
        TEST_FLOATING_EQUALITY(host_b(cell,basis),5.0,tol);
        TEST_FLOATING_EQUALITY(host_c(cell,basis),6.0,tol);
        TEST_FLOATING_EQUALITY(host_d(cell,basis),6.0,tol);
      }
    }
  }

  // Once outside unmanaged_view scope, grab the views again and make
  // sure the internally stored views in the field manager are
  // pointing to the same memory.
  PHX::View<double**> a;
  PHX::View<double**> b;
  PHX::View<const double**> c; // test const accessors
  PHX::View<const double**> d; // test const accessors

  fm.getFieldData<MyTraits::Residual>(tag_a,a);
  fm.getFieldData<MyTraits::Residual>(tag_b,b);
  fm.getFieldData<MyTraits::Residual>(tag_c,c);
  fm.getFieldData<MyTraits::Residual>(tag_d,d);

  auto host_a = Kokkos::create_mirror_view(a);
  auto host_b = Kokkos::create_mirror_view(b);
  auto host_c = Kokkos::create_mirror_view(c);
  auto host_d = Kokkos::create_mirror_view(d);
  Kokkos::deep_copy(host_a,a);
  Kokkos::deep_copy(host_b,b);
  Kokkos::deep_copy(host_c,c);
  Kokkos::deep_copy(host_d,d);

  for (int cell = 0; cell < a.extent_int(0); ++cell) {
    for (int basis = 0; basis < a.extent_int(1); ++basis) {
      TEST_FLOATING_EQUALITY(host_a(cell,basis),5.0,tol);
      TEST_FLOATING_EQUALITY(host_b(cell,basis),5.0,tol);
      TEST_FLOATING_EQUALITY(host_c(cell,basis),6.0,tol);
      TEST_FLOATING_EQUALITY(host_d(cell,basis),6.0,tol);
    }
  }

}
