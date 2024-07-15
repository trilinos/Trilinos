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

#include "Field_TestEvaluators.hpp"

TEUCHOS_UNIT_TEST(unmanaged_fields, basic)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

  // Put the field manager outside all evalaution code so that we can
  // test that memory management is handled correctly.
  FieldManager<MyTraits> fm;

  RCP<DataLayout> dl = rcp(new Layout("H-Grad<CELL,BASIS>",100,4));

  const auto tol = std::numeric_limits<double>::epsilon() * 100.0;

  {
    // Register evaluators
    {
      auto plist = Teuchos::parameterList("EvalUnmanaged");
      using EU = EvalUnmanaged<MyTraits::Residual,MyTraits>;
      RCP<EU> e = rcp(new EU(*plist));
      e->setName("EvalUnmanaged");
      fm.registerEvaluator<MyTraits::Residual>(e);
    }
    {
      using ED = EvalDummy<MyTraits::Residual,MyTraits>;
      auto plist = Teuchos::parameterList("EvalDummy");
      RCP<ED> e = rcp(new ED(*plist));
      e->setName("EvalDummy");
      fm.registerEvaluator<MyTraits::Residual>(e);
    }

    // Require fields
    {
      PHX::Tag<MyTraits::Residual::ScalarT> tag_a("a",dl);
      fm.requireField<MyTraits::Residual>(tag_a);

      PHX::Tag<MyTraits::Residual::ScalarT> tag_c("c",dl);
      fm.requireField<MyTraits::Residual>(tag_c);
    }


    Field<double,2> unmanaged_a = allocateUnmanagedField<double,2>("a",dl);
    Field<double,2> unmanaged_b = allocateUnmanagedField<double,2>("b",dl);
    Field<double,2> unmanaged_c = allocateUnmanagedField<double,2>("c",dl);
    Field<double,2> unmanaged_d = allocateUnmanagedField<double,2>("d",dl);
    unmanaged_a.deep_copy(0.0);
    unmanaged_b.deep_copy(5.0);
    unmanaged_c.deep_copy(0.0);
    unmanaged_d.deep_copy(6.0);

    auto host_a = Kokkos::create_mirror_view(unmanaged_a.get_static_view());
    auto host_b = Kokkos::create_mirror_view(unmanaged_b.get_static_view());
    auto host_c = Kokkos::create_mirror_view(unmanaged_c.get_static_view());
    auto host_d = Kokkos::create_mirror_view(unmanaged_d.get_static_view());
    Kokkos::deep_copy(host_a,unmanaged_a.get_static_view());
    Kokkos::deep_copy(host_b,unmanaged_b.get_static_view());
    Kokkos::deep_copy(host_c,unmanaged_c.get_static_view());
    Kokkos::deep_copy(host_d,unmanaged_d.get_static_view());

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
    fm.setUnmanagedField<MyTraits::Residual>(unmanaged_a);
    fm.setUnmanagedField<MyTraits::Residual>(unmanaged_b);
    fm.postRegistrationSetup(0);
    fm.setUnmanagedField<MyTraits::Residual>(unmanaged_c);
    fm.setUnmanagedField<MyTraits::Residual>(unmanaged_d);
    fm.evaluateFields<MyTraits::Residual>(0);

    Kokkos::deep_copy(host_a,unmanaged_a.get_static_view());
    Kokkos::deep_copy(host_b,unmanaged_b.get_static_view());
    Kokkos::deep_copy(host_c,unmanaged_c.get_static_view());
    Kokkos::deep_copy(host_d,unmanaged_d.get_static_view());

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
  Field<double,2> a("a",dl);
  Field<double,2> b("b",dl);
  Field<const double,2> c("c",dl); // test const version
  Field<const double,2> d("d",dl); // test const version

  fm.getFieldData<MyTraits::Residual>(a);
  fm.getFieldData<MyTraits::Residual>(b);
  fm.getFieldData<MyTraits::Residual>(c);
  fm.getFieldData<MyTraits::Residual>(d);

  auto host_a = Kokkos::create_mirror_view(a.get_static_view());
  auto host_b = Kokkos::create_mirror_view(b.get_static_view());
  auto host_c = Kokkos::create_mirror_view(c.get_static_view());
  auto host_d = Kokkos::create_mirror_view(d.get_static_view());
  Kokkos::deep_copy(host_a,a.get_static_view());
  Kokkos::deep_copy(host_b,b.get_static_view());
  Kokkos::deep_copy(host_c,c.get_static_view());
  Kokkos::deep_copy(host_d,d.get_static_view());

  for (int cell = 0; cell < a.extent_int(0); ++cell) {
    for (int basis = 0; basis < a.extent_int(1); ++basis) {
      TEST_FLOATING_EQUALITY(host_a(cell,basis),5.0,tol);
      TEST_FLOATING_EQUALITY(host_b(cell,basis),5.0,tol);
      TEST_FLOATING_EQUALITY(host_c(cell,basis),6.0,tol);
      TEST_FLOATING_EQUALITY(host_d(cell,basis),6.0,tol);
    }
  }

}
