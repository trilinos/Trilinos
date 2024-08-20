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
#include "Phalanx_MDField_UnmanagedAllocator.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
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

PHX_EXTENT(CELL)
PHX_EXTENT(BASIS)

#include "AliasedFields_TestEvaluators.hpp"

// ****************************************
// Check that alias field works correctly
// ****************************************
TEUCHOS_UNIT_TEST(aliased_fields, basic)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

  // Put the field manager outside all evalaution code so that we can
  // test that memory management is handled correctly.
  FieldManager<MyTraits> fm;

  RCP<MDALayout<CELL,BASIS>> dl = rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",100,4));

  // Register evaluators
  {
    RCP<EvalA<MyTraits::Residual,MyTraits>> a = rcp(new EvalA<MyTraits::Residual,MyTraits>(dl));
    fm.registerEvaluator<MyTraits::Residual>(a);
  }
  {
    RCP<EvalC<MyTraits::Residual,MyTraits>> c = rcp(new EvalC<MyTraits::Residual,MyTraits>(dl));
    fm.registerEvaluator<MyTraits::Residual>(c);
  }

  // Require field a
  {
    PHX::Tag<MyTraits::Residual::ScalarT> tag_a("a",dl);
    fm.requireField<MyTraits::Residual>(tag_a);
  }

  // Alias "b" to "c"
  {
    PHX::Tag<MyTraits::Residual::ScalarT> tag_b("b",dl);
    PHX::Tag<MyTraits::Residual::ScalarT> tag_c("c",dl);
    fm.aliasFieldForAllEvaluationTypes(tag_b,tag_c);
  }

  fm.postRegistrationSetup(0);
  fm.evaluateFields<MyTraits::Residual>(0);
  Kokkos::fence();

  MDField<double,CELL,BASIS> b("b",dl);
  MDField<double,CELL,BASIS> c("c",dl);
  fm.getFieldData<MyTraits::Residual,double,CELL,BASIS>(b);
  fm.getFieldData<MyTraits::Residual,double,CELL,BASIS>(c);

  auto host_b = Kokkos::create_mirror_view(b.get_static_view());
  auto host_c = Kokkos::create_mirror_view(c.get_static_view());
  Kokkos::deep_copy(host_b,b.get_static_view());
  Kokkos::deep_copy(host_c,c.get_static_view());

  const auto tol = std::numeric_limits<double>::epsilon() * 100.0;

  for (int cell = 0; cell < b.extent_int(0); ++cell) {
    for (int basis = 0; basis < c.extent_int(1); ++basis) {
      TEST_FLOATING_EQUALITY(host_b(cell,basis),4.0,tol);
      TEST_FLOATING_EQUALITY(host_c(cell,basis),4.0,tol);
    }
  }
}

// ****************************************
// Check that we catch bad data layouts in the aliased field
// ****************************************
TEUCHOS_UNIT_TEST(aliased_fields, throw_on_bad_layout)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

  // Put the field manager outside all evalaution code so that we can
  // test that memory management is handled correctly.
  FieldManager<MyTraits> fm;

  RCP<MDALayout<CELL,BASIS>> dl = rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",100,4));
  RCP<MDALayout<CELL,BASIS>> bad_dl = rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",50,4));

  // Register evaluators
  {
    RCP<EvalA<MyTraits::Residual,MyTraits>> a = rcp(new EvalA<MyTraits::Residual,MyTraits>(dl));
    fm.registerEvaluator<MyTraits::Residual>(a);
  }
  {
    RCP<EvalC<MyTraits::Residual,MyTraits>> c = rcp(new EvalC<MyTraits::Residual,MyTraits>(dl));
    fm.registerEvaluator<MyTraits::Residual>(c);
  }

  // Require field a
  {
    PHX::Tag<MyTraits::Residual::ScalarT> tag_a("a",dl);
    fm.requireField<MyTraits::Residual>(tag_a);
  }

  // Alias "b" to "c"
  {
    // use an incorrectly sized data layout
    PHX::Tag<MyTraits::Residual::ScalarT> tag_b("b",bad_dl);
    PHX::Tag<MyTraits::Residual::ScalarT> tag_c("c",dl);
    TEST_THROW(fm.aliasFieldForAllEvaluationTypes(tag_b,tag_c),std::runtime_error)
  }

}

// ****************************************
// Check that we catch bad scalar types in the aliased field
// ****************************************
TEUCHOS_UNIT_TEST(aliased_fields, throw_on_bad_scalar_type)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

  // Put the field manager outside all evalaution code so that we can
  // test that memory management is handled correctly.
  FieldManager<MyTraits> fm;

  RCP<MDALayout<CELL,BASIS>> dl = rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",100,4));

  // Register evaluators
  {
    RCP<EvalA<MyTraits::Residual,MyTraits>> a = rcp(new EvalA<MyTraits::Residual,MyTraits>(dl));
    fm.registerEvaluator<MyTraits::Residual>(a);
  }
  {
    RCP<EvalC<MyTraits::Residual,MyTraits>> c = rcp(new EvalC<MyTraits::Residual,MyTraits>(dl));
    fm.registerEvaluator<MyTraits::Residual>(c);
  }

  // Require field a
  {
    PHX::Tag<MyTraits::Residual::ScalarT> tag_a("a",dl);
    fm.requireField<MyTraits::Residual>(tag_a);
  }

  // Alias "b" to "c"
  {
    // inject a float instead of double
    PHX::Tag<float> tag_b("b",dl);
    PHX::Tag<MyTraits::Residual::ScalarT> tag_c("c",dl);
    TEST_THROW(fm.aliasFieldForAllEvaluationTypes(tag_b,tag_c),std::runtime_error)
  }

}
