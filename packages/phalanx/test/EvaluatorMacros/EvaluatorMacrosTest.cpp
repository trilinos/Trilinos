// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Phalanx_KokkosDeviceTypes.hpp"
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

#include "EvaluatorWithMacros.hpp"

TEUCHOS_UNIT_TEST(evaluator_macros, basic)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

  FieldManager<MyTraits> fm;

  // Mock evaluators for testing MACRO definitions
  using Ev1 = EvaluatorWithMacros1<MyTraits::Residual,MyTraits>;
  using Ev2 = EvaluatorWithMacros2<MyTraits::Residual,MyTraits>;
  
  {
    auto plist_a = Teuchos::parameterList("A"); 
    RCP<Ev1> a = rcp(new Ev1(*plist_a));
    a->setName("Eval_A");
    a->evaluates("A");
    a->depends("B");
    a->depends("C");
    fm.registerEvaluator<MyTraits::Residual>(a);
  }
  {
    auto plist_b = Teuchos::parameterList("B"); 
    RCP<Ev2> b = rcp(new Ev2(*plist_b));
    b->setName("Eval_B");
    b->evaluates("B");
    b->depends("D");
    fm.registerEvaluator<MyTraits::Residual>(b);
  }
  {
    auto plist_c = Teuchos::parameterList("C"); 
    RCP<Ev2> c = rcp(new Ev2(*plist_c));
    c->setName("Eval_C");
    c->evaluates("C");
    c->depends("D");
    fm.registerEvaluator<MyTraits::Residual>(c);
  }
  {
    auto plist_d = Teuchos::parameterList("D"); 
    RCP<Ev1> d = rcp(new Ev1(*plist_d));
    d->setName("Eval_D");
    d->evaluates("D");
    fm.registerEvaluator<MyTraits::Residual>(d);
  }

  RCP<PHX::MDALayout<CELL,BASIS>> dl = 
    rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",100,4));
  PHX::Tag<MyTraits::Residual::ScalarT> tag("A",dl);
  fm.requireField<MyTraits::Residual>(tag);

  fm.postRegistrationSetup(0);
  fm.preEvaluate<MyTraits::Residual>(1);
  fm.evaluateFields<MyTraits::Residual>(1);
  fm.postEvaluate<MyTraits::Residual>(1);
}
