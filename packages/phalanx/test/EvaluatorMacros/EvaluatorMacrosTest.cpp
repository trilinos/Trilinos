// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#include "Phalanx_config.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_KokkosUtilities.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldTag_Tag.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Phalanx_TypeStrings.hpp"
#include "Phalanx_DimTag.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "MyTraits.hpp"

PHX_DIM_TAG_DECLARATION(CELL)
PHX_DIM_TAG_IMPLEMENTATION(CELL)

PHX_DIM_TAG_DECLARATION(BASIS)
PHX_DIM_TAG_IMPLEMENTATION(BASIS)

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
    a->requires("B");
    a->requires("C");
    fm.registerEvaluator<MyTraits::Residual>(a);
  }
  {
    auto plist_b = Teuchos::parameterList("B"); 
    RCP<Ev2> b = rcp(new Ev2(*plist_b));
    b->setName("Eval_B");
    b->evaluates("B");
    b->requires("D");
    fm.registerEvaluator<MyTraits::Residual>(b);
  }
  {
    auto plist_c = Teuchos::parameterList("C"); 
    RCP<Ev2> c = rcp(new Ev2(*plist_c));
    c->setName("Eval_C");
    c->evaluates("C");
    c->requires("D");
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
