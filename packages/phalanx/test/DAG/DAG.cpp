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
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldTag_Tag.hpp"
#include "Phalanx_Evaluator_Manager.hpp"
#include "Phalanx_TypeStrings.hpp"
#include "Phalanx_DimTag.hpp"

// Evaluators
#include "evaluators/Evaluator_Constant.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TimeMonitor.hpp"

// From test/Utilities directory
#include "Traits.hpp"

PHX_DIM_TAG_DECLARATION(CELL)
PHX_DIM_TAG_IMPLEMENTATION(CELL)

PHX_DIM_TAG_DECLARATION(BASIS)
PHX_DIM_TAG_IMPLEMENTATION(BASIS)

#include "Evaluator_MockDAG.hpp"

void registerDagNodes(PHX::EvaluatorManager<PHX::MyTraits>& em,
		      bool addCircularDependency,
		      bool buildDuplicateEvaluator,
		      bool removeRequiredFieldEvaluatorA,
		      bool removeFieldEvaluatorC)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

  using Mock = PHX::MockDAG<PHX::MyTraits::Residual,MyTraits>;

  if (!removeRequiredFieldEvaluatorA) {
    RCP<Mock> a = rcp(new Mock);
    a->setName("Eval_A");
    a->evaluates("A");
    a->requires("B");
    a->requires("C");
    em.registerEvaluator(a);
  }

  {
    RCP<Mock> b = rcp(new Mock);
    b->setName("Eval_B");
    b->evaluates("B");
    b->evaluates("D");
    b->requires("E");
    em.registerEvaluator(b);
  }

  if (!removeFieldEvaluatorC) {
    RCP<Mock> c = rcp(new Mock);
    c->setName("Eval_C");
    c->evaluates("C");
    c->requires("E");
    em.registerEvaluator(c);
  }

  {
    RCP<Mock> e = rcp(new Mock);
    e->setName("Eval_E");
    e->evaluates("E");
    if (addCircularDependency)
      e->requires("D");
    em.registerEvaluator(e);
  }

  // Add a second evaluator that can also evaluate "C"
  if (buildDuplicateEvaluator) {
    RCP<Mock> c = rcp(new Mock);
    c->setName("DUPLICATE Eval_C");
    c->evaluates("C");
    c->requires("E");
    em.registerEvaluator(c);
  }
}

TEUCHOS_UNIT_TEST(dag, basic_dag)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  EvaluatorManager<MyTraits> em;

  registerDagNodes(em,false,false,false,false);

  RCP<PHX::MDALayout<CELL,BASIS>> dl = 
    rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",100,4));
  PHX::Tag<MyTraits::Residual::ScalarT> tag("A",dl);
  em.requireField(tag);

  TEST_ASSERT(!em.sortingCalled());
  em.sortAndOrderEvaluators();
  TEST_ASSERT(em.sortingCalled());
  //std::cout << em << std::endl;

  const auto& order_new = em.getEvaluatorInternalOrdering();

  TEST_EQUALITY(order_new[0],3);
#ifdef PHX_ENABLE_NEW_DFS_ALGORITHM
  TEST_EQUALITY(order_new[1],2);
  TEST_EQUALITY(order_new[2],1);
#else
  TEST_EQUALITY(order_new[1],1);
  TEST_EQUALITY(order_new[2],2);
#endif
  TEST_EQUALITY(order_new[3],0);

  em.sortAndOrderEvaluatorsNew();
  //std::cout << em << std::endl;

  TEST_EQUALITY(order_new[0],3);
  TEST_EQUALITY(order_new[1],2);
  TEST_EQUALITY(order_new[2],1);
  TEST_EQUALITY(order_new[3],0);

  //Teuchos::TimeMonitor::summarize();
}

// Catch cyclic dependencies (not a true DAG)
TEUCHOS_UNIT_TEST(dag, cyclic)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  EvaluatorManager<MyTraits> em("cyclic");
  em.setDefaultGraphvizFilenameForErrors("error_cyclic.dot");
  em.setWriteGraphvizFileOnError(true);

  registerDagNodes(em,true,false,false,false);

  RCP<PHX::MDALayout<CELL,BASIS>> dl = 
    rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",100,4));
  PHX::Tag<MyTraits::Residual::ScalarT> tag("A",dl);
  em.requireField(tag);

  TEST_THROW(em.sortAndOrderEvaluatorsNew(),PHX::circular_dag_exception);
}

// Catch multiple evaluators that evaluate the same field
TEUCHOS_UNIT_TEST(dag, duplicate_evaluators)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  EvaluatorManager<MyTraits> em("duplicate_evaluators");
  em.setDefaultGraphvizFilenameForErrors("error_duplicate_evaluators.dot");
  em.setWriteGraphvizFileOnError(true);

#ifndef PHX_ALLOW_MULTIPLE_EVALUATORS_FOR_SAME_FIELD
  TEST_THROW(registerDagNodes(em,false,true,false,false),
	     PHX::multiple_evaluator_for_field_exception);
#else
  registerDagNodes(em,false,true,false,false);
#endif
}

// Catch missing required field 
TEUCHOS_UNIT_TEST(dag, missing_req_field)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  EvaluatorManager<MyTraits> em("missing_req_field");
  em.setDefaultGraphvizFilenameForErrors("error_missing_req_field.dot");
  em.setWriteGraphvizFileOnError(true);

  registerDagNodes(em,false,false,true,false);

  RCP<PHX::MDALayout<CELL,BASIS>> dl = 
    rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",100,4));
  PHX::Tag<MyTraits::Residual::ScalarT> tag("A",dl);
  em.requireField(tag);

  TEST_THROW(em.sortAndOrderEvaluatorsNew(),PHX::missing_evaluator_exception);
}

// Catch missing evalautor in subtree 
TEUCHOS_UNIT_TEST(dag, missing_evaluator)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  EvaluatorManager<MyTraits> em("missing_evaluator");
  em.setDefaultGraphvizFilenameForErrors("error_missing_evaluator.dot");
  em.setWriteGraphvizFileOnError(true);

  registerDagNodes(em,false,false,false,true);

  RCP<PHX::MDALayout<CELL,BASIS>> dl = 
    rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",100,4));
  PHX::Tag<MyTraits::Residual::ScalarT> tag("A",dl);
  em.requireField(tag);

  TEST_THROW(em.sortAndOrderEvaluatorsNew(),PHX::missing_evaluator_exception);
}

// Rebuild at runtime
TEUCHOS_UNIT_TEST(dag, rebuild)
{
  // Need to remove checks
}
