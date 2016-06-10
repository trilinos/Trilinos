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
#include "Phalanx_DAG_Manager.hpp"
#include "Phalanx_TypeStrings.hpp"
#include "Phalanx_DimTag.hpp"

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

void registerDagNodes(PHX::DagManager<PHX::MyTraits>& em,
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

  DagManager<MyTraits> em;

  registerDagNodes(em,false,false,false,false);

  RCP<PHX::MDALayout<CELL,BASIS>> dl = 
    rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",100,4));
  PHX::Tag<MyTraits::Residual::ScalarT> tag("A",dl);
  em.requireField(tag);

  TEST_ASSERT(!em.sortingCalled());
  em.sortAndOrderEvaluators();
  TEST_ASSERT(em.sortingCalled());
  //std::cout << em << std::endl;

  {
    const auto& order_new = em.getEvaluatorInternalOrdering();
    TEST_EQUALITY(order_new[0],3);
    TEST_EQUALITY(order_new[1],2);
    TEST_EQUALITY(order_new[2],1);
    TEST_EQUALITY(order_new[3],0);
  }

  // Make sure we can call again
  em.sortAndOrderEvaluators();
  {
    const auto& order_new = em.getEvaluatorInternalOrdering();
    //std::cout << em << std::endl;
    TEST_EQUALITY(order_new[0],3);
    TEST_EQUALITY(order_new[1],2);
    TEST_EQUALITY(order_new[2],1);
    TEST_EQUALITY(order_new[3],0);
  }

  // Test some class methods
  const std::vector< Teuchos::RCP<PHX::FieldTag> >& tags = em.getFieldTags();
  TEST_EQUALITY(tags.size(),5);
  em.writeGraphvizFile("basic_dag.dot",true,true,false);
  cout << "\n" << em << endl;   cout << "\n" << em << endl;  
 
  //Teuchos::TimeMonitor::summarize();
  PHX::FinalizeKokkosDevice();
}

// Catch cyclic dependencies (not a true DAG)
TEUCHOS_UNIT_TEST(dag, cyclic)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  PHX::InitializeKokkosDevice();

  DagManager<MyTraits> em("cyclic");
  em.setDefaultGraphvizFilenameForErrors("error_cyclic.dot");
  em.setWriteGraphvizFileOnError(true);

  registerDagNodes(em,true,false,false,false);

  RCP<PHX::MDALayout<CELL,BASIS>> dl = 
    rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",100,4));
  PHX::Tag<MyTraits::Residual::ScalarT> tag("A",dl);
  em.requireField(tag);

  TEST_THROW(em.sortAndOrderEvaluators(),PHX::circular_dag_exception);
}

// Catch multiple evaluators that evaluate the same field
TEUCHOS_UNIT_TEST(dag, duplicate_evaluators)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

  DagManager<MyTraits> em("duplicate_evaluators");
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

  DagManager<MyTraits> em("missing_req_field");
  em.setDefaultGraphvizFilenameForErrors("error_missing_req_field.dot");
  em.setWriteGraphvizFileOnError(true);

  registerDagNodes(em,false,false,true,false);

  RCP<PHX::MDALayout<CELL,BASIS>> dl = 
    rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",100,4));
  PHX::Tag<MyTraits::Residual::ScalarT> tag("A",dl);
  em.requireField(tag);

  TEST_THROW(em.sortAndOrderEvaluators(),PHX::missing_evaluator_exception);
}

// Catch missing evalautor in subtree 
TEUCHOS_UNIT_TEST(dag, missing_evaluator)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
    
  DagManager<MyTraits> em("missing_evaluator");
  em.setDefaultGraphvizFilenameForErrors("error_missing_evaluator.dot");
  em.setWriteGraphvizFileOnError(true);

  registerDagNodes(em,false,false,false,true);

  RCP<PHX::MDALayout<CELL,BASIS>> dl = 
    rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",100,4));
  PHX::Tag<MyTraits::Residual::ScalarT> tag("A",dl);
  em.requireField(tag);

  TEST_THROW(em.sortAndOrderEvaluators(),PHX::missing_evaluator_exception);
}

// Test the analyzeGraph computation for speedup and parallelization
TEUCHOS_UNIT_TEST(dag, analyze_graph)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  using Mock = PHX::MockDAG<PHX::MyTraits::Residual,MyTraits>;

  // Perfectly parallel test
  DagManager<MyTraits> dag("analyze_graph");

  // Register evaluators
  {
    RCP<Mock> m = rcp(new Mock);
    m->setName("Eval_A");
    m->evaluates("A");
    dag.registerEvaluator(m);
  }
  {
    RCP<Mock> m = rcp(new Mock);
    m->setName("Eval_B");
    m->evaluates("B");
    dag.registerEvaluator(m);
  }
  {
    RCP<Mock> m = rcp(new Mock);
    m->setName("Eval_C");
    m->evaluates("C");
    dag.registerEvaluator(m);
  }

  // Require fields
  {
    RCP<MDALayout<CELL,BASIS>> dl = 
      rcp(new MDALayout<CELL,BASIS>("H-Grad",100,4));
    Tag<MyTraits::Residual::ScalarT> taga("A",dl);
    dag.requireField(taga);
    Tag<MyTraits::Residual::ScalarT> tagb("B",dl);
    dag.requireField(tagb);
    Tag<MyTraits::Residual::ScalarT> tagc("C",dl);
    dag.requireField(tagc);
  }

  dag.sortAndOrderEvaluators();

  const std::vector<int>& order = dag.getEvaluatorInternalOrdering();
  const std::vector<PHX::DagNode<MyTraits>>& nodes = dag.getDagNodes();

  // Insert evaluation times into the graph
  for (const auto& i : order) {
    const DagNode<MyTraits>& n = nodes[i];
    // cast away const for unit testing to set the execution times
    std::chrono::duration<double> dt(1.0);
    const_cast<DagNode<MyTraits>&>(n).setExecutionTime(dt);
  }

  double speedup = 0.0;
  double parallelizability = 0.0;
  dag.analyzeGraph(speedup, parallelizability);

  out <<  "Speedup = " << speedup << std::endl;
  out <<  "Parallelizability = " << parallelizability << std::endl;

  double tol = 1000.0 * Teuchos::ScalarTraits<double>::eps();
  TEST_FLOATING_EQUALITY(speedup,3.0,tol);
  TEST_FLOATING_EQUALITY(parallelizability,1.0,tol);
}

// Test the analyzeGraph computation for speedup and parallelization
TEUCHOS_UNIT_TEST(dag, analyze_graph2)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  using Mock = PHX::MockDAG<PHX::MyTraits::Residual,MyTraits>;

  // Perfectly parallel test
  DagManager<MyTraits> dag("analyze_graph2");

  // Register evaluators
  {
    RCP<Mock> m = rcp(new Mock);
    m->setName("Eval_A");
    m->evaluates("A");
    m->requires("B");
    m->requires("C");
    dag.registerEvaluator(m);
  }
  {
    RCP<Mock> m = rcp(new Mock);
    m->setName("Eval_B");
    m->evaluates("B");
    m->requires("D");
    dag.registerEvaluator(m);
  }
  {
    RCP<Mock> m = rcp(new Mock);
    m->setName("Eval_C");
    m->evaluates("C");
    m->requires("D");
    dag.registerEvaluator(m);
  }
  {
    RCP<Mock> m = rcp(new Mock);
    m->setName("Eval_D");
    m->evaluates("D");
    dag.registerEvaluator(m);
  }

  // Require fields
  {
    RCP<MDALayout<CELL,BASIS>> dl = 
      rcp(new MDALayout<CELL,BASIS>("H-Grad",100,4));
    Tag<MyTraits::Residual::ScalarT> taga("A",dl);
    dag.requireField(taga);
  }

  dag.sortAndOrderEvaluators();

  const std::vector<int>& order = dag.getEvaluatorInternalOrdering();
  const std::vector<PHX::DagNode<MyTraits>>& nodes = dag.getDagNodes();

  // Insert evaluation times into the graph
  for (const auto& i : order) {
    const DagNode<MyTraits>& n = nodes[i];

    // cast away const for unit testing to set the execution times
    if (n.get()->getName() == "Eval_A") {
      std::chrono::duration<double> dt(2.0);
      const_cast<DagNode<MyTraits>&>(n).setExecutionTime(dt);
    }
    else if (n.get()->getName() == "Eval_B") {
      std::chrono::duration<double> dt(2.0);
      const_cast<DagNode<MyTraits>&>(n).setExecutionTime(dt);
    }
    else if (n.get()->getName() == "Eval_C") {
      std::chrono::duration<double> dt(4.0);
      const_cast<DagNode<MyTraits>&>(n).setExecutionTime(dt);
    }
    else if (n.get()->getName() == "Eval_D") {
      std::chrono::duration<double> dt(2.0);
      const_cast<DagNode<MyTraits>&>(n).setExecutionTime(dt);
    }
  }

  double speedup = 0.0;
  double parallelizability = 0.0;
  dag.analyzeGraph(speedup, parallelizability);

  out <<  "Speedup = " << speedup << std::endl;
  out <<  "Parallelizability = " << parallelizability << std::endl;

  double tol = 1000.0 * Teuchos::ScalarTraits<double>::eps();
  double s_gold = 10.0 / 8.0;
  TEST_FLOATING_EQUALITY(speedup,s_gold,tol);
  double p_gold = (1. - 1./s_gold)/(1. - 1./4.);
  TEST_FLOATING_EQUALITY(parallelizability,p_gold,tol);

  // Now test the sum into execution time to accumulate timings
  // Insert evaluation times into the graph
  for (const auto& i : order) {
    const DagNode<MyTraits>& n = nodes[i];

    // cast away const for unit testing to set the execution times
    if (n.get()->getName() == "Eval_A") {
      std::chrono::duration<double> dt(2.0);
      const_cast<DagNode<MyTraits>&>(n).sumIntoExecutionTime(dt);
    }
    else if (n.get()->getName() == "Eval_B") {
      // force different ratio to test sumInto function
      std::chrono::duration<double> dt(8.0);
      const_cast<DagNode<MyTraits>&>(n).sumIntoExecutionTime(dt);
    }
    else if (n.get()->getName() == "Eval_C") {
      std::chrono::duration<double> dt(4.0);
      const_cast<DagNode<MyTraits>&>(n).sumIntoExecutionTime(dt);
    }
    else if (n.get()->getName() == "Eval_D") {
      std::chrono::duration<double> dt(2.0);
      const_cast<DagNode<MyTraits>&>(n).sumIntoExecutionTime(dt);
    }
  }
  
  dag.analyzeGraph(speedup, parallelizability);

  s_gold = 26.0 / 18.0;
  TEST_FLOATING_EQUALITY(speedup,s_gold,tol);
  p_gold = (1. - 1./s_gold)/(1. - 1./4.);
  TEST_FLOATING_EQUALITY(parallelizability,p_gold,tol);
}
