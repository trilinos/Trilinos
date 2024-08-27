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
#include "Phalanx_DAG_Manager.hpp"
#include "Phalanx_Evaluator_AliasField.hpp"
#include "Phalanx_Print.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TimeMonitor.hpp"

// From test/Utilities directory
#include "Traits.hpp"

PHX_EXTENT(CELL)
PHX_EXTENT(BASIS)

#include "Evaluator_MockDAG.hpp"

// *************************************************
// functions for testing
// *************************************************
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
    a->depends("B");
    a->depends("C");
    em.registerEvaluator(a);
  }

  {
    RCP<Mock> b = rcp(new Mock);
    b->setName("Eval_B");
    b->evaluates("B");
    b->evaluates("D");
    b->depends("E");
    em.registerEvaluator(b);
  }

  if (!removeFieldEvaluatorC) {
    RCP<Mock> c = rcp(new Mock);
    c->setName("Eval_C");
    c->evaluates("C");
    c->depends("E");
    em.registerEvaluator(c);
  }

  {
    RCP<Mock> e = rcp(new Mock);
    e->setName("Eval_E");
    e->evaluates("E");
    if (addCircularDependency)
      e->depends("D");
    em.registerEvaluator(e);
  }

  // Add a second evaluator that can also evaluate "C"
  if (buildDuplicateEvaluator) {
    RCP<Mock> c = rcp(new Mock);
    c->setName("DUPLICATE Eval_C");
    c->evaluates("C");
    c->depends("E");
    em.registerEvaluator(c);
  }
}

// *************************************************
// Successful DAG
// *************************************************
TEUCHOS_UNIT_TEST(dag, basic_dag)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

  DagManager<MyTraits> em;

  registerDagNodes(em,false,false,false,false);

  RCP<PHX::MDALayout<CELL,BASIS>> dl =
    rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",100,4));
  PHX::Tag<MyTraits::Residual::ScalarT> tag_a("A",dl);
  em.requireField(tag_a);

  TEST_ASSERT(!em.sortingCalled());
  em.sortAndOrderEvaluators();
  TEST_ASSERT(em.sortingCalled());

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
    TEST_EQUALITY(order_new[0],3);
    TEST_EQUALITY(order_new[1],2);
    TEST_EQUALITY(order_new[2],1);
    TEST_EQUALITY(order_new[3],0);
  }

  // Test some class methods
  const std::vector< Teuchos::RCP<PHX::FieldTag> >& tags = em.getFieldTags();
  TEST_EQUALITY(tags.size(),5);
  em.writeGraphvizFile("basic_dag.dot",true,true,false);
  std::stringstream output;
  output << em << endl;

  {
    auto& evaluators = em.getEvaluatorsBindingField(tag_a);
    TEST_EQUALITY(evaluators.size(),1);
  }

  {
    PHX::Tag<MyTraits::Residual::ScalarT> tag("B",dl);
    auto& evaluators = em.getEvaluatorsBindingField(tag);
    TEST_EQUALITY(evaluators.size(),2);
  }

  {
    PHX::Tag<MyTraits::Residual::ScalarT> tag("C",dl);
    auto& evaluators = em.getEvaluatorsBindingField(tag);
    TEST_EQUALITY(evaluators.size(),2);
  }

  {
    PHX::Tag<MyTraits::Residual::ScalarT> tag("E",dl);
    auto& evaluators = em.getEvaluatorsBindingField(tag);
    TEST_EQUALITY(evaluators.size(),3);
  }

}

// *************************************************
// Catch cyclic dependencies (not a true DAG)
// *************************************************
TEUCHOS_UNIT_TEST(dag, cyclic)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

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

// *************************************************
// Catch multiple evaluators that evaluate the same field
// *************************************************
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

// *************************************************
// Catch missing required field
// *************************************************
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

// *************************************************
// Catch missing evaluator in subtree
// *************************************************
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

// *************************************************
// Test the analyzeGraph computation for speedup and parallelization
// *************************************************
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

// *************************************************
// Test the analyzeGraph computation for speedup and parallelization
// *************************************************
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
    m->depends("B");
    m->depends("C");
    dag.registerEvaluator(m);
  }
  {
    RCP<Mock> m = rcp(new Mock);
    m->setName("Eval_B");
    m->evaluates("B");
    m->depends("D");
    dag.registerEvaluator(m);
  }
  {
    RCP<Mock> m = rcp(new Mock);
    m->setName("Eval_C");
    m->evaluates("C");
    m->depends("D");
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

// *************************************************
// Test for a field that has both an "evaluated" evaluator and
// "contributed" evaluators.
// *************************************************
TEUCHOS_UNIT_TEST(dag, contrib_and_eval_B)
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
    m->depends("B");
    m->depends("C");
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
    m->depends("D");
    dag.registerEvaluator(m);
  }
  {
    RCP<Mock> m = rcp(new Mock);
    m->setName("Eval_D");
    m->evaluates("D");
    dag.registerEvaluator(m);
  }
  { // Contributes to B
    RCP<Mock> m = rcp(new Mock);
    m->setName("Eval_B+");
    m->contributes("B");
    m->depends("D");
    dag.registerEvaluator(m);
  }
  { // Contributes to B also
    RCP<Mock> m = rcp(new Mock);
    m->setName("Eval_B++");
    m->contributes("B");
    m->depends("D");
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
  out << dag << std::endl;
  dag.writeGraphvizFile("contrib_and_eval_B.dot",true,true,true);

  // Check the ordering
  {
    const auto& order_new = dag.getEvaluatorInternalOrdering();
    //const std::vector<PHX::DagNode<MyTraits>>& nodes = dag.getDagNodes();
    TEST_EQUALITY(order_new[0],3);
    TEST_EQUALITY(order_new[1],2);
    TEST_EQUALITY(order_new[2],1);
    TEST_EQUALITY(order_new[3],4);
    TEST_EQUALITY(order_new[4],5);
    TEST_EQUALITY(order_new[5],0);
  }

  // Check that the out edges are correct.
  {
    const std::vector<PHX::DagNode<MyTraits>>& nodes = dag.getDagNodes();
    // A
    TEST_EQUALITY(nodes[0].adjacencies().size(),4);
    TEST_ASSERT(nodes[0].adjacencies().find(1) != nodes[0].adjacencies().end());
    TEST_ASSERT(nodes[0].adjacencies().find(2) != nodes[0].adjacencies().end());
    TEST_ASSERT(nodes[0].adjacencies().find(4) != nodes[0].adjacencies().end());
    TEST_ASSERT(nodes[0].adjacencies().find(5) != nodes[0].adjacencies().end());
    // B
    TEST_EQUALITY(nodes[1].adjacencies().size(),0);
    // B+
    TEST_EQUALITY(nodes[4].adjacencies().size(),2);
    TEST_ASSERT(nodes[4].adjacencies().find(1) != nodes[4].adjacencies().end());
    TEST_ASSERT(nodes[4].adjacencies().find(3) != nodes[4].adjacencies().end());
    // B++
    TEST_EQUALITY(nodes[5].adjacencies().size(),2);
    TEST_ASSERT(nodes[5].adjacencies().find(1) != nodes[5].adjacencies().end());
    TEST_ASSERT(nodes[5].adjacencies().find(3) != nodes[5].adjacencies().end());
  }
}

// *************************************************
// Test for a field that is only evaluated by contributed fields
// *************************************************
TEUCHOS_UNIT_TEST(dag, contrib_only_B)
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
    m->depends("B");
    m->depends("C");
    dag.registerEvaluator(m);
  }
  {
    RCP<Mock> m = rcp(new Mock);
    m->setName("Eval_C");
    m->evaluates("C");
    m->depends("D");
    dag.registerEvaluator(m);
  }
  {
    RCP<Mock> m = rcp(new Mock);
    m->setName("Eval_D");
    m->evaluates("D");
    dag.registerEvaluator(m);
  }
  { // Contributes to B
    RCP<Mock> m = rcp(new Mock);
    m->setName("Eval_B+");
    m->contributes("B");
    m->depends("D");
    dag.registerEvaluator(m);
  }
  { // Contributes to B also
    RCP<Mock> m = rcp(new Mock);
    m->setName("Eval_B++");
    m->contributes("B");
    m->depends("D");
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
  //out << dag << std::endl;
  dag.writeGraphvizFile("contrib_only_B.dot",true,true,true);

  // Check the ordering
  {
    const auto& order_new = dag.getEvaluatorInternalOrdering();
    //const std::vector<PHX::DagNode<MyTraits>>& nodes = dag.getDagNodes();
    TEST_EQUALITY(order_new[0],2);
    // Node 2 is first, nodes 1, 3, 4 can happen in any order, then node 0 last.
    TEST_ASSERT( (order_new[1] == 1) || (order_new[1] == 3) || (order_new[1] == 4) );
    TEST_ASSERT( (order_new[2] == 1) || (order_new[2] == 3) || (order_new[2] == 4) );
    TEST_ASSERT( (order_new[3] == 1) || (order_new[3] == 3) || (order_new[3] == 4) );
    TEST_INEQUALITY(order_new[1],order_new[3]);
    TEST_INEQUALITY(order_new[1],order_new[4]);
    TEST_INEQUALITY(order_new[3],order_new[4]);
    TEST_EQUALITY(order_new[4],0);
  }

  // Check that the out edges are correct.
  {
    const std::vector<PHX::DagNode<MyTraits>>& nodes = dag.getDagNodes();
    // A
    TEST_EQUALITY(nodes[0].adjacencies().size(),3);
    TEST_ASSERT(nodes[0].adjacencies().find(1) != nodes[0].adjacencies().end());
    TEST_ASSERT(nodes[0].adjacencies().find(3) != nodes[0].adjacencies().end());
    TEST_ASSERT(nodes[0].adjacencies().find(4) != nodes[0].adjacencies().end());
    // B+
    TEST_EQUALITY(nodes[3].adjacencies().size(),1);
    TEST_ASSERT(nodes[3].adjacencies().find(2) != nodes[3].adjacencies().end());
    // B++
    TEST_EQUALITY(nodes[4].adjacencies().size(),1);
    TEST_ASSERT(nodes[4].adjacencies().find(2) != nodes[4].adjacencies().end());
  }
}

// *************************************************
// Test for aliasing a field
// *************************************************
TEUCHOS_UNIT_TEST(dag, alias_field)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  using Mock = PHX::MockDAG<PHX::MyTraits::Residual,MyTraits>;

  DagManager<MyTraits> dag("alias_field");

  // Register evaluators
  {
    RCP<Mock> m = rcp(new Mock);
    m->setName("Eval_A");
    m->evaluates("A");
    m->depends("B");
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
  }

  // Alias field B to C
  {
    RCP<MDALayout<CELL,BASIS>> dl =
      rcp(new MDALayout<CELL,BASIS>("H-Grad",100,4));
    Tag<MyTraits::Residual::ScalarT> tag_b("B",dl);
    Tag<MyTraits::Residual::ScalarT> tag_c("C",dl);
    RCP<PHX::Evaluator<MyTraits>> e =
      rcp(new PHX::AliasField<MyTraits::Residual::ScalarT,MyTraits>(tag_b,tag_c));
    dag.registerEvaluator(e); // b points to c's memory
  }

  // This will fail if the logic for "B" evaluation is not set properly
  dag.sortAndOrderEvaluators();

  //out << dag << std::endl;
  dag.writeGraphvizFile("alias_field.dot",true,true,true);

  // Check that the out edges are correct.
  {
    const std::vector<PHX::DagNode<MyTraits>>& nodes = dag.getDagNodes();
    // A --> B (Alias)
    TEST_EQUALITY(nodes[0].adjacencies().size(),1);
    TEST_ASSERT(nodes[0].adjacencies().find(2) != nodes[0].adjacencies().end());
    // B (Alias) --> C
    TEST_EQUALITY(nodes[2].adjacencies().size(),1);
    TEST_ASSERT(nodes[2].adjacencies().find(1) != nodes[2].adjacencies().end());
    // C
    TEST_EQUALITY(nodes[1].adjacencies().size(),0);
  }
}

// *************************************************
// Test for aliasing a field
// *************************************************
TEUCHOS_UNIT_TEST(dag, use_range_and_unshared)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  using Mock = PHX::MockDAG<PHX::MyTraits::Residual,MyTraits>;

  // Since topological sort is not unique for a DAGs, we need to make
  // sure the dag construction enforces unique topological ordering
  // via dependencies.

  // A: eval:f1
  // B: eval:f2,dep:f1
  // C: eval:f3,dep:f2
  // D: contrib:f3
  // E: eval:f4,depf3

  DagManager<MyTraits> dag("use_range");
  {
    RCP<Mock> e = rcp(new Mock);
    e->setName("a");
    e->evaluates("f1");
    e->unshared("f1");
    dag.registerEvaluator(e);
  }
  {
    RCP<Mock> e = rcp(new Mock);
    e->setName("c");
    e->evaluates("f3");
    e->depends("f2");
    dag.registerEvaluator(e);
  }
  {
    RCP<Mock> e = rcp(new Mock);
    e->setName("e");
    e->evaluates("f4");
    e->depends("f3");
    dag.registerEvaluator(e);
  }
  {
    RCP<Mock> e = rcp(new Mock);
    e->setName("d");
    e->contributes("f3");
    dag.registerEvaluator(e);
  }
  {
    RCP<Mock> e = rcp(new Mock);
    e->setName("b");
    e->evaluates("f2");
    e->depends("f1");
    e->unshared("f2");
    e->unshared("f1");
    dag.registerEvaluator(e);
  }

  {
    RCP<MDALayout<CELL,BASIS>> dl =
      rcp(new MDALayout<CELL,BASIS>("H-Grad",100,4));
    Tag<MyTraits::Residual::ScalarT> tag_f4("f4",dl);
    dag.requireField(tag_f4);
  }

  dag.sortAndOrderEvaluators();

  std::stringstream output;
  dag.print(output);

  const auto& use_range = dag.getFieldUseRange();

  RCP<MDALayout<CELL,BASIS>> dl =
    rcp(new MDALayout<CELL,BASIS>("H-Grad",100,4));
  Tag<MyTraits::Residual::ScalarT> f1("f1",dl);
  Tag<MyTraits::Residual::ScalarT> f2("f2",dl);
  Tag<MyTraits::Residual::ScalarT> f3("f3",dl);
  Tag<MyTraits::Residual::ScalarT> f4("f4",dl);


  TEST_EQUALITY(use_range.at(f1.identifier()).first,0);
  TEST_EQUALITY(use_range.at(f1.identifier()).second,1);

  TEST_EQUALITY(use_range.at(f2.identifier()).first,1);
  TEST_EQUALITY(use_range.at(f2.identifier()).second,2);

  TEST_EQUALITY(use_range.at(f3.identifier()).first,2);
  TEST_EQUALITY(use_range.at(f3.identifier()).second,4);

  TEST_EQUALITY(use_range.at(f4.identifier()).first,4);
  TEST_EQUALITY(use_range.at(f4.identifier()).second,4);

  const auto& unshared = dag.getUnsharedFields();
  TEST_EQUALITY(unshared.size(),static_cast<std::size_t>(2));
  TEST_ASSERT(unshared.find(f1.identifier()) != unshared.end());
  TEST_ASSERT(unshared.find(f2.identifier()) != unshared.end());
  TEST_ASSERT(unshared.find(f3.identifier()) == unshared.end());
  TEST_ASSERT(unshared.find(f4.identifier()) == unshared.end());
}

// *************************************************
// Make sure we throw if no evaluator is found for
// a required field.
// *************************************************
TEUCHOS_UNIT_TEST(dag, missing_evaluators)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

  DagManager<MyTraits> dm;

  using Mock = PHX::MockDAG<PHX::MyTraits::Residual,MyTraits>;

  bool use_dynamic_layout = true;

  // Evaluates, B but we require A
  {
    RCP<Mock> e = rcp(new Mock);
    e->setName("B");
    e->contributes("B",use_dynamic_layout);
    dm.registerEvaluator(e);
  }

  RCP<PHX::Layout> dl = rcp(new PHX::Layout("H-Grad",100,4));
  PHX::Tag<PHX::MyTraits::Residual::ScalarT> tag_a("A",dl);
  dm.requireField(tag_a);

  TEST_THROW(dm.sortAndOrderEvaluators(),PHX::missing_evaluator_exception);
}

// *************************************************
// Checks that that required fields that are the starting node of a
// dfs search that have both an evaluated and contributed components
// in the DAG was failing to add the contributed part since it was
// outside of dfsVisit() call.
// *************************************************
TEUCHOS_UNIT_TEST(contrib, start_has_eval_and_contrib)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

  DagManager<MyTraits> dm;

  using Mock = PHX::MockDAG<PHX::MyTraits::Residual,MyTraits>;

  bool use_dynamic_layout = true;

  {
    RCP<Mock> e = rcp(new Mock);
    e->setName("A contributed 0");
    e->contributes("A",use_dynamic_layout);
    dm.registerEvaluator(e);
  }

  {
    RCP<Mock> e = rcp(new Mock);
    e->setName("A evaluated");
    e->evaluates("A",use_dynamic_layout);
    dm.registerEvaluator(e);
  }

  {
    RCP<Mock> e = rcp(new Mock);
    e->setName("A contributed 1");
    e->contributes("A",use_dynamic_layout);
    dm.registerEvaluator(e);
  }

  RCP<PHX::Layout> dl = rcp(new PHX::Layout("H-Grad",100,4));
  PHX::Tag<PHX::MyTraits::Residual::ScalarT> tag_a("A",dl);
  dm.requireField(tag_a);

  dm.sortAndOrderEvaluators();

  // Check the graph
  auto tags = dm.getFieldTags();
  TEST_EQUALITY(tags.size(),1);
  const auto& evaluators = dm.getEvaluatorInternalOrdering();
  TEST_EQUALITY(evaluators.size(),3);

  // Check the dot file graph
  std::stringstream output;
  dm.writeGraphviz(output,true,true);
  TEST_INEQUALITY(output.str().find("A contributed 0"),std::string::npos);
  TEST_INEQUALITY(output.str().find("A contributed 1"),std::string::npos);
  TEST_INEQUALITY(output.str().find("A evaluated"),std::string::npos);
}

// *************************************************
// Tests that a required field that is evaluated by all contributed
// fields collects all evaluators. The dfs covers all nodes internal,
// but was missing contributed on the starting required node.
// *************************************************
TEUCHOS_UNIT_TEST(contrib, start_has_contrib_only)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

  DagManager<MyTraits> dm("Residual");

  using Mock = PHX::MockDAG<PHX::MyTraits::Residual,MyTraits>;

  bool use_dynamic_layout = true;

  {
    RCP<Mock> e = rcp(new Mock);
    e->setName("A contributed 0");
    e->contributes("A",use_dynamic_layout);
    dm.registerEvaluator(e);
  }

  {
    RCP<Mock> e = rcp(new Mock);
    e->setName("A contributed 1");
    e->contributes("A",use_dynamic_layout);
    dm.registerEvaluator(e);
  }

  {
    RCP<Mock> e = rcp(new Mock);
    e->setName("A contributed 2");
    e->contributes("A",use_dynamic_layout);
    dm.registerEvaluator(e);
  }

  RCP<PHX::Layout> dl = rcp(new PHX::Layout("H-Grad",100,4));
  PHX::Tag<PHX::MyTraits::Residual::ScalarT> tag_a("A",dl);
  dm.requireField(tag_a);

  dm.sortAndOrderEvaluators();

  // Check the graph
  auto tags = dm.getFieldTags();
  TEST_EQUALITY(tags.size(),1);
  const auto& evaluators = dm.getEvaluatorInternalOrdering();
  TEST_EQUALITY(evaluators.size(),3);

  // Check the dot file graph
  std::stringstream output;
  dm.writeGraphviz(output,true,true);
  TEST_INEQUALITY(output.str().find("A contributed 0"),std::string::npos);
  TEST_INEQUALITY(output.str().find("A contributed 1"),std::string::npos);
  TEST_INEQUALITY(output.str().find("A contributed 2"),std::string::npos);
}

/*
// *************************************************
   Typical contributed field use case, residual is all contributed

      ScatterTag
          |
       Residual
      /   |   \
    Conv Diff Rxn
      \   |   /
          X

// *************************************************
*/
TEUCHOS_UNIT_TEST(contrib, basic_contrib_only)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

  DagManager<MyTraits> dm;
  using Mock = PHX::MockDAG<PHX::MyTraits::Residual,MyTraits>;
  bool use_dynamic_layout = true;

  {
    RCP<Mock> e = rcp(new Mock);
    e->setName("X");
    e->evaluates("X",use_dynamic_layout);
    dm.registerEvaluator(e);
  }

  {
    RCP<Mock> e = rcp(new Mock);
    e->setName("Convection Operator");
    e->contributes("Residual",use_dynamic_layout);
    e->depends("X",use_dynamic_layout);
    dm.registerEvaluator(e);
  }

  {
    RCP<Mock> e = rcp(new Mock);
    e->setName("Diffusion Operator");
    e->contributes("Residual",use_dynamic_layout);
    e->depends("X",use_dynamic_layout);
    dm.registerEvaluator(e);
  }

  {
    RCP<Mock> e = rcp(new Mock);
    e->setName("Reaction Operator");
    e->contributes("Residual",use_dynamic_layout);
    e->depends("X",use_dynamic_layout);
    dm.registerEvaluator(e);
  }

  {
    RCP<Mock> e = rcp(new Mock);
    e->setName("Scatter");
    // Important that this is "contributes" to catch writing graph
    // output correctly.
    e->contributes("Scatter",use_dynamic_layout);
    e->depends("Residual",use_dynamic_layout);
    dm.registerEvaluator(e);
  }

  RCP<PHX::Layout> dl = rcp(new PHX::Layout("H-Grad",100,4));
  PHX::Tag<PHX::MyTraits::Residual::ScalarT> tag_a("Scatter",dl);
  dm.requireField(tag_a);

  dm.sortAndOrderEvaluators();

  // Check the graph
  const auto& tags = dm.getFieldTags();
  TEST_EQUALITY(tags.size(),3);
  const auto& evaluators = dm.getEvaluatorInternalOrdering();
  TEST_EQUALITY(evaluators.size(),5);

  // Check the dot file graph
  std::stringstream output;
  dm.writeGraphviz(output,true,true);
  TEST_INEQUALITY(output.str().find("Scatter"),std::string::npos);
  TEST_INEQUALITY(output.str().find("Convection Operator"),std::string::npos);
  TEST_INEQUALITY(output.str().find("Diffusion Operator"),std::string::npos);
  TEST_INEQUALITY(output.str().find("Reaction Operator"),std::string::npos);
  TEST_EQUALITY(output.str().find("I am not in the graph!"),std::string::npos);
}

/*
// *************************************************
   Typical contributed field use case, residual is evaluated and contributed

      ScatterTag
          |
       Residual
       /  |  \
    Conv Diff Rxn
     | \  |  / |
      \  Init /
       \  |  /
          X

// *************************************************
*/
TEUCHOS_UNIT_TEST(contrib, basic_contrib_and_evalauted)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

  DagManager<MyTraits> dm;
  using Mock = PHX::MockDAG<PHX::MyTraits::Residual,MyTraits>;
  bool use_dynamic_layout = true;

  {
    RCP<Mock> e = rcp(new Mock);
    e->setName("X");
    e->evaluates("X",use_dynamic_layout);
    dm.registerEvaluator(e);
  }

  {
    RCP<Mock> e = rcp(new Mock);
    e->setName("Initialize");
    e->evaluates("Residual",use_dynamic_layout);
    e->depends("X",use_dynamic_layout);
    dm.registerEvaluator(e);
  }

  {
    RCP<Mock> e = rcp(new Mock);
    e->setName("Convection Operator");
    e->contributes("Residual",use_dynamic_layout);
    e->depends("X",use_dynamic_layout);
    dm.registerEvaluator(e);
  }

  {
    RCP<Mock> e = rcp(new Mock);
    e->setName("Diffusion Operator");
    e->contributes("Residual",use_dynamic_layout);
    e->depends("X",use_dynamic_layout);
    dm.registerEvaluator(e);
  }

  {
    RCP<Mock> e = rcp(new Mock);
    e->setName("Reaction Operator");
    e->contributes("Residual",use_dynamic_layout);
    e->depends("X",use_dynamic_layout);
    dm.registerEvaluator(e);
  }

  {
    RCP<Mock> e = rcp(new Mock);
    e->setName("Scatter");
    // Important that this is "contributes" to catch writing graph
    // output correctly.
    e->contributes("Scatter",use_dynamic_layout);
    e->depends("Residual",use_dynamic_layout);
    dm.registerEvaluator(e);
  }

  RCP<PHX::Layout> dl = rcp(new PHX::Layout("H-Grad",100,4));
  PHX::Tag<PHX::MyTraits::Residual::ScalarT> tag_a("Scatter",dl);
  dm.requireField(tag_a);

  dm.sortAndOrderEvaluators();

  // Check the graph
  const auto& tags = dm.getFieldTags();
  TEST_EQUALITY(tags.size(),3);
  const auto& evaluators = dm.getEvaluatorInternalOrdering();
  TEST_EQUALITY(evaluators.size(),6);

  // Check the dot file graph
  std::stringstream output;
  dm.writeGraphviz(output,true,true);
  TEST_INEQUALITY(output.str().find("Scatter"),std::string::npos);
  TEST_INEQUALITY(output.str().find("Convection Operator"),std::string::npos);
  TEST_INEQUALITY(output.str().find("Diffusion Operator"),std::string::npos);
  TEST_INEQUALITY(output.str().find("Reaction Operator"),std::string::npos);
  TEST_INEQUALITY(output.str().find("Initialize"),std::string::npos);
  TEST_EQUALITY(output.str().find("I am not in the graph!"),std::string::npos);
}
