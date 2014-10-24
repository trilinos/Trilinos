#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_DefaultComm.hpp"
#include <iostream>

// Prerequisites for testing
#include "Pike_Mock_ModelEvaluator.hpp"
#include "Pike_Solver_BlockGaussSeidel.hpp"

// Status tests to check
#include "Pike_StatusTest_Composite.hpp"
#include "Pike_StatusTest_MaxIterations.hpp"
#include "Pike_StatusTest_GlobalModelConvergence.hpp"
#include "Pike_StatusTest_LocalModelFailure.hpp"
#include "Pike_StatusTest_ScalarResponseRelativeTolerance.hpp"
#include "Pike_StatusTest_Factory.hpp"
#include "Pike_Mock_UserStatusTestFactory.hpp"

namespace pike {

  TEUCHOS_UNIT_TEST(status_test, MaxIterations)
  {

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    Teuchos::RCP<pike::MaxIterations> maxIters = Teuchos::rcp(new pike::MaxIterations);
    {
      Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList("ST");
      p->set("Maximum Iterations",4);
      maxIters->setParameterList(p);
    }

    Teuchos::RCP<pike_test::MockModelEvaluator> app1 = 
      pike_test::mockModelEvaluator(comm,"app1",pike_test::MockModelEvaluator::LOCAL_FAILURE,10,-1);

    Teuchos::RCP<pike_test::MockModelEvaluator> app2 = 
      pike_test::mockModelEvaluator(comm,"app2",pike_test::MockModelEvaluator::LOCAL_FAILURE,10,-1);

    Teuchos::RCP<pike::BlockGaussSeidel> solver = Teuchos::rcp(new pike::BlockGaussSeidel);
    app1->setSolver(solver);
    app2->setSolver(solver);
    solver->registerModelEvaluator(app1);
    solver->registerModelEvaluator(app2);
    solver->completeRegistration();
    solver->setStatusTests(maxIters);
    solver->solve();

    TEST_EQUALITY(solver->getStatus(),pike::FAILED);
    TEST_EQUALITY(solver->getNumberOfIterations(),4);

    // Test the reset on solver and status tests
    solver->reset();
    TEST_EQUALITY(solver->getStatus(), pike::UNCHECKED);
    TEST_EQUALITY(maxIters->getStatus(), pike::UNCHECKED);
  }

  TEUCHOS_UNIT_TEST(status_test, LocalModelConvergence)
  {

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    Teuchos::RCP<Teuchos::ParameterList> statusTestParams = Teuchos::parameterList("My Status Tests");
    {
      statusTestParams->set("Type","Composite OR");
      Teuchos::ParameterList& converged = statusTestParams->sublist("CONVERGED AND");
      converged.set("Type","Composite AND");
      Teuchos::ParameterList& app1LocalConv = converged.sublist("Local Convergence");      
      app1LocalConv.set("Type","Local Model Convergence");
      app1LocalConv.set("Model Name","app1");
      Teuchos::ParameterList& app2GlobalConv = converged.sublist("Converged");
      app2GlobalConv.set("Type","Global Model Convergence");
      app2GlobalConv.set("Model Name","app2");
      Teuchos::ParameterList& failure = statusTestParams->sublist("Failure 1");
      failure.set("Type","Maximum Iterations");
      failure.set("Maximum Iterations",10);
    }

    pike::StatusTestFactory stFactory;
    Teuchos::RCP<pike::StatusTest> tests = stFactory.buildStatusTests(statusTestParams);

    Teuchos::RCP<pike_test::MockModelEvaluator> app1 = 
      pike_test::mockModelEvaluator(comm,"app1",pike_test::MockModelEvaluator::LOCAL_FAILURE,7,-1);

    Teuchos::RCP<pike_test::MockModelEvaluator> app2 = 
      pike_test::mockModelEvaluator(comm,"app2",pike_test::MockModelEvaluator::GLOBAL_CONVERGENCE,10,-1);

    Teuchos::RCP<pike::BlockGaussSeidel> solver = Teuchos::rcp(new pike::BlockGaussSeidel);
    app1->setSolver(solver);
    app2->setSolver(solver);
    solver->registerModelEvaluator(app1);
    solver->registerModelEvaluator(app2);
    solver->completeRegistration();
    solver->setStatusTests(tests);
    solver->solve();

    TEST_EQUALITY(solver->getStatus(),pike::CONVERGED);
    TEST_EQUALITY(solver->getNumberOfIterations(),10);

    // Test the reset on solver and status tests
    solver->reset();
    TEST_EQUALITY(solver->getStatus(), pike::UNCHECKED);
    TEST_EQUALITY(tests->getStatus(), pike::UNCHECKED);
  }

  TEUCHOS_UNIT_TEST(status_test, LocalModelFailure)
  {

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    Teuchos::RCP<Teuchos::ParameterList> statusTestParams = Teuchos::parameterList("My Status Tests");
    {
      statusTestParams->set("Type","Composite OR");
      Teuchos::ParameterList& failure = statusTestParams->sublist("Failure 1");
      failure.set("Type","Maximum Iterations");
      failure.set("Maximum Iterations",10);
      Teuchos::ParameterList& app1LocalConv = statusTestParams->sublist("Failure 2");      
      app1LocalConv.set("Type","Local Model Failure");
      app1LocalConv.set("Model Name","app1");
      Teuchos::ParameterList& app2GlobalConv = statusTestParams->sublist("Converged");
      app2GlobalConv.set("Type","Global Model Convergence");
      app2GlobalConv.set("Model Name","app2");
    }

    pike::StatusTestFactory stFactory;
    Teuchos::RCP<pike::StatusTest> tests = stFactory.buildStatusTests(statusTestParams);

    Teuchos::RCP<pike_test::MockModelEvaluator> app1 = 
      pike_test::mockModelEvaluator(comm,"app1",pike_test::MockModelEvaluator::LOCAL_FAILURE,7,-1);

    Teuchos::RCP<pike_test::MockModelEvaluator> app2 = 
      pike_test::mockModelEvaluator(comm,"app2",pike_test::MockModelEvaluator::GLOBAL_CONVERGENCE,10,-1);

    Teuchos::RCP<pike::BlockGaussSeidel> solver = Teuchos::rcp(new pike::BlockGaussSeidel);
    app1->setSolver(solver);
    app2->setSolver(solver);
    solver->registerModelEvaluator(app1);
    solver->registerModelEvaluator(app2);
    solver->completeRegistration();
    solver->setStatusTests(tests);
    solver->solve();

    TEST_EQUALITY(solver->getStatus(),pike::FAILED);
    TEST_EQUALITY(solver->getNumberOfIterations(),7);

    // Test the reset on solver and status tests
    solver->reset();
    TEST_EQUALITY(solver->getStatus(), pike::UNCHECKED);
    TEST_EQUALITY(tests->getStatus(), pike::UNCHECKED);
  }

  TEUCHOS_UNIT_TEST(status_test, GlobalModelConvergence)
  {

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    Teuchos::RCP<Teuchos::ParameterList> statusTestParams = Teuchos::parameterList("My Status Tests");
    {
      statusTestParams->set("Type","Composite OR");
      Teuchos::ParameterList& failure = statusTestParams->sublist("Failure 1");
      failure.set("Type","Maximum Iterations");
      failure.set("Maximum Iterations",10);
      Teuchos::ParameterList& app1LocalConv = statusTestParams->sublist("Failure 2");      
      app1LocalConv.set("Type","Local Model Failure");
      app1LocalConv.set("Model Name","app1");
      Teuchos::ParameterList& app2GlobalConv = statusTestParams->sublist("Converged");
      app2GlobalConv.set("Type","Global Model Convergence");
      app2GlobalConv.set("Model Name","app2");
    }

    pike::StatusTestFactory stFactory;
    Teuchos::RCP<pike::StatusTest> tests = stFactory.buildStatusTests(statusTestParams);

    Teuchos::RCP<pike_test::MockModelEvaluator> app1 = 
      pike_test::mockModelEvaluator(comm,"app1",pike_test::MockModelEvaluator::LOCAL_FAILURE,11,-1);

    Teuchos::RCP<pike_test::MockModelEvaluator> app2 = 
      pike_test::mockModelEvaluator(comm,"app2",pike_test::MockModelEvaluator::GLOBAL_CONVERGENCE,8,-1);

    Teuchos::RCP<pike::BlockGaussSeidel> solver = Teuchos::rcp(new pike::BlockGaussSeidel);
    app1->setSolver(solver);
    app2->setSolver(solver);
    solver->registerModelEvaluator(app1);
    solver->registerModelEvaluator(app2);
    solver->completeRegistration();
    solver->setStatusTests(tests);
    solver->solve();

    TEST_EQUALITY(solver->getStatus(),pike::CONVERGED);
    TEST_EQUALITY(solver->getNumberOfIterations(),8);

    // Test the reset on solver and status tests
    solver->reset();
    TEST_EQUALITY(solver->getStatus(), pike::UNCHECKED);
    TEST_EQUALITY(tests->getStatus(), pike::UNCHECKED);
  }

  TEUCHOS_UNIT_TEST(status_test, ScalarResponseRelativeError)
  {

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    Teuchos::RCP<pike::ScalarResponseRelativeTolerance> relTol = 
      Teuchos::rcp(new pike::ScalarResponseRelativeTolerance);
    {
      Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList("ST");
      p->set("Application Name","app1");
      p->set("Response Name","Mock Response");
      p->set("Tolerance",1.0e-3);
      relTol->setParameterList(p);
    }

    Teuchos::RCP<pike_test::MockModelEvaluator> app1 = 
      pike_test::mockModelEvaluator(comm,"app1",pike_test::MockModelEvaluator::LOCAL_FAILURE,10,5);

    Teuchos::RCP<pike_test::MockModelEvaluator> app2 = 
      pike_test::mockModelEvaluator(comm,"app2",pike_test::MockModelEvaluator::LOCAL_FAILURE,10,5);

    Teuchos::RCP<pike::BlockGaussSeidel> solver = Teuchos::rcp(new pike::BlockGaussSeidel);
    app1->setSolver(solver);
    app2->setSolver(solver);
    solver->registerModelEvaluator(app1);
    solver->registerModelEvaluator(app2);
    solver->completeRegistration();
    solver->setStatusTests(relTol);
    solver->solve();

    TEST_EQUALITY(solver->getStatus(),pike::CONVERGED);
    TEST_EQUALITY(solver->getNumberOfIterations(),5);

    // Test the reset on solver and status tests
    solver->reset();
    TEST_EQUALITY(solver->getStatus(), pike::UNCHECKED);
    TEST_EQUALITY(relTol->getStatus(), pike::UNCHECKED);
    TEST_EQUALITY(solver->getNumberOfIterations(),0);
  }

  TEUCHOS_UNIT_TEST(status_test, Composite_AND)
  {

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    Teuchos::RCP<pike::ScalarResponseRelativeTolerance> relTol1 = 
      Teuchos::rcp(new pike::ScalarResponseRelativeTolerance);
    {
      Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList("ST");
      p->set("Application Name","app1");
      p->set("Response Name","Mock Response");
      p->set("Tolerance",1.0e-3);
      relTol1->setParameterList(p);
    }

    Teuchos::RCP<pike::ScalarResponseRelativeTolerance> relTol2 = 
      Teuchos::rcp(new pike::ScalarResponseRelativeTolerance);
    {
      Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList("ST");
      p->set("Application Name","app2");
      p->set("Response Name","Mock Response");
      p->set("Tolerance",1.0e-3);
      relTol2->setParameterList(p);
    }
    
    Teuchos::RCP<pike::Composite> converged = pike::composite(pike::Composite::AND);
    converged->addTest(relTol1);
    converged->addTest(relTol2);

    Teuchos::RCP<pike_test::MockModelEvaluator> app1 = 
      pike_test::mockModelEvaluator(comm,"app1",pike_test::MockModelEvaluator::LOCAL_FAILURE,10,5);

    Teuchos::RCP<pike_test::MockModelEvaluator> app2 = 
      pike_test::mockModelEvaluator(comm,"app2",pike_test::MockModelEvaluator::LOCAL_FAILURE,10,7);

    Teuchos::RCP<pike::BlockGaussSeidel> solver = Teuchos::rcp(new pike::BlockGaussSeidel);
    app1->setSolver(solver);
    app2->setSolver(solver);
    solver->registerModelEvaluator(app1);
    solver->registerModelEvaluator(app2);
    solver->completeRegistration();
    solver->setStatusTests(converged);
    solver->solve();

    TEST_EQUALITY(solver->getStatus(),pike::CONVERGED);
    TEST_EQUALITY(solver->getNumberOfIterations(),7);

    // Test solver and status test reset
    TEST_EQUALITY(solver->getStatus(),pike::CONVERGED);
    solver->reset();
    TEST_EQUALITY(converged->getStatus(),pike::UNCHECKED);
  }

  TEUCHOS_UNIT_TEST(status_test, Composite_OR)
  {

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    Teuchos::RCP<pike::ScalarResponseRelativeTolerance> relTol1 = 
      Teuchos::rcp(new pike::ScalarResponseRelativeTolerance);
    {
      Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList("ST");
      p->set("Application Name","app1");
      p->set("Response Name","Mock Response");
      p->set("Tolerance",1.0e-3);
      relTol1->setParameterList(p);
    }

    Teuchos::RCP<pike::ScalarResponseRelativeTolerance> relTol2 = 
      Teuchos::rcp(new pike::ScalarResponseRelativeTolerance);
    {
      Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList("ST");
      p->set("Application Name","app2");
      p->set("Response Name","Mock Response");
      p->set("Tolerance",1.0e-3);
      relTol2->setParameterList(p);
    }

    Teuchos::RCP<pike::Composite> converged = pike::composite(pike::Composite::OR);
    converged->addTest(relTol1);
    converged->addTest(relTol2);

    Teuchos::RCP<pike_test::MockModelEvaluator> app1 = 
      pike_test::mockModelEvaluator(comm,"app1",pike_test::MockModelEvaluator::LOCAL_FAILURE,10,5);

    Teuchos::RCP<pike_test::MockModelEvaluator> app2 = 
      pike_test::mockModelEvaluator(comm,"app2",pike_test::MockModelEvaluator::LOCAL_FAILURE,10,7);

    Teuchos::RCP<pike::BlockGaussSeidel> solver = Teuchos::rcp(new pike::BlockGaussSeidel);
    app1->setSolver(solver);
    app2->setSolver(solver);
    solver->registerModelEvaluator(app1);
    solver->registerModelEvaluator(app2);
    solver->completeRegistration();
    solver->setStatusTests(converged);
    solver->solve();

    TEST_EQUALITY(solver->getStatus(),pike::CONVERGED);
    TEST_EQUALITY(solver->getNumberOfIterations(),5);

    // Test solver and status test reset
    TEST_EQUALITY(solver->getStatus(),pike::CONVERGED);
    solver->reset();
    TEST_EQUALITY(converged->getStatus(),pike::UNCHECKED);
  }

  TEUCHOS_UNIT_TEST(status_test, Composite_NESTED)
  {

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    Teuchos::RCP<pike::ScalarResponseRelativeTolerance> relTol1 = 
      Teuchos::rcp(new pike::ScalarResponseRelativeTolerance);
    {
      Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList("ST");
      p->set("Application Name","app1");
      p->set("Response Name","Mock Response");
      p->set("Tolerance",1.0e-3);
      relTol1->setParameterList(p);
    }

    Teuchos::RCP<pike::ScalarResponseRelativeTolerance> relTol2 = 
      Teuchos::rcp(new pike::ScalarResponseRelativeTolerance);
    {
      Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList("ST");
      p->set("Application Name","app2");
      p->set("Response Name","Mock Response");
      p->set("Tolerance",1.0e-3);
      relTol2->setParameterList(p);
    }

    Teuchos::RCP<pike::Composite> converged = pike::composite(pike::Composite::AND);
    converged->addTest(relTol1);
    converged->addTest(relTol2);

    Teuchos::RCP<pike::MaxIterations> maxIters = Teuchos::rcp(new pike::MaxIterations);
    {
      Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList("ST");
      p->set("Maximum Iterations",6);
      maxIters->setParameterList(p);
    }

    Teuchos::RCP<pike::Composite> tests = pike::composite(pike::Composite::OR);
    tests->addTest(converged);
    tests->addTest(maxIters);

    Teuchos::RCP<pike_test::MockModelEvaluator> app1 = 
      pike_test::mockModelEvaluator(comm,"app1",pike_test::MockModelEvaluator::LOCAL_FAILURE,10,5);

    Teuchos::RCP<pike_test::MockModelEvaluator> app2 = 
      pike_test::mockModelEvaluator(comm,"app2",pike_test::MockModelEvaluator::LOCAL_FAILURE,10,7);

    Teuchos::RCP<pike::BlockGaussSeidel> solver = Teuchos::rcp(new pike::BlockGaussSeidel);
    app1->setSolver(solver);
    app2->setSolver(solver);
    solver->registerModelEvaluator(app1);
    solver->registerModelEvaluator(app2);
    solver->completeRegistration();
    solver->setStatusTests(tests);
    solver->solve();

    TEST_EQUALITY(solver->getStatus(),pike::FAILED);
    TEST_EQUALITY(solver->getNumberOfIterations(),6); 

    // Test solver and status test reset
    TEST_EQUALITY(solver->getStatus(),pike::FAILED);
    solver->reset();
    TEST_EQUALITY(tests->getStatus(),pike::UNCHECKED);
 }

  TEUCHOS_UNIT_TEST(status_test, Factory)
  {

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    Teuchos::RCP<Teuchos::ParameterList> statusTestParams = Teuchos::parameterList("My Status Tests");
    {
      statusTestParams->set("Type","Composite OR");
      Teuchos::ParameterList& failure = statusTestParams->sublist("Failure");
      failure.set("Type","Maximum Iterations");
      failure.set("Maximum Iterations",6);
      Teuchos::ParameterList& app1LocalConv = statusTestParams->sublist("Failure 2");      
      app1LocalConv.set("Type","Local Model Failure");
      app1LocalConv.set("Model Name","app1");
      Teuchos::ParameterList& app2LocalConv = statusTestParams->sublist("Failure 3");      
      app2LocalConv.set("Type","Local Model Failure");
      app2LocalConv.set("Model Name","app2");
      Teuchos::ParameterList& converged = statusTestParams->sublist("Converged");
      converged.set("Type","Composite AND");
      Teuchos::ParameterList& relTolApp1 = converged.sublist("App 1");
      relTolApp1.set("Type","Scalar Response Relative Tolerance");
      relTolApp1.set("Application Name","app1");
      relTolApp1.set("Response Name","Mock Response");
      relTolApp1.set("Tolerance",1.0e-3);
      Teuchos::ParameterList& relTolApp2 = converged.sublist("App 2"); 
      relTolApp2.set("Type","Scalar Response Relative Tolerance");
      relTolApp2.set("Application Name","app2");
      relTolApp2.set("Response Name","Mock Response");
      relTolApp2.set("Tolerance",1.0e-3);
      Teuchos::ParameterList& app1GlobalConv = converged.sublist("Global App 1");
      app1GlobalConv.set("Type","Global Model Convergence");
      app1GlobalConv.set("Model Name","app1");
      Teuchos::ParameterList& app2GlobalConv = converged.sublist("Global App 2");
      app2GlobalConv.set("Type","Global Model Convergence");
      app2GlobalConv.set("Model Name","app2");
    }

    pike::StatusTestFactory stFactory;
    Teuchos::RCP<pike::StatusTest> tests = stFactory.buildStatusTests(statusTestParams);
    TEST_ASSERT(nonnull(tests));

    Teuchos::RCP<pike_test::MockModelEvaluator> app1 = 
      pike_test::mockModelEvaluator(comm,"app1",pike_test::MockModelEvaluator::LOCAL_FAILURE,10,5);

    Teuchos::RCP<pike_test::MockModelEvaluator> app2 = 
      pike_test::mockModelEvaluator(comm,"app2",pike_test::MockModelEvaluator::LOCAL_FAILURE,10,7);

    Teuchos::RCP<pike::BlockGaussSeidel> solver = Teuchos::rcp(new pike::BlockGaussSeidel);
    app1->setSolver(solver);
    app2->setSolver(solver);
    solver->registerModelEvaluator(app1);
    solver->registerModelEvaluator(app2);
    solver->completeRegistration();
    solver->setStatusTests(tests);
    solver->solve();

    TEST_EQUALITY(solver->getStatus(),pike::FAILED);
    TEST_EQUALITY(solver->getNumberOfIterations(),6);
  }

  TEUCHOS_UNIT_TEST(status_test, UserFactory)
  {

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    Teuchos::RCP<Teuchos::ParameterList> statusTestParams = Teuchos::parameterList("My Status Tests");
    {
      statusTestParams->set("Type","Composite OR");
      Teuchos::ParameterList& failure = statusTestParams->sublist("Failure");
      failure.set("Type","My Super Special MaxIterations");
      failure.set("Maximum Iterations",6);
      Teuchos::ParameterList& app1LocalConv = statusTestParams->sublist("Failure 2");      
      app1LocalConv.set("Type","Local Model Failure");
      app1LocalConv.set("Model Name","app1");
      Teuchos::ParameterList& app2LocalConv = statusTestParams->sublist("Failure 3");      
      app2LocalConv.set("Type","Local Model Failure");
      app2LocalConv.set("Model Name","app2");
      Teuchos::ParameterList& failure2 = statusTestParams->sublist("Failure 4");
      failure2.set("Type","My Other Super Special MaxIterations");
      failure2.set("Maximum Iterations",6);
      Teuchos::ParameterList& converged = statusTestParams->sublist("Converged");
      converged.set("Type","Composite AND");
      Teuchos::ParameterList& relTolApp1 = converged.sublist("App 1");
      relTolApp1.set("Type","Scalar Response Relative Tolerance");
      relTolApp1.set("Application Name","app1");
      relTolApp1.set("Response Name","Mock Response");
      relTolApp1.set("Tolerance",1.0e-3);
      Teuchos::ParameterList& relTolApp2 = converged.sublist("App 2"); 
      relTolApp2.set("Type","Scalar Response Relative Tolerance");
      relTolApp2.set("Application Name","app2");
      relTolApp2.set("Response Name","Mock Response");
      relTolApp2.set("Tolerance",1.0e-3);
      Teuchos::ParameterList& app1GlobalConv = converged.sublist("Global App 1");
      app1GlobalConv.set("Type","Global Model Convergence");
      app1GlobalConv.set("Model Name","app1");
      Teuchos::ParameterList& app2GlobalConv = converged.sublist("Global App 2");
      app2GlobalConv.set("Type","Global Model Convergence");
      app2GlobalConv.set("Model Name","app2");
    }

    pike::StatusTestFactory stFactory;
    Teuchos::RCP<pike_test::UserStatusTestFactory> user1 = 
      Teuchos::rcp(new pike_test::UserStatusTestFactory("My Super Special MaxIterations"));
    Teuchos::RCP<pike_test::UserStatusTestFactory> user2 = 
      Teuchos::rcp(new pike_test::UserStatusTestFactory("My Other Super Special MaxIterations"));
    stFactory.addFactory(user1);
    stFactory.addFactory(user2);
    Teuchos::RCP<pike::StatusTest> tests = stFactory.buildStatusTests(statusTestParams);
    TEST_ASSERT(nonnull(tests));

    Teuchos::RCP<pike_test::MockModelEvaluator> app1 = 
      pike_test::mockModelEvaluator(comm,"app1",pike_test::MockModelEvaluator::LOCAL_FAILURE,10,5);

    Teuchos::RCP<pike_test::MockModelEvaluator> app2 = 
      pike_test::mockModelEvaluator(comm,"app2",pike_test::MockModelEvaluator::LOCAL_FAILURE,10,7);

    Teuchos::RCP<pike::BlockGaussSeidel> solver = Teuchos::rcp(new pike::BlockGaussSeidel);
    app1->setSolver(solver);
    app2->setSolver(solver);
    solver->registerModelEvaluator(app1);
    solver->registerModelEvaluator(app2);
    solver->completeRegistration();
    solver->setStatusTests(tests);
    solver->solve();

    TEST_EQUALITY(solver->getStatus(),pike::FAILED);
    TEST_EQUALITY(solver->getNumberOfIterations(),6);
  }

}
