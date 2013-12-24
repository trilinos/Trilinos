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
#include "Pike_StatusTest_ScalarResponseRelativeTolerance.hpp"

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
  }

}
