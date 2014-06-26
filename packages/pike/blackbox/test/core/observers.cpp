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

// Observers
#include "Pike_SolverObserver_Logger.hpp"

namespace pike {

  TEUCHOS_UNIT_TEST(observers, logger)
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
      p->set("Maximum Iterations",10);
      maxIters->setParameterList(p);
    }

    Teuchos::RCP<pike::Composite> tests = pike::composite(pike::Composite::OR);
    tests->addTest(converged);
    tests->addTest(maxIters);

    Teuchos::RCP<pike_test::MockModelEvaluator> app1 = 
      pike_test::mockModelEvaluator(comm,"app1",pike_test::MockModelEvaluator::LOCAL_FAILURE,10,5);

    Teuchos::RCP<pike_test::MockModelEvaluator> app2 = 
      pike_test::mockModelEvaluator(comm,"app2",pike_test::MockModelEvaluator::LOCAL_FAILURE,10,7);

    Teuchos::RCP<pike::LoggerObserver> logger = pike::loggerObserver();

    Teuchos::RCP<pike::BlockGaussSeidel> solver = Teuchos::rcp(new pike::BlockGaussSeidel);
    app1->setSolver(solver);
    app2->setSolver(solver);
    solver->registerModelEvaluator(app1);
    solver->registerModelEvaluator(app2);
    solver->completeRegistration();
    solver->setStatusTests(tests);
    solver->addObserver(logger);
    solver->initialize();
    solver->solve();
    solver->finalize();

    TEST_EQUALITY(solver->getStatus(),pike::CONVERGED);
    TEST_EQUALITY(solver->getNumberOfIterations(),7);

    Teuchos::RCP<const std::vector<std::string> > log = logger->getLog();
    
    TEST_EQUALITY(log->size(), 19);

    for (std::vector<std::string>::const_iterator l=log->begin(); l != log->end(); ++l)
      out << *l << std::endl;

    Teuchos::RCP<std::vector<std::string> > nonConstLog = logger->getNonConstLog();

    nonConstLog->push_back("test");
    TEST_EQUALITY(log->size(), 20);
  }

}
