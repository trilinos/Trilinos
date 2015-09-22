#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_DefaultComm.hpp"
#include <iostream>

// Prerequisites for testing
#include "Pike_Mock_ModelEvaluator.hpp"
#include "Pike_Mock_DataTransfer.hpp"
#include "Pike_Solver_BlockGaussSeidel.hpp"

// Status tests to check
#include "Pike_StatusTest_Composite.hpp"
#include "Pike_StatusTest_MaxIterations.hpp"
#include "Pike_StatusTest_ScalarResponseRelativeTolerance.hpp"

// Observers and decorators
#include "Pike_SolverObserver_Logger.hpp"
#include "Pike_BlackBoxModelEvaluator_Logger.hpp"
#include "Pike_DataTransfer_Logger.hpp"

namespace pike {

  TEUCHOS_UNIT_TEST(decorators, loggers)
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
    Teuchos::RCP<pike::ModelEvaluatorLogger> app1Logged = pike::modelEvaluatorLogger(app1);

    Teuchos::RCP<pike_test::MockModelEvaluator> app2 = 
      pike_test::mockModelEvaluator(comm,"app2",pike_test::MockModelEvaluator::LOCAL_FAILURE,10,7);
    Teuchos::RCP<pike::ModelEvaluatorLogger> app2Logged = pike::modelEvaluatorLogger(app2);

    std::vector<std::string> app1StringVec;
    app1StringVec.push_back("app1");
    std::vector<std::string> app2StringVec;
    app2StringVec.push_back("app2");
    Teuchos::RCP<pike_test::MockDataTransfer> trans1To2 = 
      pike_test::mockDataTransfer(comm,"app1 to app2",app1StringVec,app2StringVec);
    Teuchos::RCP<pike_test::MockDataTransfer> trans2To1 = 
      pike_test::mockDataTransfer(comm,"app2 to app1",app2StringVec,app1StringVec);
    Teuchos::RCP<pike::DataTransferLogger> trans1To2Logged =
      pike::dataTransferLogger(trans1To2);
    Teuchos::RCP<pike::DataTransferLogger> trans2To1Logged =
      pike::dataTransferLogger(trans2To1);

    Teuchos::RCP<pike::LoggerObserver> logger = pike::loggerObserver();
    Teuchos::RCP<std::vector<std::string> > log = Teuchos::rcp(new std::vector<std::string>);
    logger->setLog(log);
    app1Logged->setLog(log);
    app2Logged->setLog(log);
    trans1To2Logged->setLog(log);
    trans2To1Logged->setLog(log);
    
    Teuchos::RCP<pike::BlockGaussSeidel> solver = Teuchos::rcp(new pike::BlockGaussSeidel);
    app1->setSolver(solver);
    app2->setSolver(solver);
    solver->registerModelEvaluator(app1Logged);
    solver->registerModelEvaluator(app2Logged);
    solver->registerDataTransfer(trans1To2Logged);
    solver->registerDataTransfer(trans2To1Logged);
    solver->completeRegistration();
    solver->setStatusTests(tests);
    solver->addObserver(logger);
    solver->initialize();
    solver->solve();
    solver->finalize();

    TEST_EQUALITY(solver->getStatus(),pike::CONVERGED);
    TEST_EQUALITY(solver->getNumberOfIterations(),7);

    TEST_EQUALITY(log->size(), 63);

    for (std::vector<std::string>::const_iterator l=log->begin(); l != log->end(); ++l)
      out << *l << std::endl;

    // Test extra logger functions for coverage testing.  These are
    // just pass through functions to the underlying ME
    Teuchos::RCP<std::vector<std::string> > logRCP = app1Logged->getNonConstLog();
    Teuchos::RCP<const std::vector<std::string> > nonconstLogRCP = app1Logged->getLog();
    TEST_EQUALITY(app1Logged->isLocallyConverged(), true);
    TEST_EQUALITY(app1Logged->isGloballyConverged(), false);
    TEST_EQUALITY(app1Logged->getResponseName(0), "Mock Response");
    TEST_EQUALITY(app1Logged->supportsResponse("Mock Response"), true);
    TEST_EQUALITY(app1Logged->getNumberOfResponses(),1);
    TEST_EQUALITY(app1Logged->supportsParameter("Mock Parameter"), true);
    TEST_EQUALITY(app1Logged->getNumberOfParameters(), 1);
    TEST_EQUALITY(app1Logged->getParameterIndex("Mock Parameter"), 0);
    Teuchos::Array<double> a(1);
    app1Logged->setParameter(0,a);
    TEST_EQUALITY(app1Logged->isTransient(), false);
    TEST_EQUALITY(app1Logged->getCurrentTime(), 0.0);
    TEST_EQUALITY(app1Logged->getTentativeTime(), 0.0);
    TEST_EQUALITY(app1Logged->solvedTentativeStep(), false);
    TEST_EQUALITY(app1Logged->getCurrentTimeStepSize(), 1.0);
    TEST_EQUALITY(app1Logged->getDesiredTimeStepSize(), 1.0);
    TEST_EQUALITY(app1Logged->getMaxTimeStepSize(), 1.0);
    app1Logged->setNextTimeStepSize(1.0);
    app1Logged->acceptTimeStep();
    // data transfer logger
    trans1To2Logged->getLog();
    trans1To2Logged->getNonConstLog();
    TEST_ASSERT(trans1To2Logged->transferSucceeded());
    TEST_EQUALITY(*(trans1To2Logged->getSourceModelNames().begin()), "app1");
    TEST_EQUALITY(*(trans1To2Logged->getTargetModelNames().begin()), "app2");
  }

}
