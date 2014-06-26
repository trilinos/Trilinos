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
  }

}
