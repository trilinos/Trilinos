#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Pike_BlackBox_config.hpp"
#include "Pike_MultiphysicsDistributor.hpp"

// Solvers
#include "Pike_Solver_BlockGaussSeidel.hpp"
#include "Pike_Solver_BlockJacobi.hpp"
#include "Pike_Solver_Factory.hpp"
#include "Pike_Mock_UserSolverFactory.hpp"

// Models
#include "Pike_Rxn_ModelEvaluator_All.hpp"
#include "Pike_Rxn_ModelEvaluator_SingleEq1.hpp"
#include "Pike_Rxn_ModelEvaluator_SingleEq2.hpp"
#include "Pike_Rxn_ModelEvaluator_SingleEq3.hpp"
#include "Pike_Rxn_DataTransfer_Eq1ToEq2.hpp"
#include "Pike_Rxn_DataTransfer_Eq1ToEq3.hpp"
#include "Pike_BlackBoxModelEvaluator_SolverAdapter.hpp"

// Status tests
#include "Pike_StatusTest_Factory.hpp"
#include "Pike_StatusTest_Composite.hpp"
#include "Pike_StatusTest_MaxIterations.hpp"
#include "Pike_StatusTest_ScalarResponseRelativeTolerance.hpp"

#include <vector>
#include <iostream>

namespace pike_test {

  /* Transient unit tests for chemical reaction in a hierarchical
     solve

     This will demonstrate order of accuracy for split system.
  */

  double evaluateOrder(const std::vector<std::pair<double,double>>& error);

  void runTransientSolve(const double& startTime,
			 const double& endTime,
			 const double& stepSize,
			 RxnAll& rxnME,
			 std::vector<std::pair<double,double>>& error);

  void runTransientSolveSingleME(const double& startTime,
				 const double& endTime,
				 const double& stepSize,
				 pike::SolverAdapterModelEvaluator& rxnME,
				 pike_test::RxnSingleEq1& rxnME1,
				 pike_test::RxnSingleEq2& rxnME2,
				 pike_test::RxnSingleEq3& rxnME3,
				 std::vector<std::pair<double,double>>& error);

  TEUCHOS_UNIT_TEST(rxn, monolithic)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using MPD = pike::MultiphysicsDistributor;

    RCP<const Teuchos::Comm<int> > globalComm = Teuchos::DefaultComm<int>::getComm();
    
    RCP<MPD> mpd = rcp(new MPD("global"));
    RCP<RxnAll> rxnME = rxnAll(mpd);

    // Drive the time stepping algorithm
    
    const double startTime = 0.0;
    const double endTime = 0.1;

    std::vector<std::pair<double,double>> error;
    runTransientSolve(startTime,endTime,1e-1,*rxnME,error);
    runTransientSolve(startTime,endTime,5e-2,*rxnME,error);
    runTransientSolve(startTime,endTime,1e-2,*rxnME,error);
    runTransientSolve(startTime,endTime,5e-3,*rxnME,error);
    runTransientSolve(startTime,endTime,1e-3,*rxnME,error);

    std::cout << std::endl;
    for (auto&& e : error)
      std::cout << "h=" << e.first << ", error=" << e.second << std::endl;

    // compute order
    double order = evaluateOrder(error);

    TEST_ASSERT( std::abs(order-4.0) < 1.0e-1);
  }

  TEUCHOS_UNIT_TEST(rxn, hierarchic)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using MPD = pike::MultiphysicsDistributor;

    RCP<const Teuchos::Comm<int> > globalComm = Teuchos::DefaultComm<int>::getComm();
    
    RCP<MPD> mpd = rcp(new MPD("global"));
    auto eq1 = mpd->addApplication("Eq1",0,0);
    auto eq2 = mpd->addApplication("Eq2",0,0);
    auto eq3 = mpd->addApplication("Eq3",0,0);
    mpd->addTransfer("Eq1->Eq2",eq1,eq2);
    mpd->addTransfer("Eq1->Eq3",eq1,eq3);
    mpd->setup(globalComm,true);
    RCP<RxnSingleEq1> rxnME1 = rxnSingleEq1(mpd);
    RCP<RxnSingleEq2> rxnME2 = rxnSingleEq2(mpd);
    RCP<RxnSingleEq3> rxnME3 = rxnSingleEq3(mpd);
    RCP<RxnDT1To2> rxnDT1t2 = rxnDT1To2(mpd,rxnME1,rxnME2);
    RCP<RxnDT1To3> rxnDT1t3 = rxnDT1To3(mpd,rxnME1,rxnME3);

    Teuchos::RCP<Teuchos::ParameterList> solverParams = Teuchos::parameterList();
    Teuchos::updateParametersFromXmlFileAndBroadcast("rxn_params.xml",
						     solverParams.ptr(),
						     *globalComm);
    pike::SolverFactory factory;
    RCP<pike::Solver> solver = factory.buildSolver(solverParams);
    solver->registerModelEvaluator(rxnME1);
    solver->registerModelEvaluator(rxnME2);
    solver->registerModelEvaluator(rxnME3);
    solver->registerDataTransfer(rxnDT1t2);
    solver->registerDataTransfer(rxnDT1t3);
    solver->completeRegistration();
    Teuchos::RCP<Teuchos::ParameterList> statusTestParams = Teuchos::sublist(solverParams,"My Stopping Criteria");
    pike::StatusTestFactory stFactory;
    Teuchos::RCP<pike::StatusTest> tests = stFactory.buildStatusTests(statusTestParams);
    solver->setStatusTests(tests);

    // wrap solver in ME
    RCP<pike::SolverAdapterModelEvaluator> rxnME = rcp(new pike::SolverAdapterModelEvaluator("Rxn"));
    rxnME->setSolver(solver);

    // Drive the time stepping algorithm
    
    const double startTime = 0.0;
    const double endTime = 0.1;

    std::vector<std::pair<double,double>> error;
    runTransientSolveSingleME(startTime,endTime,1e-1,*rxnME,*rxnME1,*rxnME2,*rxnME3,error);
    runTransientSolveSingleME(startTime,endTime,5e-2,*rxnME,*rxnME1,*rxnME2,*rxnME3,error);
    runTransientSolveSingleME(startTime,endTime,1e-2,*rxnME,*rxnME1,*rxnME2,*rxnME3,error);
    runTransientSolveSingleME(startTime,endTime,5e-3,*rxnME,*rxnME1,*rxnME2,*rxnME3,error);
    runTransientSolveSingleME(startTime,endTime,1e-3,*rxnME,*rxnME1,*rxnME2,*rxnME3,error);

    std::cout << std::endl;
    for (auto&& e : error)
      std::cout << "h=" << e.first << ", error=" << e.second << std::endl;

    // compute order
    double order = evaluateOrder(error);

    TEST_ASSERT( std::abs(order-4.0) < 1.0e-1);
  }

  double evaluateOrder(const std::vector<std::pair<double,double>>& error)
  {
    const std::size_t size = error.size();
    std::vector<double> log_x(size);
    std::vector<double> log_y(size);
    double avg_log_x = 0.0;
    double avg_log_y = 0.0;
    for (std::size_t i=0; i < size; ++i) {
      log_x[i] = std::log(error[i].first); 
      log_y[i] = std::log(error[i].second);
      avg_log_x += log_x[i];
      avg_log_y += log_y[i];
    }
    avg_log_x /= static_cast<double>(log_x.size());
    avg_log_y /= static_cast<double>(log_y.size());

    double sd_x = 0.0;
    double sd_y = 0.0;
    for (std::size_t i=0; i < log_x.size(); ++i) {
      sd_x += (log_x[i]-avg_log_x)*(log_x[i]-avg_log_x);
      sd_y += (log_y[i]-avg_log_y)*(log_y[i]-avg_log_y);
    }

    sd_x = std::sqrt(sd_x/static_cast<double>(log_x.size()));
    sd_y = std::sqrt(sd_y/static_cast<double>(log_y.size()));

    std::cout << "sd_x=" << sd_x << std::endl;
    std::cout << "sd_y=" << sd_y << std::endl;

    double sum_dev_x_dev_y = 0.0;
    double sum_dev_x2 = 0.0;
    double sum_dev_y2 = 0.0;
    for (std::size_t i=0; i < log_x.size(); ++i) {
      sum_dev_x_dev_y += (log_x[i]-avg_log_x)*(log_y[i]-avg_log_y);
      sum_dev_x2 += (log_x[i]-avg_log_x)*(log_x[i]-avg_log_x);
      sum_dev_y2 += (log_y[i]-avg_log_y)*(log_y[i]-avg_log_y);
    }
    double r = sum_dev_x_dev_y / std::sqrt(sum_dev_x2*sum_dev_y2);
    double slope = r * sd_y / sd_x;

    std::cout << "order = " << slope << std::endl;

    return slope;
  }

  void runTransientSolve(const double& startTime,
			 const double& endTime,
			 const double& stepSize,
			 RxnAll& rxnME,
			 std::vector<std::pair<double,double>>& error)
  {
    int numSteps = (endTime - startTime) / stepSize;
    TEUCHOS_ASSERT(std::fabs(numSteps*stepSize - (endTime-startTime) ) < 1.0e-10);
    rxnME.reset();
    rxnME.setNextTimeStepSize(stepSize);
    for (int i=0;i<numSteps;++i) {
      rxnME.solve();
      rxnME.acceptTimeStep();
    }
    error.push_back(std::make_pair(stepSize,rxnME.evaluateError()));
  }

  void runTransientSolveSingleME(const double& startTime,
				 const double& endTime,
				 const double& stepSize,
				 pike::SolverAdapterModelEvaluator& rxnME,
				 pike_test::RxnSingleEq1& rxnME1,
				 pike_test::RxnSingleEq2& rxnME2,
				 pike_test::RxnSingleEq3& rxnME3,
				 std::vector<std::pair<double,double>>& error)
  {
    int numSteps = (endTime - startTime) / stepSize;
    TEUCHOS_ASSERT(std::fabs(numSteps*stepSize - (endTime-startTime) ) < 1.0e-10);
    rxnME1.reset();
    rxnME2.reset();
    rxnME3.reset();
    rxnME.setNextTimeStepSize(stepSize);
    for (int i=0;i<numSteps;++i) {
      rxnME.solve();
      rxnME.acceptTimeStep();
    }
    
    double error_value = 0.0;
    error_value += rxnME1.evaluateError()*rxnME1.evaluateError();
    error_value += rxnME2.evaluateError()*rxnME2.evaluateError();
    error_value += rxnME3.evaluateError()*rxnME3.evaluateError();
    error.push_back(std::make_pair(stepSize,std::sqrt(error_value)));
  }

}
