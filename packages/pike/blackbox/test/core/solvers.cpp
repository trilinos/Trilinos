#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Pike_BlackBox_config.hpp"
#include <iostream>

// Solvers
#include "Pike_Solver_BlockGaussSeidel.hpp"
#include "Pike_Solver_BlockJacobi.hpp"

// Models
#include "Pike_LinearHeatConduction_ModelEvaluator.hpp"
#include "Pike_LinearHeatConduction_DataTransfer.hpp"

// Status tests
#include "Pike_StatusTest_Composite.hpp"
#include "Pike_StatusTest_MaxIterations.hpp"
#include "Pike_StatusTest_ScalarResponseRelativeTolerance.hpp"

namespace pike_test {

  TEUCHOS_UNIT_TEST(solvers, ostream_overload)
  {
    pike::BlockGaussSeidel solver;

    std::cout << solver << std::endl;
  }

  TEUCHOS_UNIT_TEST(solvers, block_gauss_seidel)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    Teuchos::RCP<const Teuchos::Comm<int> > globalComm = Teuchos::DefaultComm<int>::getComm();

    RCP<LinearHeatConductionModelEvaluator> leftWall;
    {
      leftWall = linearHeatConductionModelEvaluator(globalComm,"left wall",pike_test::LinearHeatConductionModelEvaluator::T_RIGHT_IS_RESPONSE);
      leftWall->set_T_left(7.0);
      leftWall->set_T_right(5.0);  // final solution is 6.0
      leftWall->set_k(1.0);
      leftWall->set_q(1.0);
    }

    RCP<LinearHeatConductionModelEvaluator> middleWall;
    {
      middleWall = linearHeatConductionModelEvaluator(globalComm,"middle wall",pike_test::LinearHeatConductionModelEvaluator::T_RIGHT_IS_RESPONSE);
      middleWall->set_T_left(6.0);
      middleWall->set_T_right(3.0);  // final solution is 4.0
      middleWall->set_k(1.0/2.0);
      middleWall->set_q(1.0);
    }

    RCP<LinearHeatConductionModelEvaluator> rightWall;
    {
      rightWall =linearHeatConductionModelEvaluator(globalComm,"right wall",pike_test::LinearHeatConductionModelEvaluator::Q_IS_RESPONSE);
      rightWall->set_T_left(4.0);
      rightWall->set_T_right(1.0);
      rightWall->set_k(1.0/3.0);
      rightWall->set_q(1.5); // final solution is 1.0
    }

    RCP<LinearHeatConductionDataTransfer> transferQ = 
      linearHeatConductionDataTransfer(globalComm,"tranfers q: right->{left,middle}",pike_test::LinearHeatConductionDataTransfer::TRANSFER_Q);
    transferQ->setSource(rightWall);
    transferQ->addTarget(leftWall);
    transferQ->addTarget(middleWall);
    
    RCP<LinearHeatConductionDataTransfer> transferLeftToMiddle =
      linearHeatConductionDataTransfer(globalComm,"tranfer T: left->middle",pike_test::LinearHeatConductionDataTransfer::TRANSFER_T);
    transferLeftToMiddle->setSource(leftWall);
    transferLeftToMiddle->addTarget(middleWall);

    RCP<LinearHeatConductionDataTransfer> transferMiddleToRight =
      linearHeatConductionDataTransfer(globalComm,"tranfer T: middle->right",pike_test::LinearHeatConductionDataTransfer::TRANSFER_T);
    transferMiddleToRight->setSource(middleWall);
    transferMiddleToRight->addTarget(rightWall);

    Teuchos::RCP<pike::Composite> status = pike::composite(pike::Composite::OR);
    Teuchos::RCP<pike::MaxIterations> maxIters =
      Teuchos::rcp(new pike::MaxIterations(20));
    Teuchos::RCP<pike::Composite> convergedTests = 
      pike::composite(pike::Composite::AND);
    Teuchos::RCP<pike::ScalarResponseRelativeTolerance> t1 = 
      Teuchos::rcp(new pike::ScalarResponseRelativeTolerance);
    {
      Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
      p->set("Application Name","left wall");
      p->set("Response Name","T_right");
      p->set("Tolerance",1.0e-5);
      t1->setParameterList(p);
    }
    Teuchos::RCP<pike::ScalarResponseRelativeTolerance> t2 = 
      Teuchos::rcp(new pike::ScalarResponseRelativeTolerance);
    {
      Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
      p->set("Application Name","middle wall");
      p->set("Response Name","T_right");
      p->set("Tolerance",1.0e-5);
      t2->setParameterList(p);
    }
    Teuchos::RCP<pike::ScalarResponseRelativeTolerance> t3 = 
      Teuchos::rcp(new pike::ScalarResponseRelativeTolerance);
    {
      Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
      p->set("Application Name","right wall");
      p->set("Response Name","q");
      p->set("Tolerance",1.0e-5);
      t3->setParameterList(p);
    }
    convergedTests->addTest(t1);
    convergedTests->addTest(t2);
    convergedTests->addTest(t3);
    status->addTest(maxIters);
    status->addTest(convergedTests);

    pike::BlockGaussSeidel solver;
    solver.registerModelEvaluator(leftWall);
    solver.registerModelEvaluator(middleWall);
    solver.registerModelEvaluator(rightWall);
    solver.registerDataTransfer(transferQ);
    solver.registerDataTransfer(transferLeftToMiddle);
    solver.registerDataTransfer(transferMiddleToRight);
    solver.completeRegistration();
    solver.setStatusTests(status);
    solver.solve();

    TEST_EQUALITY(solver.getNumberOfIterations(),10);
    TEST_EQUALITY(solver.getStatus(),pike::CONVERGED);
  }

  TEUCHOS_UNIT_TEST(solvers, block_jacobi)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    Teuchos::RCP<const Teuchos::Comm<int> > globalComm = Teuchos::DefaultComm<int>::getComm();

    RCP<LinearHeatConductionModelEvaluator> leftWall;
    {
      leftWall = linearHeatConductionModelEvaluator(globalComm,"left wall",pike_test::LinearHeatConductionModelEvaluator::T_RIGHT_IS_RESPONSE);
      leftWall->set_T_left(7.0);
      leftWall->set_T_right(5.0);  // final solution is 6.0
      leftWall->set_k(1.0);
      leftWall->set_q(1.0);
    }

    RCP<LinearHeatConductionModelEvaluator> middleWall;
    {
      middleWall = linearHeatConductionModelEvaluator(globalComm,"middle wall",pike_test::LinearHeatConductionModelEvaluator::T_RIGHT_IS_RESPONSE);
      middleWall->set_T_left(6.0);
      middleWall->set_T_right(3.0);  // final solution is 4.0
      middleWall->set_k(1.0/2.0);
      middleWall->set_q(1.0);
    }

    RCP<LinearHeatConductionModelEvaluator> rightWall;
    {
      rightWall =linearHeatConductionModelEvaluator(globalComm,"right wall",pike_test::LinearHeatConductionModelEvaluator::Q_IS_RESPONSE);
      rightWall->set_T_left(4.0);
      rightWall->set_T_right(1.0);
      rightWall->set_k(1.0/3.0);
      rightWall->set_q(1.5); // final solution is 1.0
    }

    RCP<LinearHeatConductionDataTransfer> transferQ = 
      linearHeatConductionDataTransfer(globalComm,"tranfers q: right->{left,middle}",pike_test::LinearHeatConductionDataTransfer::TRANSFER_Q);
    transferQ->setSource(rightWall);
    transferQ->addTarget(leftWall);
    transferQ->addTarget(middleWall);
    
    RCP<LinearHeatConductionDataTransfer> transferLeftToMiddle =
      linearHeatConductionDataTransfer(globalComm,"tranfer T: left->middle",pike_test::LinearHeatConductionDataTransfer::TRANSFER_T);
    transferLeftToMiddle->setSource(leftWall);
    transferLeftToMiddle->addTarget(middleWall);

    RCP<LinearHeatConductionDataTransfer> transferMiddleToRight =
      linearHeatConductionDataTransfer(globalComm,"tranfer T: middle->right",pike_test::LinearHeatConductionDataTransfer::TRANSFER_T);
    transferMiddleToRight->setSource(middleWall);
    transferMiddleToRight->addTarget(rightWall);

    Teuchos::RCP<pike::Composite> status = pike::composite(pike::Composite::OR);
    Teuchos::RCP<pike::MaxIterations> maxIters =
      Teuchos::rcp(new pike::MaxIterations(20));
    Teuchos::RCP<pike::Composite> convergedTests = 
      pike::composite(pike::Composite::AND);
    Teuchos::RCP<pike::ScalarResponseRelativeTolerance> t1 = 
      Teuchos::rcp(new pike::ScalarResponseRelativeTolerance);
    {
      Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
      p->set("Application Name","left wall");
      p->set("Response Name","T_right");
      p->set("Tolerance",1.0e-5);
      t1->setParameterList(p);
    }
    Teuchos::RCP<pike::ScalarResponseRelativeTolerance> t2 = 
      Teuchos::rcp(new pike::ScalarResponseRelativeTolerance);
    {
      Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
      p->set("Application Name","middle wall");
      p->set("Response Name","T_right");
      p->set("Tolerance",1.0e-5);
      t2->setParameterList(p);
    }
    Teuchos::RCP<pike::ScalarResponseRelativeTolerance> t3 = 
      Teuchos::rcp(new pike::ScalarResponseRelativeTolerance);
    {
      Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
      p->set("Application Name","right wall");
      p->set("Response Name","q");
      p->set("Tolerance",1.0e-5);
      t3->setParameterList(p);
    }
    convergedTests->addTest(t1);
    convergedTests->addTest(t2);
    convergedTests->addTest(t3);
    status->addTest(maxIters);
    status->addTest(convergedTests);

    pike::BlockJacobi solver;
    solver.registerModelEvaluator(leftWall);
    solver.registerModelEvaluator(middleWall);
    solver.registerModelEvaluator(rightWall);
    solver.registerDataTransfer(transferQ);
    solver.registerDataTransfer(transferLeftToMiddle);
    solver.registerDataTransfer(transferMiddleToRight);
    solver.completeRegistration();
    solver.setStatusTests(status);
    solver.solve();

    TEST_EQUALITY(solver.getNumberOfIterations(),18);
    TEST_EQUALITY(solver.getStatus(),pike::CONVERGED);
  }

}
