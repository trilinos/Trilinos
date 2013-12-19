#include "Teuchos_UnitTestHarness.hpp"
#include "Pike_BlackBox_config.hpp"
#include "Pike_Solver_BlockGaussSeidel.hpp"
#include "Pike_LinearHeatConductionModelEvaluator.hpp"
#include "Pike_LinearHeatConduction_DataTransfer.hpp"

#include "Teuchos_DefaultMpiComm.hpp"

namespace pike_test {

  TEUCHOS_UNIT_TEST(app, LinearHeatConduction_BlockJacobi)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    RCP<Teuchos::MpiComm<int> > globalComm = rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));

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

    bool converged = false;
    bool failed = false;
    double tol = 1.0e-5;
    int maxIters = 20;
    int iter = 0;

    std::vector<double> oldValues(3);

    pike::BlockGaussSeidel solver;

    while (!converged && !failed) {

      oldValues[0] = leftWall->get_T_right();
      oldValues[1] = middleWall->get_T_right();
      oldValues[2] = rightWall->get_q();

      transferQ->doTransfer(solver);
      leftWall->solve();
      std::cout << "\nq = " << rightWall->get_q() << std::endl;
      std::cout << "leftWall->T_right = " << leftWall->get_T_right() << std::endl;
      transferLeftToMiddle->doTransfer(solver);
      std::cout << "middleWall->T_left = " << middleWall->get_T_left() << std::endl;
      middleWall->solve();
      transferMiddleToRight->doTransfer(solver);
      rightWall->solve();
      
      ++iter;
      
      if ( (fabs(leftWall->get_T_right() - oldValues[0]) < tol) &&
	   (fabs(middleWall->get_T_right() - oldValues[1]) < tol) &&
	   (fabs(rightWall->get_q() - oldValues[2]) < tol)  )
	converged = true;

      if (iter >= maxIters)
	failed = true;
      
      std::cout << "iter " << iter 
		<< "\n  leftWall(" << leftWall->get_T_right() << "," << oldValues[0] << "," << fabs(leftWall->get_T_right() - oldValues[0]) << ")"
		<< "\n  middleWall(" << middleWall->get_T_right() << "," << oldValues[1] << "," << fabs(middleWall->get_T_right() - oldValues[1]) << ")" 		
		<< "\n  rightWall(" << rightWall->get_q() << "," << oldValues[2] << "," << fabs(rightWall->get_q() - oldValues[2]) << ")\n";
    }

    TEUCHOS_ASSERT(converged);

  }

}
