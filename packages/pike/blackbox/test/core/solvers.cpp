#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Pike_BlackBox_config.hpp"

// Solvers
#include "Pike_Solver_BlockGaussSeidel.hpp"
#include "Pike_Solver_BlockJacobi.hpp"
#include "Pike_Solver_Factory.hpp"
#include "Pike_Mock_UserSolverFactory.hpp"

// Models
#include "Pike_LinearHeatConduction_ModelEvaluator.hpp"
#include "Pike_LinearHeatConduction_DataTransfer.hpp"
#include "Pike_BlackBoxModelEvaluator_SolverAdapter.hpp"

// Status tests
#include "Pike_StatusTest_Composite.hpp"
#include "Pike_StatusTest_MaxIterations.hpp"
#include "Pike_StatusTest_ScalarResponseRelativeTolerance.hpp"

#include <iostream>

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
    Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
    p->set("Type","Block Jacobi");
    p->set("MPI Barrier Transfers",true);
    p->set("MPI Barrier Solves",true);
    solver.setParameterList(p);
    solver.registerModelEvaluator(leftWall);
    solver.registerModelEvaluator(middleWall);
    solver.registerModelEvaluator(rightWall);
    solver.registerDataTransfer(transferQ);
    solver.registerDataTransfer(transferLeftToMiddle);
    solver.registerDataTransfer(transferMiddleToRight);
    solver.registerComm(globalComm);
    solver.completeRegistration();
    solver.setStatusTests(status);
    solver.solve();

    TEST_EQUALITY(solver.getNumberOfIterations(),18);
    TEST_EQUALITY(solver.getStatus(),pike::CONVERGED);
  }

  TEUCHOS_UNIT_TEST(solvers, factory)
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

    Teuchos::RCP<Teuchos::ParameterList> solverParams = Teuchos::parameterList();
    Teuchos::updateParametersFromXmlFileAndBroadcast("solver_factory_test_params.xml",
						     solverParams.ptr(),
						     *globalComm);

    pike::SolverFactory factory;
    Teuchos::RCP<pike_test::UserSolverFactory> userFactory1 = 
      Teuchos::rcp(new pike_test::UserSolverFactory("My Super-Special Solver"));
    factory.addFactory(userFactory1);
    Teuchos::RCP<pike_test::UserSolverFactory> userFactory2 = 
      Teuchos::rcp(new pike_test::UserSolverFactory("My Other Super-Special Solver"));
    factory.addFactory(userFactory2);
    
    Teuchos::RCP<pike::Solver> solver = factory.buildSolver(solverParams);
    solver->registerModelEvaluator(leftWall);
    solver->registerModelEvaluator(middleWall);
    solver->registerModelEvaluator(rightWall);
    solver->registerDataTransfer(transferQ);
    solver->registerDataTransfer(transferLeftToMiddle);
    solver->registerDataTransfer(transferMiddleToRight);
    solver->completeRegistration();
    solver->setStatusTests(status);
    solver->solve();

    TEST_EQUALITY(solver->getNumberOfIterations(),18);
    TEST_EQUALITY(solver->getStatus(),pike::CONVERGED);

    // for coverage testing
    TEST_ASSERT(factory.supportsType("Block Jacobi"));
    TEST_ASSERT(factory.supportsType("My Super-Special Solver"));
    TEST_ASSERT(factory.supportsType("My Other Super-Special Solver"));
    TEST_ASSERT(!factory.supportsType("My Vaporware Solver"));
    { // build a block jacobi solver
      Teuchos::RCP<Teuchos::ParameterList> jpl = Teuchos::parameterList("Test Block Jacobi");
      jpl->set("Solver Sublist Name","My Block Jacobi Solver");
      jpl->sublist("My Block Jacobi Solver").set("Type","Block Jacobi");
      Teuchos::RCP<pike::Solver> js = factory.buildSolver(jpl);
      TEST_ASSERT(nonnull(js));
    }
  }

  TEUCHOS_UNIT_TEST(solvers, hierarchic)
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

    // Inner solve is middle and right wall.  Outer solve is left and
    // inner solve object.  Currently there is no view of the global
    // connectivity graph.  constructing this complicates the library
    // by forcing a hierarchical mix-in interface.  We may want to add
    // this in the future just to autmate things for users.  But for
    // now it is not needed.

    /* 
       Model Connectivity Graph:
      
             S_outer
             /     \
          Left    S_inner
                  /     \
               Middle   Right
    
       Model Transfer Graph (4 transfers)
      
          Left  <->  S_inner
                     /     \
                Middle <-> Right
    */

    RCP<LinearHeatConductionDataTransfer> transferQRightToLeft = 
      linearHeatConductionDataTransfer(globalComm,"tranfers q: right->left",pike_test::LinearHeatConductionDataTransfer::TRANSFER_Q);
    transferQRightToLeft->setSource(rightWall);
    transferQRightToLeft->addTarget(leftWall);

    RCP<LinearHeatConductionDataTransfer> transferQRightToMiddle = 
      linearHeatConductionDataTransfer(globalComm,"tranfers q: right->middle",pike_test::LinearHeatConductionDataTransfer::TRANSFER_Q);
    transferQRightToMiddle->setSource(rightWall);
    transferQRightToMiddle->addTarget(middleWall);
    
    RCP<LinearHeatConductionDataTransfer> transferTLeftToMiddle =
      linearHeatConductionDataTransfer(globalComm,"tranfer T: left->middle",pike_test::LinearHeatConductionDataTransfer::TRANSFER_T);
    transferTLeftToMiddle->setSource(leftWall);
    transferTLeftToMiddle->addTarget(middleWall, "Inner Solver");

    RCP<LinearHeatConductionDataTransfer> transferTMiddleToRight =
      linearHeatConductionDataTransfer(globalComm,"tranfer T: middle->right",pike_test::LinearHeatConductionDataTransfer::TRANSFER_T);
    transferTMiddleToRight->setSource(middleWall);
    transferTMiddleToRight->addTarget(rightWall);

    Teuchos::RCP<pike::BlockGaussSeidel> innerSolver = 
      Teuchos::rcp(new pike::BlockGaussSeidel);
    Teuchos::RCP<Teuchos::ParameterList> innerSolverPL = 
      Teuchos::parameterList();
    innerSolverPL->set("Name","Inner Solver");
    innerSolver->setParameterList(innerSolverPL);
    Teuchos::RCP<Teuchos::FancyOStream> innerOS =
      Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcpFromRef(std::cout)));
    innerOS = Teuchos::tab(innerOS,9);
    innerSolver->setOStream(innerOS);
    {
      Teuchos::RCP<pike::Composite> status = pike::composite(pike::Composite::OR);
      Teuchos::RCP<pike::MaxIterations> maxIters =Teuchos::rcp(new pike::MaxIterations(20));
      Teuchos::RCP<pike::Composite> convergedTests = pike::composite(pike::Composite::AND);

      Teuchos::RCP<pike::ScalarResponseRelativeTolerance> t1 = 
	Teuchos::rcp(new pike::ScalarResponseRelativeTolerance);
      {
	Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
	p->set("Application Name","middle wall");
	p->set("Response Name","T_right");
	p->set("Tolerance",1.0e-5);
	t1->setParameterList(p);
      }
      Teuchos::RCP<pike::ScalarResponseRelativeTolerance> t2 = 
	Teuchos::rcp(new pike::ScalarResponseRelativeTolerance);
      {
	Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
	p->set("Application Name","right wall");
	p->set("Response Name","q");
	p->set("Tolerance",1.0e-5);
	t2->setParameterList(p);
      }
      convergedTests->addTest(t1);
      convergedTests->addTest(t2);
      status->addTest(maxIters);
      status->addTest(convergedTests);

      innerSolver->registerModelEvaluator(middleWall);
      innerSolver->registerModelEvaluator(rightWall);
      innerSolver->registerDataTransfer(transferQRightToMiddle);
      innerSolver->registerDataTransfer(transferTMiddleToRight);
      innerSolver->completeRegistration();
      innerSolver->setStatusTests(status);
    }

    Teuchos::RCP<pike::BlockGaussSeidel> outerSolver = 
      Teuchos::rcp(new pike::BlockGaussSeidel);
    Teuchos::RCP<Teuchos::ParameterList> outerSolverPL = 
      Teuchos::parameterList();
    outerSolverPL->set("Name","Outer Solver");
    outerSolver->setParameterList(outerSolverPL);
    Teuchos::RCP<Teuchos::FancyOStream> outerOS = 
      Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcpFromRef(std::cout)));
    outerSolver->setOStream(outerOS);
    {
      Teuchos::RCP<pike::Composite> status = pike::composite(pike::Composite::OR);
      Teuchos::RCP<pike::MaxIterations> maxIters =Teuchos::rcp(new pike::MaxIterations(20));
      Teuchos::RCP<pike::Composite> convergedTests = pike::composite(pike::Composite::AND);

      Teuchos::RCP<pike::ScalarResponseRelativeTolerance> t0 = 
	Teuchos::rcp(new pike::ScalarResponseRelativeTolerance);
      {
	Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
	p->set("Application Name","left wall");
	p->set("Response Name","T_right");
	p->set("Tolerance",1.0e-5);
	t0->setParameterList(p);
      }
      Teuchos::RCP<pike::ScalarResponseRelativeTolerance> t1 = 
	Teuchos::rcp(new pike::ScalarResponseRelativeTolerance);
      {
	Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
	p->set("Application Name","Inner Solver");
	p->set("Response Name","T_right");
	p->set("Tolerance",1.0e-5);
	t1->setParameterList(p);
      }
      Teuchos::RCP<pike::ScalarResponseRelativeTolerance> t2 = 
	Teuchos::rcp(new pike::ScalarResponseRelativeTolerance);
      {
	Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
	p->set("Application Name","Inner Solver");
	p->set("Response Name","q");
	p->set("Tolerance",1.0e-5);
	t2->setParameterList(p);
      }
      convergedTests->addTest(t0);
      convergedTests->addTest(t1);
      convergedTests->addTest(t2);
      status->addTest(maxIters);
      status->addTest(convergedTests);

      Teuchos::RCP<pike::SolverAdapterModelEvaluator> solver2ModelAdapter = 
	Teuchos::rcp(new pike::SolverAdapterModelEvaluator("Inner Solver"));
      solver2ModelAdapter->setSolver(innerSolver);

      outerSolver->registerModelEvaluator(leftWall);
      outerSolver->registerModelEvaluator(solver2ModelAdapter);
      outerSolver->registerDataTransfer(transferQRightToLeft);
      outerSolver->registerDataTransfer(transferTLeftToMiddle);
      outerSolver->completeRegistration();
      outerSolver->setStatusTests(status);
    }

    outerSolver->solve();

    TEST_EQUALITY(outerSolver->getNumberOfIterations(),6);
    TEST_EQUALITY(outerSolver->getStatus(),pike::CONVERGED);

    // Coverage testing extra checks
    TEST_ASSERT(nonnull(Teuchos::rcp_dynamic_cast<const pike::SolverAdapterModelEvaluator>(outerSolver->getModelEvaluator("Inner Solver"))->getSolver()));
    TEST_ASSERT(nonnull(Teuchos::rcp_dynamic_cast<const pike::SolverAdapterModelEvaluator>(outerSolver->getModelEvaluator("Inner Solver"))->getNonconstSolver()));
    Teuchos::RCP<const pike::BlackBoxModelEvaluator> testInnerSolver =  outerSolver->getModelEvaluator("Inner Solver");
    TEST_ASSERT(testInnerSolver->isLocallyConverged());
    TEST_ASSERT(testInnerSolver->isGloballyConverged());
    TEST_EQUALITY(testInnerSolver->getNumberOfParameters(),2);
    TEST_ASSERT(testInnerSolver->supportsParameter("q"));
    TEST_EQUALITY(testInnerSolver->getParameterIndex("q"),0);
    TEST_EQUALITY(testInnerSolver->getParameterName(0),"q");
    Teuchos::Array<double> val(1);
    val[0] = 1.0;
    Teuchos::rcp_const_cast<pike::BlackBoxModelEvaluator>(testInnerSolver)->setParameter(0,val);
    TEST_EQUALITY(testInnerSolver->getNumberOfResponses(),2);
    TEST_ASSERT(testInnerSolver->supportsResponse("T_right"));
    TEST_EQUALITY(testInnerSolver->getResponseName(0),"T_right");
    TEST_EQUALITY(testInnerSolver->getResponseIndex("T_right"),0);
    double tol = 1.0e-2;
    TEST_FLOATING_EQUALITY(testInnerSolver->getResponse(0)[0],3.24,tol);
  }

}
