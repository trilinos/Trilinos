#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Pike_BlackBox_config.hpp"
#include "Pike_MultiphysicsDistributor.hpp"
#include "Pike_StatusTest_Factory.hpp"
#include "Pike_SolverObserver_Logger.hpp"

// Solvers
#include "Pike_Solver_BlockGaussSeidel.hpp"
#include "Pike_Solver_BlockJacobi.hpp"
#include "Pike_Solver_Factory.hpp"
#include "Pike_Mock_UserSolverFactory.hpp"

// Models
#include "Pike_VanderPol_ModelEvaluator_Eq1.hpp"
#include "Pike_VanderPol_ModelEvaluator_Eq2.hpp"
#include "Pike_VanderPol_DataTransfer_Eq1ToEq2.hpp"
#include "Pike_VanderPol_DataTransfer_Eq2ToEq1.hpp"

#include <iostream>

namespace pike_test {

  TEUCHOS_UNIT_TEST(TransientStepper, VanDerPol)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::MpiComm;
    using Teuchos::ParameterList;
    
    RCP<MpiComm<int> > globalComm = 
      rcp(new MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));

    // Run this on 2 processes only
    TEST_EQUALITY(globalComm->getSize(), 2);

    RCP<pike::MultiphysicsDistributor> dist = rcp(new pike::MultiphysicsDistributor);
    typedef pike::MultiphysicsDistributor::ApplicationIndex AppIndex;
    AppIndex EQ1 = dist->addApplication("Eq1",0,0);
    AppIndex EQ2 = dist->addApplication("Eq2",1,1);
    dist->addTransfer("Eq1->Eq2",EQ1,EQ2);
    dist->addTransfer("Eq2->Eq1",EQ2,EQ1);
    dist->setup(globalComm,true);

    RCP<VanderPolME1> me1 = vanderPolME1(dist);
    RCP<VanderPolME2> me2 = vanderPolME2(dist);
    RCP<VanderPolDT1To2> dt1To2 = vanderPolDT1To2(dist,me1,me2);
    RCP<VanderPolDT2To1> dt2To1 = vanderPolDT2To1(dist,me1,me2);

    // Special hooks to insert failures to test logic paths in
    // transient solver
    me1->setUnitTestFailureForTimeStep(3);

    RCP<ParameterList> p = Teuchos::parameterList("Transient Solver");
    {
      p->set("Solver Sublist Name", "My Transient Test");

      ParameterList& pt = p->sublist("My Transient Test");
      pt.set("Type","Transient Stepper");
      pt.set("Maximum Number of Time Steps",100);
      pt.set("Begin Time", 0.0);
      pt.set("End Time", 1.0);
      pt.set("Initial Time Step Size",1.0e-3);
      pt.set("Minimum Time Step Size",1.0e-5);
      pt.set("Maximum Time Step Size",10.0);
      Teuchos::Array<double> checkPoints;
      checkPoints.push_back(0.51);
      pt.set("Check Points",checkPoints);      
      pt.set("Internal Solver Sublist","My Jacobi Solver");

      ParameterList& pj = p->sublist("My Jacobi Solver");
      pj.set("Type","Block Gauss Seidel");
      pj.set("MPI Barrier Transfers",true);
      pj.set("Name","My Special Solver");

      p->print(*dist->getSerialOStream(),0,true);
    }

    RCP<pike::StatusTest> tests;
    {
      RCP<ParameterList> stp = Teuchos::parameterList("Status Test Builder");
      stp->set("Type","Composite OR");
      Teuchos::ParameterList& failure = stp->sublist("Failure");
      failure.set("Type","Maximum Iterations");
      failure.set("Maximum Iterations",6);

      // For unit testing to trigger failure of a complete time step
      Teuchos::ParameterList& app1LocalFail = stp->sublist("Local Failure 1");      
      app1LocalFail.set("Type","Local Model Failure");
      app1LocalFail.set("Model Name","Eq1");

      Teuchos::ParameterList& converged = stp->sublist("Converged");
      converged.set("Type","Composite AND");
      Teuchos::ParameterList& relTolApp1 = converged.sublist("Van der Pol Equation 1");
      relTolApp1.set("Type","Scalar Response Relative Tolerance");
      relTolApp1.set("Application Name","Eq1");
      relTolApp1.set("Response Name","x1");
      relTolApp1.set("Tolerance",1.0e-3);
      Teuchos::ParameterList& relTolApp2 = converged.sublist("Van der Pol Equation 2"); 
      relTolApp2.set("Type","Scalar Response Relative Tolerance");
      relTolApp2.set("Application Name","Eq2");
      relTolApp2.set("Response Name","x2");
      relTolApp2.set("Tolerance",1.0e-3);
      Teuchos::ParameterList& app1LocalConv = converged.sublist("Local Convergence 1");      
      app1LocalConv.set("Type","Local Model Convergence");
      app1LocalConv.set("Model Name","Eq1");
      Teuchos::ParameterList& app2LocalConv = converged.sublist("Local Convergence 2");      
      app2LocalConv.set("Type","Local Model Convergence");
      app2LocalConv.set("Model Name","Eq2");

      pike::StatusTestFactory stFactory;
      tests = stFactory.buildStatusTests(stp);
    }

    pike::SolverFactory factory;
    RCP<pike::Solver> solver = factory.buildSolver(p);
    solver->registerComm(globalComm);
    solver->registerModelEvaluator(me1);
    solver->registerModelEvaluator(me2);
    solver->registerDataTransfer(dt1To2);
    solver->registerDataTransfer(dt2To1);
    solver->completeRegistration();
    solver->setStatusTests(tests);
    solver->addObserver(pike::loggerObserver());
    solver->initialize();
    solver->solve();
    solver->finalize();
    
    // Takes 15 time steps to converge.  Should be 12 if we don't
    // introduce artificial time step failure.
    TEST_EQUALITY(solver->getNumberOfIterations(),15);
    TEST_EQUALITY(solver->getStatus(),pike::CONVERGED);

    // Test accessors for code coverage
    TEST_ASSERT(nonnull(solver->getModelEvaluator("Eq1")));
    TEST_EQUALITY(solver->getModelEvaluators().size(),2)
    TEST_ASSERT(nonnull(solver->getDataTransfer("Eq1->Eq2")));
    TEST_EQUALITY(solver->getDataTransfers().size(),2);
    TEST_EQUALITY(solver->name(), "My Special Solver");
    RCP<pike::SolverObserver> logger = (solver->getObservers())[0];
    TEST_EQUALITY((Teuchos::rcp_dynamic_cast<pike::LoggerObserver>(logger))->getLog()->size(),82);
    TEST_THROW(solver->reset(), std::logic_error);
  }

}
