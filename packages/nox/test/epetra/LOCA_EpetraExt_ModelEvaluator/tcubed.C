//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include <iostream>

#include "LOCA.H"
#include "LOCA_Epetra.H"
#include "NOX_TestCompare.H"

// FEApp is defined in Trilinos/packages/sacado/example/FEApp
#include "FEApp_ModelEvaluator.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

int main(int argc, char *argv[]) {
  unsigned int nelem = 100;
  double h = 1.0/nelem;
  double alpha = 1.0;
  double leftBC = 0.0;
  double rightBC = 0.1;
  int maxNewtonIters = 15;
  int ierr = 0;
  int MyPID;

  try {

    // Initialize MPI
#ifdef HAVE_MPI
    MPI_Init(&argc,&argv);
#endif

    // Create a communicator for Epetra objects
    Teuchos::RefCountPtr<Epetra_Comm> Comm;
#ifdef HAVE_MPI
    Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    Comm = Teuchos::rcp(new Epetra_SerialComm);
#endif

    MyPID = Comm->MyPID();

    // Check for verbose output
    bool verbose = false;
    if (argc>1) 
      if (argv[1][0]=='-' && argv[1][1]=='v') 
        verbose = true;
    
    // Create mesh
    vector<double> x(nelem+1);
    for (unsigned int i=0; i<=nelem; i++)
      x[i] = h*i;

    // Set up application parameters
    Teuchos::RefCountPtr<Teuchos::ParameterList> appParams = 
      Teuchos::rcp(new Teuchos::ParameterList);
    Teuchos::ParameterList& problemParams = 
      appParams->sublist("Problem");
    problemParams.set("Name", "Heat Nonlinear Source");
    problemParams.set("Left BC", leftBC);
    problemParams.set("Right BC", rightBC);
    Teuchos::ParameterList& sourceParams = 
      problemParams.sublist("Source Function");
    sourceParams.set("Name", "Cubic");
    sourceParams.set("Nonlinear Factor", alpha);
    

    // Create set of free parameters for LOCA to use
    Teuchos::RefCountPtr< Teuchos::Array<std::string> > free_param_names =
      Teuchos::rcp(new Teuchos::Array<std::string>);
    free_param_names->push_back("Constant Node BC 1");
    free_param_names->push_back("Constant Node BC 2");
    free_param_names->push_back("Cubic Source Function Nonlinear Factor");

    // Create application
    Teuchos::RefCountPtr<FEApp::Application> app = 
      Teuchos::rcp(new FEApp::Application(x, Comm, appParams, true));

    // Create initial guess
    Teuchos::RefCountPtr<const Epetra_Vector> u = app->getInitialSolution();
    NOX::Epetra::Vector nox_u(*u);

    // Set up LOCA parameters
    Teuchos::RefCountPtr<Teuchos::ParameterList> locaParams =
      Teuchos::rcp(&(appParams->sublist("LOCA")),false);

    // Create the stepper sublist and set the stepper parameters
    Teuchos::ParameterList& stepperList = locaParams->sublist("Stepper");
    stepperList.set("Bordered Solver Method", "Householder");
    stepperList.set("Continuation Parameter", "Constant Node BC 2");
    stepperList.set("Initial Value", rightBC);
    stepperList.set("Max Value", 100.0);
    stepperList.set("Min Value", 0.05);
    stepperList.set("Max Steps", 30);
    stepperList.set("Max Nonlinear Iterations", maxNewtonIters);

#ifdef HAVE_LOCA_ANASAZI
    // Create Anasazi Eigensolver sublist (needs --with-loca-anasazi)
    stepperList.set("Compute Eigenvalues",true);
    Teuchos::ParameterList& aList = stepperList.sublist("Eigensolver");
    aList.set("Method", "Anasazi");
    if (!verbose)
      aList.set("Verbosity", Anasazi::Errors);
#else
    stepperList.set("Compute Eigenvalues",false);
#endif
  
    // Create predictor sublist
    Teuchos::ParameterList& predictorList = locaParams->sublist("Predictor");
    predictorList.set("Method", "Tangent");

    // Create step size sublist
    Teuchos::ParameterList& stepSizeList = locaParams->sublist("Step Size");
    stepSizeList.set("Initial Step Size", 0.1);
    stepSizeList.set("Min Step Size", 1.0e-3);
    stepSizeList.set("Max Step Size", 2000.0);
    stepSizeList.set("Aggressiveness", 0.1);
    
    // Set up NOX parameters
    Teuchos::RefCountPtr<Teuchos::ParameterList> noxParams =
      Teuchos::rcp(&(appParams->sublist("NOX")),false);

    // Set the nonlinear solver method
    noxParams->set("Nonlinear Solver", "Line Search Based");

    // Set the printing parameters in the "Printing" sublist
    Teuchos::ParameterList& printParams = noxParams->sublist("Printing");
    printParams.set("MyPID", MyPID); 
    printParams.set("Output Precision", 3);
    printParams.set("Output Processor", 0);
    if (verbose)
      printParams.set("Output Information", 
                      NOX::Utils::OuterIteration + 
                      NOX::Utils::OuterIterationStatusTest + 
                      NOX::Utils::InnerIteration +
//                       NOX::Utils::Details + 
//                       NOX::Utils::LinearSolverDetails +
                      NOX::Utils::StepperIteration + 
                      NOX::Utils::StepperDetails +
                      NOX::Utils::StepperParameters + 
                      NOX::Utils::TestDetails + 
                      NOX::Utils::Warning + 
                      NOX::Utils::Error);
    else
      printParams.set("Output Information", NOX::Utils::Error);

    // Create printing utilities
    NOX::Utils utils(printParams);

    // Sublist for line search 
    Teuchos::ParameterList& searchParams = noxParams->sublist("Line Search");
    searchParams.set("Method", "Full Step");

    // Sublist for direction
    Teuchos::ParameterList& dirParams = noxParams->sublist("Direction");
    dirParams.set("Method", "Newton");
    Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
    newtonParams.set("Forcing Term Method", "Constant");

    // Sublist for linear solver for the Newton method
    Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
    lsParams.set("Aztec Solver", "GMRES");  
    lsParams.set("Max Iterations", 100);  
    lsParams.set("Tolerance", 1e-4); 
    if (verbose)
      lsParams.set("Output Frequency", 1);
    else
      lsParams.set("Output Frequency", 0);
    lsParams.set("Preconditioner", "Ifpack");

    // Create model evaluator
    Teuchos::RefCountPtr<FEApp::ModelEvaluator> model = 
      Teuchos::rcp(new FEApp::ModelEvaluator(app, free_param_names));

    // Create Epetra factory
    Teuchos::RefCountPtr<LOCA::Abstract::Factory> epetraFactory =
      Teuchos::rcp(new LOCA::Epetra::Factory);

    // Create global data object
    Teuchos::RefCountPtr<LOCA::GlobalData> globalData = 
      LOCA::createGlobalData(appParams, epetraFactory);

    // Create LOCA interface
    Teuchos::RefCountPtr<LOCA::Epetra::ModelEvaluatorInterface> interface =
      Teuchos::rcp(new LOCA::Epetra::ModelEvaluatorInterface(globalData,
                                                             model));

    // Get LOCA parameter vector
    LOCA::ParameterVector pVector = interface->getLOCAParameterVector();

    // Create the Jacobian matrix
    Teuchos::RefCountPtr<Epetra_Operator> A = model->create_W(); 

    // Create the linear system
    Teuchos::RefCountPtr<LOCA::Epetra::Interface::Required> iReq = interface;
    Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian> iJac = interface;
    Teuchos::RefCountPtr<NOX::Epetra::LinearSystemAztecOO> linsys = 
      Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, 
                                                        lsParams,
                                                        iReq, iJac, A, nox_u));

    // Create the Group
    Teuchos::RefCountPtr<LOCA::Epetra::Group> grp =
      Teuchos::rcp(new LOCA::Epetra::Group(globalData, printParams, iReq, 
                                           nox_u, linsys, pVector)); 
    grp->setDerivUtils(interface);

    // Create the Solver convergence test
    Teuchos::RefCountPtr<NOX::StatusTest::NormF> wrms = 
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
    Teuchos::RefCountPtr<NOX::StatusTest::MaxIters> maxiters = 
      Teuchos::rcp(new NOX::StatusTest::MaxIters(maxNewtonIters));
    Teuchos::RefCountPtr<NOX::StatusTest::Combo> combo = 
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    combo->addStatusTest(wrms);
    combo->addStatusTest(maxiters);

    // Create the stepper  
    LOCA::Stepper stepper(globalData, grp, combo, appParams);

    // Run the stepper
    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

    if (status != LOCA::Abstract::Iterator::Finished) {
      if (globalData->locaUtils->isPrintType(NOX::Utils::Error))
        globalData->locaUtils->out() 
          << "Stepper failed to converge!" << std::endl;
    }

    // Get the final solution from the stepper
    Teuchos::RCP<const LOCA::Epetra::Group> finalGroup = 
      Teuchos::rcp_dynamic_cast<const LOCA::Epetra::Group>(stepper.getSolutionGroup());
    const NOX::Epetra::Vector& finalSolution = 
      dynamic_cast<const NOX::Epetra::Vector&>(finalGroup->getX());

    // Output the parameter list
    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperParameters)) {
      globalData->locaUtils->out() 
        << std::endl << "Final Parameters" << std::endl
        << "****************" << std::endl;
      stepper.getList()->print(globalData->locaUtils->out());
      globalData->locaUtils->out() << std::endl;
    }

    // Check some statistics on the solution
    NOX::TestCompare testCompare(globalData->locaUtils->out(), 
				 *(globalData->locaUtils));
  
    if (globalData->locaUtils->isPrintType(NOX::Utils::TestDetails))
      globalData->locaUtils->out() 
        << std::endl 
        << "***** Checking solution statistics *****" 
        << std::endl;

    // Check number of steps
    int numSteps = stepper.getStepNumber();
    int numSteps_expected = 23;
    ierr += testCompare.testValue(numSteps, numSteps_expected, 0.0,
				  "number of continuation steps",
				  NOX::TestCompare::Absolute);

    // Check number of failed steps
    int numFailedSteps = stepper.getNumFailedSteps();
    int numFailedSteps_expected = 0;
    ierr += testCompare.testValue(numFailedSteps, numFailedSteps_expected, 0.0,
                                  "number of failed continuation steps",
                                  NOX::TestCompare::Absolute);

    // Check final value of continuation parameter
    double right_bc_final = finalGroup->getParam("Constant Node BC 2");
    double right_bc_expected = 0.05;
    ierr += testCompare.testValue(right_bc_final, right_bc_expected, 1.0e-14,
                                  "final value of continuation parameter", 
                                  NOX::TestCompare::Relative);
 
    // Check norm of solution
    double norm_x = finalSolution.norm();
    double norm_x_expected = 25.00498021;
    ierr += testCompare.testValue(norm_x, norm_x_expected, 1.0e-7,
                                  "norm of final solution",
                                  NOX::TestCompare::Relative);

    LOCA::destroyGlobalData(globalData);

  }
  
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
    ierr = 1;
  }
  catch (string& s) {
    std::cout << s << std::endl;
    ierr = 1;
  }
  catch (const char *s) {
    std::cout << s << std::endl;
    ierr = 1;
  }
  catch (...) {
    std::cout << "Caught unknown exception!" <<std:: endl;
    ierr = 1;
  }

  if (MyPID == 0) {
    if (ierr == 0)
      std::cout << "All tests passed!" << std::endl;
    else
      std::cout << ierr << " test(s) failed!" << std::endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return ierr;
}
