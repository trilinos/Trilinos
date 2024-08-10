// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>
#include <string>

#include "MockModelEval_B.hpp"
#include "ObserveSolution_Epetra.hpp"

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "Piro_ConfigDefs.hpp"

#include "Piro_Epetra_VelocityVerletSolver.hpp"
#include "Piro_Epetra_NOXSolver.hpp"
#include "Piro_Epetra_TrapezoidRuleSolver.hpp"
#include "Piro_Epetra_NewmarkSolver.hpp"


int main(int argc, char *argv[]) {

  int status=0; // 0 = pass, failures are incremented
  int overall_status=0; // 0 = pass, failures are incremented over multiple tests
  bool success=true;

  // Initialize MPI 
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  int Proc=mpiSession.getRank();
#ifdef HAVE_MPI
  MPI_Comm appComm = MPI_COMM_WORLD;
#else
  int appComm=0;
#endif

  using Teuchos::RCP;
  using Teuchos::rcp;
  std::string inputFile;

  bool doAll = (argc==1);
  if (argc>1) doAll = !strcmp(argv[1],"-v");

  for (int iTest=0; iTest<3; iTest++) {

    if (doAll) {
      switch (iTest) {
       case 0: inputFile="input_Solve_VV.xml"; break;
       case 1: inputFile="input_Solve_TR.xml"; break;
       case 2: inputFile="input_Solve_NB.xml"; break;
       default : std::cout << "iTest logic error " << std::endl; exit(-1);
      }
    }
     else {
      inputFile=argv[1];
      iTest = 999;
    }

    if (Proc==0)
     std::cout << "===================================================\n"
          << "======  Running input file: "<< inputFile <<"\n"
          << "===================================================\n"
          << std::endl;
    
    try {

      // Create (1) a Model Evaluator and (2) a ParameterList
      RCP<EpetraExt::ModelEvaluator> Model = rcp(new MockModelEval_B(appComm));

      RCP<Teuchos::ParameterList> piroParams =
         rcp(new Teuchos::ParameterList("Piro Parameters"));
      Teuchos::updateParametersFromXmlFile(inputFile, piroParams.ptr());

      // Use these two objects to construct a Piro solved application 
      //   EpetraExt::ModelEvaluator is  base class of all Piro::Epetra solvers
      RCP<EpetraExt::ModelEvaluator> piro;

      std::string& solver = piroParams->get("Solver Type","NOX");
      RCP<NOX::Epetra::Observer> observer = rcp(new ObserveSolution_Epetra());

      if (solver=="NOX") {
        piro = rcp(new Piro::Epetra::NOXSolver(piroParams, Model, observer));
      }
      else if (solver=="Velocity Verlet") {
        piro = rcp(new Piro::Epetra::VelocityVerletSolver(
                       piroParams, Model, observer));
      }
      else if (solver=="Trapezoid Rule") {
        piro = rcp(new Piro::Epetra::TrapezoidRuleSolver(
                       piroParams, Model, observer));
      }
      else if (solver=="Newmark") {
        piro = rcp(new Piro::Epetra::NewmarkSolver(
                       piroParams, Model, observer));
      }
      else
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Error: Unknown Piro Solver : " << solver);

      // Now the (somewhat cumbersome) setting of inputs and outputs
      EpetraExt::ModelEvaluator::InArgs inArgs = piro->createInArgs();
      TEUCHOS_ASSERT(inArgs.Np() > 0); // Number of *vectors* of parameters
      RCP<Epetra_Vector> p1 = rcp(new Epetra_Vector(*(piro->get_p_init(0))));
      inArgs.set_p(0,p1);

      // Set output arguments to evalModel call
      EpetraExt::ModelEvaluator::OutArgs outArgs = piro->createOutArgs();
      TEUCHOS_ASSERT(outArgs.Ng() >= 2); // Number of *vectors* of responses
      RCP<Epetra_Vector> g1 = rcp(new Epetra_Vector(*(piro->get_g_map(0))));
      outArgs.set_g(0,g1);
      // Solution vector is returned as extra respons vector
      RCP<Epetra_Vector> gx = rcp(new Epetra_Vector(*(piro->get_g_map(1))));
      outArgs.set_g(1,gx);

      // Now, solve the problem and return the responses
      piro->evalModel(inArgs, outArgs);

      // Print out everything
      if (Proc == 0)
        std::cout << "Finished Model Evaluation: Printing everything {Exact in brackets}" 
             << "\n-----------------------------------------------------------------"
             << std::setprecision(9) << std::endl;

      p1->Print(std::cout << "\nParameters! {1}\n");
      g1->Print(std::cout << "\nResponses! {0.0}\n");
      gx->Print(std::cout << "\nSolution! {0.0}\n");

      if (Proc == 0)
        std::cout <<
          "\n-----------------------------------------------------------------\n";
    }
    TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
    if (!success) status=10; else status=0;

    overall_status += status;

  }

  if (Proc==0) {
    if (overall_status==0) std::cout << "\nTEST PASSED\n" << std::endl;
    else std::cout << "\nTEST Failed:  " << overall_status << std::endl;
  }

  return status;
}
