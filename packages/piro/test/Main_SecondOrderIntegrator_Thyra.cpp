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

#include "MockModelEval_B_Tpetra_legacy.hpp"

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "Piro_ConfigDefs.hpp"
#include "Piro_Test_MockObserver.hpp"
#include "Piro_SolverFactory.hpp"
#include "Piro_StratimikosUtils.hpp"


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

  for (int iTest=0; iTest<2; iTest++) {

    if (doAll) {
      switch (iTest) {
       case 1: inputFile="input_Solve_VV_mod.xml"; break;
       case 0: inputFile="input_Solve_TR_mod.xml"; break;
       //case 1: inputFile="input_Solve_NB_mod.xml"; break;
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
      auto comm = Tpetra::getDefaultComm();
      auto model = rcp(new MockModelEval_B_Tpetra_legacy(comm));

      RCP<Teuchos::ParameterList> piroParams =
         rcp(new Teuchos::ParameterList("Piro Parameters"));
      Teuchos::updateParametersFromXmlFile(inputFile, piroParams.ptr());
      auto stratParams = Piro::extractStratimikosParams(piroParams);
      std::cout << "Piro Params" <<std::endl;
      piroParams->print();
      Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
      std::cout << "\n\n\nStrat Params" <<std::endl;
      stratParams->print();
      std::cout << "\n\n\n" <<std::endl;
      linearSolverBuilder.setParameterList(stratParams);
      auto lowsFactory = createLinearSolveStrategy(linearSolverBuilder);
      auto modelWithSolve = rcp(new Thyra::DefaultModelEvaluatorWithSolveFactory<double>(model, lowsFactory));

      auto observer = rcp(new Piro::Test::MockObserver<double>());

      Piro::SolverFactory solverFactory;
      auto solver =  solverFactory.createSolver<double>(piroParams, modelWithSolve, Teuchos::null, observer);

      // Now the (somewhat cumbersome) setting of inputs and outputs

      auto inArgs = solver->getNominalValues();

      auto outArgs = solver->createOutArgs();
      TEUCHOS_ASSERT(solver->Ng() == 2); // Number of *vectors* of responses

      const RCP<Thyra::VectorBase<double> > g0 = Thyra::createMember(solver->get_g_space(0));
      outArgs.set_g(0, g0);

      const int solutionResponseIndex = solver->Ng() - 1;  //1 in this case
      const RCP<Thyra::VectorBase<double> > solution = Thyra::createMember(solver->get_g_space(solutionResponseIndex));
      outArgs.set_g(solutionResponseIndex, solution);

      // Now, solve the problem and return the responses
      solver->evalModel(inArgs, outArgs);

      // Print out everything
      if (Proc == 0)
        std::cout << "Finished Model Evaluation: Printing everything {Exact in brackets}" 
             << "\n-----------------------------------------------------------------"
             << std::setprecision(9) << std::endl;

      //outArgs.get_p(0)->Print(std::cout << "\nParameters! {1}\n");
      //g0->Print(std::cout << "\nResponses! {0.0}\n");
     // solution->Print(std::cout << "\nSolution! {0.0}\n");

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
