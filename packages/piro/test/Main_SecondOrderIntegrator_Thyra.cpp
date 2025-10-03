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

#include "MockModelEval_C_Tpetra.hpp"

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Thyra_DetachedVectorView.hpp"

#include "Piro_ConfigDefs.hpp"
#include "Piro_PerformSolve.hpp"
#include "Piro_Test_MockObserver.hpp"
#include "Piro_SolverFactory.hpp"
#include "Piro_StratimikosUtils.hpp"


int main(int argc, char *argv[]) {

  int status=0; // 0 = pass, failures are incremented
  int overall_status=0; // 0 = pass, failures are incremented over multiple tests
  bool success=true;

  // Initialize MPI 
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  int Proc = mpiSession.getRank();

  Teuchos::oblackholestream blackHole;
  std::ostream &out = (Proc == 0) ? std::cout : blackHole;

  auto appComm = Tpetra::getDefaultComm();

  using Teuchos::RCP;
  using Teuchos::rcp;
  std::string inputFile;

  bool doAll = (argc==1);
  if (argc>1) doAll = !strcmp(argv[1],"-v");

  double tol = 1e-12;

  for (int iTest=0; iTest<1; iTest++) {

    if (doAll) {
      switch (iTest) {
       case 0: inputFile="input_Solve_TR_Thyra.xml"; break;
       default : out << "iTest logic error " << std::endl; exit(-1);
      }
    }
     else {
      inputFile=argv[1];
      iTest = 999;
    }

    out << "===================================================\n"
        << "======  Running input file: "<< inputFile <<"\n"
        << "===================================================\n"
        << std::endl;
    
    try {

      // Create (1) a Model Evaluator and (2) a ParameterList
      auto comm = Tpetra::getDefaultComm();
      auto model = rcp(new MockModelEval_C_Tpetra(comm));

      RCP<Teuchos::ParameterList> piroParams =
         rcp(new Teuchos::ParameterList("Piro Parameters"));
      Teuchos::updateParametersFromXmlFile(inputFile, piroParams.ptr());
      auto stratParams = Piro::extractStratimikosParams(piroParams);
      out << "\nPiro Params" <<std::endl;
      piroParams->print();
      Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
      out << "\nStratimikos Params" <<std::endl;
      stratParams->print();
      linearSolverBuilder.setParameterList(stratParams);
      auto lowsFactory = createLinearSolveStrategy(linearSolverBuilder);
      auto modelWithSolve = rcp(new Thyra::DefaultModelEvaluatorWithSolveFactory<double>(model, lowsFactory));

      auto observer = rcp(new Piro::Test::MockObserver<double>());

      Piro::SolverFactory solverFactory;
      auto solver =  solverFactory.createSolver<double>(piroParams, modelWithSolve, Teuchos::null, observer);

      Teuchos::ParameterList emptyList;
      Teuchos::Array<RCP<Thyra::VectorBase<double>>> responses;
      Teuchos::Array<Teuchos::Array<RCP<Thyra::MultiVectorBase<double>>>> sensitivities;
      Piro::PerformSolve(*solver, emptyList, responses, sensitivities);

      // Print out everything
      out << "Finished Model Evaluation: Printing everything {Exact in brackets}" 
          << "\n-----------------------------------------------------------------"
          << std::setprecision(9) << std::endl;
            
       // Extract default input parameters
      const RCP<const Thyra::VectorBase<double>> p = solver->getNominalValues().get_p(0);

      // Extract output arguments
      const RCP<Thyra::VectorBase<double>> g = responses[0];
      const RCP<Thyra::VectorBase<double>> x = responses[1];      
      
      if (Teuchos::is_null(p)) {
        out << "\nError: parameters pointer is null" << std::endl;
        status += 33;
      } else {          
        const Thyra::DetachedVectorView<double> p_view(p->clone_v());
        out << "\nParameter! {1.0}\n" << p_view(0) << std::endl;
        double p_exact = 1;

        double diff =  std::abs(p_view(0) - p_exact);
        if (diff > tol) {
              status += 100;
              out << "\nExpected parameter value is: "
                  << p_exact  << ", but the value is: "
                  << p_view(0) << ".\n"
                  << "Absolute Difference: " << diff << " > tol: " << tol << std::endl;
          }
      }
      if (Teuchos::is_null(g)) {
        out << "\nError: Responses pointer is null" << std::endl;
        status += 33;
      } else {
          const Thyra::DetachedVectorView<double> g_view(g);
          double g_exact = 0.18;
          out << "\nResponse! {0.18}\n" <<  g_view(0) << std::endl;

          double diff = std::abs(g_exact -  g_view(0));
          if (diff > tol) {
              status += 100;
              out << "\nResponse value is: {"
                  << g_exact << "}, but the computed value is: {"
                  <<  g_view(0) << "}.\n"
                  << "Absolute difference: " << diff << " > tol: " << tol << std::endl;
          }
      }
      if (Teuchos::is_null(x)) {
        out << "\nError: solution pointer is null" << std::endl;
        status += 33;
      } else {
        const Thyra::DetachedVectorView<double> x_view(x);
        out << "\nSolution! {0.18}\n" << x_view(0) << std::endl;
        
        double x_exact = 0.18;
        double diff =  std::abs(x_view(0) - x_exact);
        if (diff > tol) {
            status += 100;
            out << "\nPiro_AnalysisDrvier:  Expected solution vector is : "
                << x_exact << ", but computed solution vector is: "
                << x_view(0) << "\n"
                << "Absolute difference: " << diff << " > tol: " << tol << std::endl;
        }
      }

      out <<
          "\n-----------------------------------------------------------------\n";
    }
    TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
    if (!success) status += 1000; 

    overall_status += status;

  }

  if (overall_status==0) out << "\nTEST PASSED\n" << std::endl;
  else out << "\nTEST Failed:  " << overall_status << std::endl;


  return status;
}
