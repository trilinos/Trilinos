// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#include <iostream>
#include <string>

#include "MockModelEval_A_Tpetra.hpp"
//#include "ObserveSolution_Epetra.hpp"

#include "Piro_SolverFactory.hpp"
#include "Piro_ProviderHelpers.hpp"

#include "Piro_PerformAnalysis.hpp"

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Thyra_Ifpack2PreconditionerFactory.hpp"
#include "Stratimikos_MueLuHelpers.hpp"
#include "Piro_StratimikosUtils.hpp"
#include "Thyra_DetachedVectorView.hpp"





#include "Piro_ConfigDefs.hpp"

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

  Piro::SolverFactory solverFactory;

  for (int iTest=0; iTest<2; iTest++) {

    if (doAll) {
      switch (iTest) {
       case 0: inputFile="input_Analysis_ROL_Tpetra.xml"; break;
       case 1: inputFile="input_Analysis_ROL_AdjointSensitivities_Tpetra.xml"; break;
       case 2: inputFile="input_Analysis_OptiPack.xml"; break;
       case 3: inputFile="input_Analysis_MOOCHO.xml"; break;
       default : std::cout << "iTest logic error " << std::endl; exit(-1);
      }
    }
    else {
      inputFile=argv[1];
      iTest = 999;
    }

    if (Proc==0)
     std::cout << "===================================================\n"
          << "======  Running input file "<< iTest <<": "<< inputFile <<"\n"
          << "===================================================\n"
          << std::endl;

    try {

      // Create (1) a Model Evaluator and (2) a ParameterList
      const RCP<Thyra::ModelEvaluator<double>> Model = rcp(new MockModelEval_A_Tpetra(appComm));

      // BEGIN Builder
      const RCP<Teuchos::ParameterList> appParams = rcp(new Teuchos::ParameterList("Application Parameters"));
      Teuchos::updateParametersFromXmlFile(inputFile, Teuchos::ptr(appParams.get()));

      const RCP<Teuchos::ParameterList>  piroParams = Teuchos::sublist(appParams,"Piro");
      Teuchos::ParameterList& analysisParams = appParams->sublist("Analysis");


      Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
      typedef Thyra::PreconditionerFactoryBase<double>              Base;
      typedef Thyra::Ifpack2PreconditionerFactory<Tpetra_CrsMatrix> Impl;
      linearSolverBuilder.setPreconditioningStrategyFactory(
          Teuchos::abstractFactoryStd<Base, Impl>(), "Ifpack2");
      Stratimikos::enableMueLu<int, int>(linearSolverBuilder);

      const Teuchos::RCP<Teuchos::ParameterList> stratList = Piro::extractStratimikosParams(piroParams);
      linearSolverBuilder.setParameterList(stratList);


      const RCP<Thyra::LinearOpWithSolveFactoryBase<double>> lowsFactory =
          createLinearSolveStrategy(linearSolverBuilder);

      RCP<Thyra::ModelEvaluator<double>> ModelWithSolve = rcp(new Thyra::DefaultModelEvaluatorWithSolveFactory<double>(
          Model, lowsFactory));

      // Use these two objects to construct a Piro solved application
      // Thyra::ModelEvaluator is the base class of all Piro::Epetra solvers
      const RCP<Thyra::ModelEvaluatorDefaultBase<double>> piro = solverFactory.createSolver(piroParams, ModelWithSolve);
      // END Builder

      // Call the analysis routine
      RCP<Thyra::VectorBase<double>> p;
      status = Piro::PerformAnalysis(*piro, analysisParams, p);

      Thyra::DetachedVectorView<double> p_view(p);
      double p_exact[2] = {1,3};
      double tol = 1e-6;

      double l2_diff = std::sqrt(std::pow(p_view(0)-p_exact[0],2) + std::pow(p_view(1)-p_exact[1],2));
      if(l2_diff > tol) {
        status+=100;
        if (Proc==0) {
          std::cout << "\nPiro_AnalysisDrvier:  Optimum parameter values are: {"
              <<  p_exact[0] << ", " << p_exact[1] << "}, but computed values are: {"
              <<  p_view(0) << ", " << p_view(1) << "}." <<
              "\n                      Difference in l2 norm: " << l2_diff << " > tol: " << tol <<   std::endl;
        }
      }

    }
    TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
    if (!success) status+=1000;

    overall_status += status;
  }  // End loop over tests

  if (Proc==0) {
    if (overall_status==0)
      std::cout << "\nTEST PASSED\n" << std::endl;
    else
      std::cout << "\nTEST Failed: " << overall_status << "\n" << std::endl;
  }

  return status;
}
