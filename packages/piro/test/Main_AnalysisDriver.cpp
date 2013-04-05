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

#include "MockModelEval_A.hpp"
#include "ObserveSolution_Epetra.hpp"

#include "Piro_Epetra_SolverFactory.hpp"
#include "Piro_ProviderHelpers.hpp"

#include "Piro_Epetra_PerformAnalysis.hpp"

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

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

  Piro::Epetra::SolverFactory solverFactory;

  for (int iTest=0; iTest<3; iTest++) {

    if (doAll) {
      switch (iTest) {
       case 0: inputFile="input_Analysis_Dakota.xml"; break;
       case 1: inputFile="input_Analysis_OptiPack.xml"; break;
       case 2: inputFile="input_Analysis_MOOCHO.xml"; break;
       default : cout << "iTest logic error " << endl; exit(-1);
      }
    }
    else {
      inputFile=argv[1];
      iTest = 999;
    }

    if (Proc==0)
     cout << "===================================================\n"
          << "======  Running input file "<< iTest <<": "<< inputFile <<"\n"
          << "===================================================\n"
          << endl;

    try {

      // Create (1) a Model Evaluator and (2) a ParameterList
      RCP<EpetraExt::ModelEvaluator> Model = rcp(new MockModelEval_A(appComm));

      // BEGIN Builder
      Teuchos::ParameterList appParams("Application Parameters");
      Teuchos::updateParametersFromXmlFile(inputFile, Teuchos::ptr(&appParams));

      Teuchos::ParameterList piroParams = appParams.sublist("Piro");
      Teuchos::ParameterList& analysisParams = appParams.sublist("Analysis");

#ifdef Piro_ENABLE_NOX
      solverFactory.setSource<NOX::Epetra::Observer>(
          Piro::providerFromDefaultConstructor<ObserveSolution_Epetra>());
#endif

      // Use these two objects to construct a Piro solved application
      // EpetraExt::ModelEvaluator is the base class of all Piro::Epetra solvers
      const RCP<Teuchos::ParameterList> piroParamsRCP = rcp(&piroParams, false);
      const RCP<EpetraExt::ModelEvaluator> piro = solverFactory.createSolver(piroParamsRCP, Model);
      // END Builder

      // Call the analysis routine
      RCP<Epetra_Vector> p;
      status = Piro::Epetra::PerformAnalysis(*piro, analysisParams, p);

      if (Teuchos::nonnull(p)) {
        // Can post-process results here
        if (Proc==0) {
          cout << "\nPiro_AnalysisDrvier:  Optimum printed above has exact soln = {1,3}" << endl;
        }
      }

    }
    TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
    if (!success) status+=1000;

    overall_status += status;
  }  // End loop over tests

  if (Proc==0) {
    if (overall_status==0)
      cout << "\nTEST PASSED\n" << endl;
    else
      cout << "\nTEST Failed: " << overall_status << "\n" << endl;
  }

  return status;
}
