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

#include "Piro_Epetra_SolverFactory.hpp"
#include "Piro_ProviderHelpers.hpp"
#include "Piro_Epetra_PerformSolve.hpp"

#include "MockModelEval_A.hpp"

#include "Piro_ConfigDefs.hpp"

#ifdef HAVE_PIRO_NOX
#include "SaveEigenData_Epetra.hpp"
#endif /* HAVE_PIRO_NOX */

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

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

#ifdef HAVE_PIRO_NOX
  solverFactory.setSource<LOCA::SaveEigenData::AbstractStrategy>(
      Piro::providerFromReferenceAcceptingConstructor<SaveEigenData_Epetra>());
#endif /* HAVE_PIRO_NOX */

  int numTests=3;
  for (int iTest=0; iTest<numTests; iTest++) {

    if (doAll) {
      switch (iTest) {
       case 0: inputFile="input_Solve_NOX_1.xml"; break;
       case 1: inputFile="input_Solve_NOX_Adjoint.xml"; break;
       case 2: inputFile="input_Solve_LOCA_1.xml"; break;
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
      RCP<EpetraExt::ModelEvaluator> Model = rcp(new MockModelEval_A(appComm));

      RCP<Teuchos::ParameterList> piroParams =
         rcp(new Teuchos::ParameterList("Piro Parameters"));
      Teuchos::updateParametersFromXmlFile(inputFile, piroParams.ptr());

      // Wrap original model into a Piro solver to build a response-only application
      const RCP<EpetraExt::ModelEvaluator> piro = solverFactory.createSolver(piroParams, Model);

      const Teuchos::RCP<Teuchos::ParameterList> solveParams =
        Teuchos::sublist(Teuchos::sublist(piroParams, "Analysis"), "Solve");

      Teuchos::Array<Teuchos::RCP<const Epetra_Vector> > responses;
      Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Epetra_MultiVector> > > sensitivities;
      Piro::Epetra::PerformSolve(*piro, *solveParams, responses, sensitivities);

      // Extract
      const RCP<const Epetra_Vector> g1 = responses[0];
      const RCP<const Epetra_Vector> gx = responses[1];
      const RCP<const Epetra_MultiVector> dgdp = sensitivities[0][0];

      const RCP<const Epetra_Vector> p1 = piro->get_p_init(0);

      // Print out everything
      if (Proc == 0)
        std::cout << "Finished Model Evaluation: Printing everything {Exact in brackets}"
             << "\n-----------------------------------------------------------------"
             << std::setprecision(9) << std::endl;

      p1->Print(std::cout << "\nParameters! {1,1}\n");
      g1->Print(std::cout << "\nResponses! {8.0}\n");
      gx->Print(std::cout << "\nSolution! {1,2,3,4}\n");
      if (Teuchos::nonnull(dgdp))
        dgdp->Print(std::cout <<"\nSensitivities {2.0, -8.0}\n");

      if (Proc == 0)
        std::cout <<
          "\n-----------------------------------------------------------------\n";
    }
    TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
    if (!success) status+=1000;

    overall_status += status;
  }

  if (Proc==0) {
    if (overall_status==0)
      std::cout << "\nTEST PASSED\n" << std::endl;
    else
      std::cout << "\nTEST Failed:  " << overall_status << std::endl;
  }

  return status;
}
