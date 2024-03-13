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

#include "Piro_SolverFactory.hpp"
#include "Piro_StratimikosUtils.hpp"
#include "Piro_PerformSolve.hpp"

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"

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

  int numTests=2;
  for (int iTest=0; iTest<numTests; iTest++) {

    if (doAll) {
      switch (iTest) {
       case 0: inputFile="input_Solve_NOX_3.xml"; break;
       case 1: inputFile="input_Solve_LOCA_1.xml"; break;
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
      RCP<EpetraExt::ModelEvaluator> epetraModel = rcp(new MockModelEval_A(appComm));

      RCP<Teuchos::ParameterList> piroParams =
         rcp(new Teuchos::ParameterList("Piro Parameters"));
      Teuchos::updateParametersFromXmlFile(inputFile, piroParams.ptr());

      // Use these two objects to construct a Piro solved application
      RCP<const Thyra::ResponseOnlyModelEvaluatorBase<double> > piro;
      {
        const RCP<Teuchos::ParameterList> stratParams = Piro::extractStratimikosParams(piroParams);

        Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
        linearSolverBuilder.setParameterList(stratParams);

        RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory =
          createLinearSolveStrategy(linearSolverBuilder);

        const RCP<Thyra::ModelEvaluator<double> > thyraModel =
          Thyra::epetraModelEvaluator(epetraModel,lowsFactory);

        piro = solverFactory.createSolver(piroParams, thyraModel);
      }

      const Teuchos::RCP<Teuchos::ParameterList> solveParams =
        Teuchos::sublist(Teuchos::sublist(piroParams, "Analysis"), "Solve");

      Teuchos::Array<RCP<const Thyra::VectorBase<double> > > responses;
      Teuchos::Array<Teuchos::Array<RCP<const Thyra::MultiVectorBase<double> > > > sensitivities;
      Piro::PerformSolve(*piro, *solveParams, responses, sensitivities);

      // Extract default input parameters
      const RCP<const Thyra::VectorBase<double> > p1 = piro->getNominalValues().get_p(0);

      // Extract output arguments
      const RCP<const Thyra::VectorBase<double> > g1 = responses[0];
      const RCP<const Thyra::VectorBase<double> > gx = responses[1];
      const RCP<const Thyra::MultiVectorBase<double> > dgdp = sensitivities[0][0];

      // Print out everything
      if (Proc == 0)
        std::cout << "Finished Model Evaluation: Printing everything {Exact in brackets}"
             << "\n-----------------------------------------------------------------"
             << std::setprecision(9) << std::endl;

      std::cout << "\nParameters! {1,1}\n" << *p1 << std::endl;
      std::cout << "\nResponses! {8.0}\n" << *g1 << std::endl;
      std::cout << "\nSolution! {1,2,3,4}\n" << *gx << std::endl;
      if (Teuchos::nonnull(dgdp))
        std::cout <<"\nSensitivities {2.0, -8.0}\n" << *dgdp << std::endl;

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
