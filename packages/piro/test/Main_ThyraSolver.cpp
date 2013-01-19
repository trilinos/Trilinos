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

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"

#include "Piro_ConfigDefs.hpp"

#ifdef Piro_ENABLE_NOX
#include "Piro_NOXSolver.hpp"
//#include "Piro_LOCASolver.hpp"
#endif
#ifdef Piro_ENABLE_Rythmos
#include "Piro_RythmosSolver.hpp"
#endif


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

#ifdef Piro_ENABLE_Rythmos
  int numTests=3;
#else
  int numTests=2;
#endif

#ifdef NO_LOCA_THYRA_YET
  for (int iTest=0; iTest<numTests; iTest++) {
#else
  for (int iTest=0; iTest<numTests; iTest+=2) {
#endif

    if (doAll) {
      switch (iTest) {
       case 0: inputFile="input_Solve_NOX_3.xml"; break;
       case 1: inputFile="input_Solve_LOCA_1.xml"; break;
       case 2: inputFile="input_Solve_Rythmos_2.xml"; break;
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
      RCP<EpetraExt::ModelEvaluator> epetraModel = rcp(new MockModelEval_A(appComm));

      RCP<Teuchos::ParameterList> piroParams =
         rcp(new Teuchos::ParameterList("Piro Parameters"));
      Teuchos::updateParametersFromXmlFile(inputFile, piroParams.ptr());

      // Use these two objects to construct a Piro solved application
      RCP<Thyra::ModelEvaluatorDefaultBase<double> > piro;
      RCP<Thyra::ModelEvaluatorDefaultBase<double> > thyraModel;


      std::string& solver = piroParams->get("Piro Solver","NOX");
      RCP<Teuchos::ParameterList> stratParams;

      if (solver=="NOX" || solver=="LOCA") {
        stratParams = Teuchos::rcp(&(piroParams->sublist("NOX").sublist("Direction").
          sublist("Newton").sublist("Stratimikos Linear Solver").sublist("Stratimikos")),false);
      }
      else if (solver=="Rythmos") {
	piroParams->sublist("Rythmos").set("Nonlinear Solver Type", "NOX");
        stratParams = Teuchos::rcp(&(piroParams->sublist("Rythmos").sublist("Stratimikos")),false);
      }
      Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
      linearSolverBuilder.setParameterList(stratParams);
      RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory =
        createLinearSolveStrategy(linearSolverBuilder);
      thyraModel = Thyra::epetraModelEvaluator(epetraModel,lowsFactory);


#ifdef Piro_ENABLE_NOX
      if (solver=="NOX") {
        piro = rcp(new Piro::NOXSolver<double>(piroParams, thyraModel));
      }
      else
#ifdef NO_LOCA_YET
 if (solver=="LOCA") {
        piro = rcp(new Piro::LOCASolver<double>( piroParams, thyraModel, Teuchos::null));
      }
      else
#endif
#endif
#ifdef Piro_ENABLE_Rythmos
      if (solver=="Rythmos") {
        piro = rcp(new Piro::RythmosSolver<double>(piroParams, thyraModel));
      }
      else
#endif
      {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Error: Unknown Piro Solver : " << solver);
      }

      const bool computeSens = piroParams->get("Compute Sensitivities", false);

      // Now the (somewhat cumbersome) setting of inputs and outputs
      const Thyra::ModelEvaluatorBase::InArgs<double> inArgs = piro->getNominalValues();
      const int num_p = inArgs.Np();     // Number of *vectors* of parameters
      assert (num_p == 1);  // Logic needs to be generalized -- hardwire to 1 p vector in model
      const RCP<const Thyra::VectorBase<double> > p1 = inArgs.get_p(0);
      const int numParams = p1->space()->dim(); // Number of parameters in p1 vector

      // Set output arguments to evalModel call
      Thyra::ModelEvaluatorBase::OutArgs<double> outArgs = piro->createOutArgs();
      const int num_g = outArgs.Ng(); // Number of *vectors* of responses
      assert (num_g == 2);  // Logic needs to be generalized -- hardwire to 1 g vector in model

      const RCP<Thyra::VectorBase<double> > g1 = Thyra::createMember(*piro->get_g_space(0));
      outArgs.set_g(0,g1);

      // Solution vector is returned as extra respons vector
      const RCP<Thyra::VectorBase<double> > gx = Thyra::createMember(*thyraModel->get_x_space());
      outArgs.set_g(1,gx);

      const RCP<Thyra::MultiVectorBase<double> > dgdp =
        Thyra::createMembers(*piro->get_g_space(0),numParams);
      if (computeSens) outArgs.set_DgDp(0, 0, dgdp);

      // Now, solve the problem and return the responses
      piro->evalModel(inArgs, outArgs);

      // Print out everything
      if (Proc == 0)
        cout << "Finished Model Evaluation: Printing everything {Exact in brackets}"
             << "\n-----------------------------------------------------------------"
             << std::setprecision(9) << endl;

      cout << "\nParameters! {1,1}\n" << *p1 << endl;
      cout << "\nResponses! {8.0}\n" << *g1 << endl;
      cout << "\nSolution! {1,2,3,4}\n" << *gx << endl;
      if (computeSens)
        cout <<"\nSensitivities {2.0, -8.0}\n" << *dgdp << endl;

      if (Proc == 0)
        cout <<
          "\n-----------------------------------------------------------------\n";
    }
    TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
    if (!success) status+=1000;

    overall_status += status;
  }

  if (Proc==0) {
    if (overall_status==0)
      cout << "\nTEST PASSED\n" << endl;
    else
      cout << "\nTEST Failed:  " << overall_status << endl;
  }

  return status;
}
