// @HEADER
// ************************************************************************
// 
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#include <iostream>
#include <string>

#include "MockModelEval_A.hpp"

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TestForException.hpp"
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
       case 0: inputFile="input_Solve_NOX_2.xml"; break;
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
      Teuchos::updateParametersFromXmlFile(inputFile, piroParams.get());

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
        TEST_FOR_EXCEPTION(true, std::logic_error,
          "Error: Unknown Piro Solver : " << solver);
      }

      bool computeSens = piroParams->get("Compute Sensitivities", false);

      // Now the (somewhat cumbersome) setting of inputs and outputs
      Thyra::ModelEvaluatorBase::InArgs<double> inArgs = piro->createInArgs();
      const Thyra::ModelEvaluatorBase::InArgs<double> inArgsNominal = piro->getNominalValues();
      int num_p = inArgs.Np();     // Number of *vectors* of parameters
      assert (num_p == 1);  // Logic needs to be generalized -- hardwire to 1 p vector in model
      RCP<Thyra::VectorBase<double> > p1 = Thyra::createMember(*piro->get_p_space(0));
      Thyra::copy(*inArgsNominal.get_p(0), p1.ptr());
      int numParams = p1->space()->dim(); // Number of parameters in p1 vector
      
      inArgs.set_p(0,p1);

      // Set output arguments to evalModel call
      Thyra::ModelEvaluatorBase::OutArgs<double> outArgs = piro->createOutArgs();
      int num_g = outArgs.Ng(); // Number of *vectors* of responses
      assert (num_g == 2);  // Logic needs to be generalized -- hardwire to 1 g vector in model

      RCP<Thyra::VectorBase<double> > g1 = Thyra::createMember(*piro->get_g_space(0));
      outArgs.set_g(0,g1);

      // Solution vector is returned as extra respons vector
      RCP<Thyra::VectorBase<double> > gx = Thyra::createMember(*thyraModel->get_x_space());
      outArgs.set_g(1,gx);

      RCP<Thyra::MultiVectorBase<double> > dgdp =
        Thyra::createMembers(*thyraModel->get_x_space(),numParams);
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
