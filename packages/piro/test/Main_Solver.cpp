#include <iostream>
#include <string>

#include "MockModelEval_A.hpp"
#include "SaveEigenData_Epetra.hpp"

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TestForException.hpp"

#include "Piro_ConfigDefs.hpp"

#ifdef Piro_ENABLE_NOX
#include "Piro_Epetra_NOXSolver.hpp"
#include "Piro_Epetra_LOCASolver.hpp"
#endif
#ifdef Piro_ENABLE_Rythmos
#include "Piro_Epetra_RythmosSolver.hpp"
#endif


int main(int argc, char *argv[]) {

  int status=0; // 0 = pass, failures are incremented
  int overall_status=0; // 0 = pass, failures are incremented over multiple tests

  // Initialize MPI and timer
  int Proc=0;
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  double total_time = -MPI_Wtime();
  (void) MPI_Comm_rank(MPI_COMM_WORLD, &Proc);
  MPI_Comm appComm = MPI_COMM_WORLD;
#else
  int appComm=0;
#endif

  using Teuchos::RCP;
  using Teuchos::rcp;
  char* inputFile;

  bool doAll = (argc==1);
  if (argc>1) doAll = !strcmp(argv[1],"-v");
 
  for (int iTest=0; iTest<3; iTest++) {

    if (doAll) {
      switch (iTest) {
       case 0: inputFile="input_Solve_NOX_1.xml"; break;
       case 1: inputFile="input_Solve_LOCA_1.xml"; break;
       case 2: inputFile="input_Solve_Rythmos_1.xml"; break;
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

      RCP<Teuchos::ParameterList> piroParams =
         rcp(new Teuchos::ParameterList("Piro Parameters"));
      Teuchos::updateParametersFromXmlFile(inputFile, piroParams.get());

      // Use these two objects to construct a Piro solved application 
      //   EpetraExt::ModelEvaluator is  base class of all Piro::Epetra solvers
      RCP<EpetraExt::ModelEvaluator> piro;

      std::string& solver = piroParams->get("Piro Solver","NOX");
#ifdef Piro_ENABLE_NOX
      if (solver=="NOX")
        piro = rcp(new Piro::Epetra::NOXSolver(piroParams, Model));
      else if (solver=="LOCA") {
        RCP<LOCA::SaveEigenData::AbstractStrategy> saveEigs =
            rcp(new SaveEigenData_Epetra(piroParams->sublist("LOCA")));
        piro = rcp(new Piro::Epetra::LOCASolver(
                       piroParams, Model, Teuchos::null, saveEigs));
      }
      else
#endif
#ifdef Piro_ENABLE_Rythmos
      if (solver=="Rythmos")
        piro = rcp(new Piro::Epetra::RythmosSolver(piroParams, Model));
      else 
#endif
        TEST_FOR_EXCEPTION(true, std::logic_error,
          "Error: Unknown Piro Solver : " << solver);

      bool computeSens = piroParams->get("Compute Sensitivities", false);

      // Now the (somewhat cumbersome) setting of inputs and outputs
      EpetraExt::ModelEvaluator::InArgs inArgs = piro->createInArgs();
      int num_p = inArgs.Np();     // Number of *vectors* of parameters
      RCP<Epetra_Vector> p1 = rcp(new Epetra_Vector(*(piro->get_p_init(0))));
      int numParams = p1->MyLength(); // Number of parameters in p1 vector
      inArgs.set_p(0,p1);

      // Set output arguments to evalModel call
      EpetraExt::ModelEvaluator::OutArgs outArgs = piro->createOutArgs();
      int num_g = outArgs.Ng(); // Number of *vectors* of responses
      RCP<Epetra_Vector> g1 = rcp(new Epetra_Vector(*(piro->get_g_map(0))));
      outArgs.set_g(0,g1);
      // Solution vector is returned as extra respons vector
      RCP<Epetra_Vector> gx = rcp(new Epetra_Vector(*(piro->get_g_map(1))));
      outArgs.set_g(1,gx);

      RCP<Epetra_MultiVector> dgdp = rcp(new Epetra_MultiVector(g1->Map(), numParams));
      if (computeSens) outArgs.set_DgDp(0, 0, dgdp);

      // Now, solve the problem and return the responses
      piro->evalModel(inArgs, outArgs);

      // Print out everything
      if (Proc == 0)
        cout << "Finished Model Evaluation: Printing everything {Exact in brackets}" 
             << "\n-----------------------------------------------------------------"
             << std::setprecision(9) << endl;

      p1->Print(cout << "\nParameters! {1,1}\n");
      g1->Print(cout << "\nResponses! {8.0}\n");
      gx->Print(cout << "\nSolution! {1,2,3,4}\n");
      if (computeSens)
        dgdp->Print(cout <<"\nSensitivities {2.0, -8.0}\n");

      if (Proc == 0)
        cout <<
          "\n-----------------------------------------------------------------\n";
    }

    catch (std::exception& e) {
      cout << e.what() << endl;
      status = 10;
    }
    catch (string& s) {
      cout << s << endl;
      status = 20;
    }
    catch (char *s) {
      cout << s << endl;
      status = 30;
    }
    catch (...) {
      cout << "Caught unknown exception!" << endl;
      status = 40;
    }

    overall_status += status;
  }

#ifdef HAVE_MPI
  total_time +=  MPI_Wtime();
  MPI_Barrier(MPI_COMM_WORLD);
  if (Proc==0) cout << "\n\nTOTAL TIME     " 
                    << total_time << endl;
  MPI_Finalize() ;
#endif

  if (Proc==0) {
    if (overall_status==0) 
      cout << "\nTEST PASSED\n" << endl;
    else 
      cout << "\nTEST Failed:  " << overall_status << endl;
  }

  return status;
}
