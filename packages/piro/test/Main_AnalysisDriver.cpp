#include <iostream>
#include <string>

#include "MockModelEval_A.hpp"

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

#include "Piro_Thyra_PerformAnalysis.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"


int main(int argc, char *argv[]) {

  int status=0; // 0 = pass, failures are incremented
  int overall_status=0; // 0 = pass, failures are incremented over multiple tests

  // Initialize MPI and timer
  int Proc=0;
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  double total_time = -MPI_Wtime();
  (void) MPI_Comm_rank(MPI_COMM_WORLD, &Proc);
#endif
  MPI_Comm appComm = MPI_COMM_WORLD;

  using Teuchos::RCP;
  using Teuchos::rcp;

  char* inputFile;
  bool doAll = (argc==1);
  if (argc>1) doAll = !strcmp(argv[1],"-v");


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
    Teuchos::updateParametersFromXmlFile(inputFile, &appParams);

    Teuchos::ParameterList piroParams = appParams.sublist("Piro");
    Teuchos::ParameterList& analysisParams = appParams.sublist("Analysis");

    // Use these two objects to construct a Piro solved application 
    //   EpetraExt::ModelEvaluator is  base class of all Piro::Epetra solvers
    RCP<EpetraExt::ModelEvaluator> piro;

    std::string& solver = piroParams.get("Piro Solver","NOX");
    const RCP<Teuchos::ParameterList> piroParamsRCP = rcp(&piroParams, false);
#ifdef Piro_ENABLE_NOX
    if (solver=="NOX")
      piro = rcp(new Piro::Epetra::NOXSolver(piroParamsRCP, Model));
    else if (solver=="LOCA")
      piro = rcp(new Piro::Epetra::LOCASolver(piroParamsRCP, Model));
    else
#endif
#ifdef Piro_ENABLE_Rythmos
    if (solver=="Rythmos")
      piro = rcp(new Piro::Epetra::RythmosSolver(piroParamsRCP, Model));
    else 
#endif
      TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error: Unknown Piro Solver : " << solver);
    // END Builder

    Thyra::EpetraModelEvaluator piroThyra;
    piroThyra.initialize(piro, Teuchos::null);

    // Now call the analysis routine
    Piro::Thyra::PerformAnalysis(piroThyra, analysisParams);

    if (Proc == 0) cout << 
      "\nFinished Analysis: Optimum printed above exact soln = {1,3}" << endl;
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
  } catch (...) {
    cout << "Caught unknown exception!" << endl;
    status = 40;
  }

  overall_status += status;
 }  // End loop over tests

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
      cout << "\nTEST Failed\n" << endl;
  }

  return status;
}
