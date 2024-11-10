// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
                                                                                
// 1D Finite Element Test Problem
/* Solves the nonlinear equation:
 *
 * d2u 
 * --- - k * u**2 = 0
 * dx2
 *
 * subject to @ x=0, u=1
 */
// see documentation for the ml-class ML_NOX::ML_Nox_Preconditioner.
// compiling and running this example requires a bunch of packages:
// --enable-nox --enable-prerelease
// --enable-nox-epetra
// --enable-epetra
// --enable-epetraext
// --enable-teuchos
// --enable-ml --with-ml_nox
// --enable-aztecoo
// --enable-amesos
// due to cross-dependencies with nox, it is necessary to 
// configure Trilinos and make install WITHOUT the --with-ml_nox option, then 
// configure Trilinos and make install WITH the --with-ml_nox option again.
// The example should now work, if not, contact Michael Gee mwgee@sandia.gov.
//
// usage:
// ml_nox_1Delasticity_example.exe <number_of_elements>
// Try a large number (otherwise You don't get several levels), let's say
// at least 50000

// ml objects
#include "ml_common.h"
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

// ml objects
#include "nlnml_preconditioner.H"
#include "nlnml_linearsystem.H"

// NOX Objects
#include "NOX.H"
#include "NOX_Epetra.H"

// Trilinos Objects
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"

// User's application specific files 
#include "Problem_Interface.H"    // Interface to application
#include "FiniteElementProblem.H" // the application            

using namespace std;
using namespace Teuchos;

int main(int argc, char *argv[])
{
  // Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  // Get the process ID and the total number of processors
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  // Get the number of elements from the command line
  int NumGlobalElements;
  switch(argc) {
    case 2:
      NumGlobalElements = atoi(argv[1]) + 1;
      break;
    case 1:
      NumGlobalElements = 5001;
      break;
    default:
      cout << "Usage: " << argv[0] << " <number_of_elements>" << endl;
#     ifdef HAVE_MPI
      MPI_Finalize() ;
#     endif
      exit(EXIT_FAILURE);
      break;
  }
  cout << NumGlobalElements << endl;

  // The number of unknowns must be at least equal to the 
  // number of processors.
  if (NumGlobalElements < NumProc) {
    cout << "numGlobalBlocks = " << NumGlobalElements 
	 << " cannot be < number of processors = " << NumProc << endl;
#   ifdef HAVE_MPI
    MPI_Finalize() ;
#   endif
    exit(EXIT_FAILURE);
  }

  // Create the FiniteElementProblem class.  This creates all required
  // Epetra objects for the problem and allows calls to the 
  // function (RHS) and Jacobian evaluation routines.
  FiniteElementProblem nlnproblem(NumGlobalElements, Comm);

  // Get the vector from the Problem
  Epetra_Vector& soln = nlnproblem.getSolution();

  // Initialize Solution
  soln.PutScalar(1.0);
  
  // evaluate the nonlinear function once to be safe
  {
     Epetra_Vector* rhs = new Epetra_Vector(Copy,soln,0);
     nlnproblem.evaluate(ALL,&soln,rhs,NULL);
     delete rhs; rhs = 0;
  }
  
  
  // Create the interface between the test problem and the nonlinear solver
  // This is created by the user using inheritance of the abstract base class:
  // Ml_Nox_Fineinterface (which inherits from NOX::EpetraNew:: - interfaces)
  Teuchos::RefCountPtr<Problem_Interface> fineinterface = 
    rcp(new Problem_Interface(nlnproblem));

  // Begin Nonlinear Solver ************************************
   // Create the top level nox parameter list
   Teuchos::ParameterList nlParams;
   RefCountPtr<Teuchos::ParameterList> rcpparams = rcp(&nlParams);
   rcpparams.release();

   // Set the printing parameters in the "Printing" sublist
   Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
   printParams.set("MyPID", MyPID); 
   printParams.set("Output Precision", 9);
   printParams.set("Output Processor", 0);
   printParams.set("Output Information", 
  			    NOX::Utils::OuterIteration + 
			    NOX::Utils::Warning
                            );

   //----------------- Newton------------------------------
   Teuchos::ParameterList* lsparams = 0;
   bool isnlnCG = true;
   if (!isnlnCG)
   {
     // Set the nonlinear solver method
     nlParams.set("Nonlinear Solver", "Line Search Based");
     
     // Sublist for line search 
     Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
     searchParams.set("Method", "Full Step");
     
     // Sublist for direction
     Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
     dirParams.set("Method", "Newton");
     Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
     newtonParams.set("Forcing Term Method", "Constant");
     //newtonParams.set("Forcing Term Method", "Type 1");
     //newtonParams.set("Forcing Term Method", "Type 2");
     newtonParams.set("Forcing Term Minimum Tolerance", 1.0e-6);
     newtonParams.set("Forcing Term Maximum Tolerance", 0.1);
     
     Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
     lsParams.set("Aztec Solver", "CG"); 
     lsParams.set("Max Iterations", 100);  
     lsParams.set("Tolerance", 1e-7);
     lsParams.set("Output Frequency", 50);   
     lsparams = &lsParams;
     
     //lsParams.set("Preconditioning", "AztecOO: Jacobian Matrix");
     lsParams.set("Preconditioning", "User Supplied Preconditioner");
     lsParams.set("Preconditioner","User Defined");

     lsParams.set("Aztec Preconditioner", "ilu");
     lsParams.set("Graph Fill", 2);
     lsParams.set("Fill Factor", 1);
     //-----------------------------------------------------------
   }
   else
   {
     //-----------------nonlinearCG------------------------------
     // Set the nonlinear solver method as line search
     nlParams.set("Nonlinear Solver", "Line Search Based");

     // get sublist for type of linesearch
     Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
     
     // set the nonlinearCG method
     searchParams.set("Method", "NonlinearCG");

     // Sublist for direction
     Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
     dirParams.set("Method", "NonlinearCG");
     
     // sublist for nlnCG params
     Teuchos::ParameterList& nlcgParams = dirParams.sublist("Nonlinear CG");
     nlcgParams.set("Restart Frequency", 500);
     //nlcgParams.set("Precondition", "Off");
     nlcgParams.set("Precondition", "On");
     nlcgParams.set("Orthogonalize", "Polak-Ribiere");
     //nlcgParams.set("Orthogonalize", "Fletcher-Reeves");
     nlcgParams.set("Restart Frequency", 25);
     
     Teuchos::ParameterList& lsParams = nlcgParams.sublist("Linear Solver");
     lsParams.set("Aztec Solver", "CG"); 
     lsParams.set("Max Iterations", 1);  
     lsParams.set("Tolerance", 1e-11);
     lsParams.set("Output Frequency", 50);   
     //lsParams.set("Preconditioning", "None");
     lsParams.set("Preconditioning", "User Supplied Preconditioner");
     
     // EpetraNew takes this as user supplied preconditioner
     lsParams.set("Preconditioner","User Defined");
     
     //lsParams.set("Preconditioning", "AztecOO: Jacobian Matrix");
     lsParams.set("Aztec Preconditioner", "ilu");
     lsParams.set("Graph Fill", 0);
     lsParams.set("Fill Factor", 1);
   }
  // End Nonlinear Solver ************************************


  // Begin Preconditioner ************************************
   ParameterList mlparams;
   mlparams.set("nlnML output",                                      6         ); // ML-output-level (0-10)
   mlparams.set("nlnML max levels",                                  10         ); // max. # levels (minimum = 2 !)
   mlparams.set("nlnML coarse: max size",                            5000        ); // the size ML stops generating coarser levels
   mlparams.set("nlnML is linear preconditioner",                    false       );
   mlparams.set("nlnML apply constraints",                           false       );
   mlparams.set("nlnML is matrixfree",                               false       ); 
   mlparams.set("nlnML finite difference fine level",                false       );
   mlparams.set("nlnML finite difference alpha",                     1.0e-08    );    
   mlparams.set("nlnML finite difference beta",                      1.0e-07    );    
   mlparams.set("nlnML finite difference centered",                  false      );     
   mlparams.set("nlnML Jacobian fix diagonal",                       false       );     

   mlparams.set("nlnML absolute residual tolerance",                 1.0e-06    );
   mlparams.set("nlnML max cycles",                                  500        );
   mlparams.set("nlnML adaptive recompute",                          0.0        ); // recompute if residual is larger then this value
   mlparams.set("nlnML offset recompute",                            0          ); // every offset this preconditioner is recomputed     
   mlparams.set("nlnML additional adaptive nullspace",               0          ); // compute adaptive nullspace (additional kernel vectors)
   mlparams.set("nlnML PDE equations",                               1          ); // dof per node
   mlparams.set("nlnML null space: dimension",                       1          ); // dimension of nullspace
   mlparams.set("nlnML spatial dimension",                           1          );
   mlparams.set("nlnML coarse: type",                                "Uncoupled"); // Uncoupled METIS VBMETIS
   mlparams.set("nlnML nodes per aggregate",                         3          ); // # nodes per agg for coarsening METIS and VBMETIS

   mlparams.set("nlnML use nlncg on fine level",                     true); // use nlnCG or mod. Newton's method   
   mlparams.set("nlnML use nlncg on medium level",                   true);    
   mlparams.set("nlnML use nlncg on coarsest level",                 false);    
   
   mlparams.set("nlnML max iterations newton-krylov fine level",     0); // # iterations of lin. CG in mod. Newton's method    
   mlparams.set("nlnML max iterations newton-krylov medium level" ,  0);    
   mlparams.set("nlnML max iterations newton-krylov coarsest level", 15);    

   mlparams.set("nlnML linear smoother type fine level",             "MLS"); // SGS BSGS Jacobi MLS Bcheby AmesosKLU   
   mlparams.set("nlnML linear smoother type medium level",           "MLS"); 
   mlparams.set("nlnML linear smoother type coarsest level",         "AmesosKLU"); 
   mlparams.set("nlnML linear smoother sweeps fine level",           24);
   mlparams.set("nlnML linear smoother sweeps medium level",         24);
   mlparams.set("nlnML linear smoother sweeps coarsest level",       1);

   mlparams.set("nlnML nonlinear presmoothing sweeps fine level",    1);
   mlparams.set("nlnML nonlinear presmoothing sweeps medium level",  0);
   mlparams.set("nlnML nonlinear smoothing sweeps coarse level",     2);
   mlparams.set("nlnML nonlinear postsmoothing sweeps medium level", 2);
   mlparams.set("nlnML nonlinear postsmoothing sweeps fine level",   3);

   RefCountPtr<NLNML::NLNML_Preconditioner> Prec = 
     rcp(new NLNML::NLNML_Preconditioner(fineinterface,mlparams,Comm));
 
  // End Preconditioner **************************************


  // run the preconditioner as a solver **********************
#if 0
   // the preconditioner can also act as a multigrid solver (without outer Krylov method)
   if (mlparams.get("nlnML is linear preconditioner",true)==false)
   {
      int printlevel = mlparams.get("nlnML output",6);
      int notconverged = Prec->solve();
      double appltime = fineinterface->getsumtime();
      if (printlevel>0 && Comm.MyPID()==0)
      {
         cout << "NOX/ML :===========of which time in ccarat : " << appltime << " sec\n";
         cout << "NOX/ML :======number calls to computeF in this solve : " << fineinterface->getnumcallscomputeF() << "\n\n\n";
      }
      fineinterface->resetsumtime();
      fineinterface->setnumcallscomputeF(0);
      return notconverged;
   }
#endif
  // End run the preconditioner as a solver *******************

   // creat initial guess
   NOX::Epetra::Vector initialGuess(soln);


   // for nlnCG:
   NLNML::NLNML_LinearSystem*              linSys       = 0;
   NOX::Epetra::LinearSystemAztecOO*       azlinSys     = 0;

   if (isnlnCG)
   {
     RefCountPtr<NOX::Epetra::MatrixFree> B = 
       rcp(new NOX::Epetra::MatrixFree(printParams,fineinterface,soln,false));
     // iPrec  = &Prec;
     // iJac   = B;
     // iReq   = fineinterface.get();
     bool matrixfree    = mlparams.get("nlnML is matrixfree",false);
     int  ml_printlevel = mlparams.get("nlnML output",6);
     linSys = new NLNML::NLNML_LinearSystem(B,B,Prec,null,Prec,matrixfree,0,ml_printlevel);
   }
   else
   {
     // for Newton:
     RefCountPtr<Epetra_CrsMatrix> A = rcp(fineinterface->getJacobian());
     A.release();
     azlinSys = new NOX::Epetra::LinearSystemAztecOO(printParams,*lsparams,
                                                     fineinterface,A,Prec,
                                                     Prec,initialGuess);
   }

   
   RefCountPtr<NLNML::NLNML_LinearSystem> rcplinsys = rcp(linSys);
   RefCountPtr<NOX::Epetra::LinearSystemAztecOO> rcpazlinsys = rcp(azlinSys);
   
   // Create the Group
   NOX::Epetra::Group* grp = 0;
   if (isnlnCG)
      grp = new NOX::Epetra::Group(printParams,fineinterface,initialGuess,rcplinsys); 
   else 
      grp = new NOX::Epetra::Group(printParams,fineinterface,initialGuess,rcpazlinsys); 

   RefCountPtr<NOX::Epetra::Group> rcpgrp = rcp(grp);

   // Create the convergence tests
   double FAS_normF = mlparams.get("nlnML absolute residual tolerance",1.0e-06);
   RefCountPtr<NOX::StatusTest::NormF> absresid 
     = rcp(new NOX::StatusTest::NormF(FAS_normF));
   RefCountPtr<NOX::StatusTest::NormUpdate> nupdate 
     = rcp(new NOX::StatusTest::NormUpdate(FAS_normF));
   RefCountPtr<NOX::StatusTest::Combo> converged 
     = rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
   converged->addStatusTest(absresid);
   converged->addStatusTest(nupdate);
   
   RefCountPtr<NOX::StatusTest::FiniteValue> fv = 
     rcp(new NOX::StatusTest::FiniteValue());
   RefCountPtr<NOX::StatusTest::MaxIters> maxiters = 
     rcp(new NOX::StatusTest::MaxIters(250));
   RefCountPtr<NOX::StatusTest::Combo> combo = 
     rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
   combo->addStatusTest(maxiters);
   combo->addStatusTest(converged);
   combo->addStatusTest(fv);

   // Create the method
   RCP<NOX::Solver::Generic> rcpsolver = NOX::Solver::buildSolver(rcpgrp,combo,rcpparams);
#if 0
   NOX::Solver::Manager solver(rcpgrp, combo, rcpparams);
   RefCountPtr<NOX::Solver::Manager> rcpsolver = 
     rcp(&solver);
   rcpsolver.release();
#endif
   // register the solver class in the preconditioner in case of nonlinear preconditioning
   Prec->SetNoxSolver(rcpsolver);


   // solve
   double t0 = GetClock();
   NOX::StatusTest::StatusType status = rcpsolver->solve();
   double t1 = GetClock();
   int  ml_printlevel = mlparams.get("nlnML output",6);
   if (ml_printlevel>0 && Comm.MyPID()==0)
      cout << "NOX/ML :============solve time incl. setup : " << (t1-t0) << " sec\n";
   double appltime = fineinterface->getsumtime();
   if (ml_printlevel>0 && Comm.MyPID()==0)
   {
      cout << "NOX/ML :===========of which time in application : " << appltime << " sec\n";
      cout << "NOX/ML :======number calls to computeF in this solve : " 
           << fineinterface->getnumcallscomputeF() << "\n\n\n";
   }
   fineinterface->resetsumtime();
   fineinterface->setnumcallscomputeF(0);
   
   if (status != NOX::StatusTest::Converged) 
   {
      if (Comm.MyPID()==0)
         cout << "***WRN***: ML/NOX not converged!";
#     ifdef HAVE_MPI
      MPI_Finalize() ;
#     endif
      exit(EXIT_FAILURE);
   }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  exit(EXIT_SUCCESS);

} /* end main */



#else  
// ml objects
#include <iostream>
#include "ml_common.h"
int main(int argc, char *argv[])
{
  std::cout << "running ml_nox_1Delasticity_example.exe needs: \n"
       << "defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)\n"
       << "see documentation for the ml-class ML_NOX::ML_Nox_Preconditioner\n";
  std::cout.flush();
  return 0;
}
#endif 
