//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER
                                                                                
// 1D Finite Element Test Problem
/* Solves continuation problem (Parameter c="Right BC")
 *
 * d2u 
 * --- + a * u**3 = 0
 * dx2
 *
 * subject to @ x=0, u=b
 * subject to @ x=1, u=c
 */

// LOCA Objects
#include "LOCA.H"
#include "LOCA_Epetra.H"
#include "LOCA_MultiStepper.H"

// Trilinos Objects
#ifdef HAVE_MPI
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
#include "Problem_Interface.H" // Interface file to NOX
#include "FiniteElementProblem.H"              

using namespace std;

int main(int argc, char *argv[])
{
  int ierr = 0;
  
  // scale factor to test arc-length scaling
  double scale = 1.0;

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
  if (argc!=2) { 
    cout << "Usage: " << argv[0] << " number_of_elements" << endl;
    exit(1);
  }
  int NumGlobalElements = atoi(argv[1]) + 1;

  // The number of unknowns must be at least equal to the 
  // number of processors.
  if (NumGlobalElements < NumProc) {
    cout << "numGlobalBlocks = " << NumGlobalElements 
	 << " cannot be < number of processors = " << NumProc << endl;
    exit(1);
  }

  // Create the FiniteElementProblem class.  This creates all required
  // Epetra objects for the problem and allows calls to the 
  // function (RHS) and Jacobian evaluation routines.
  FiniteElementProblem Problem(NumGlobalElements, Comm, scale);

  // Get the vector from the Problem
  Epetra_Vector& soln = Problem.getSolution();

  // Initialize Solution
  soln.PutScalar(1.0);
  
  // Begin LOCA Solver ************************************

  // Create parameter list
  Teuchos::RefCountPtr<Teuchos::ParameterList> paramList = 
    Teuchos::rcp(new Teuchos::ParameterList);

  // Create LOCA sublist
  Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

  // Create the stepper sublist and set the stepper parameters
  Teuchos::ParameterList& locaStepperList = locaParamsList.sublist("Stepper");
  locaStepperList.set("Continuation Method", "Arc Length");
  locaStepperList.set("Bordered Solver Method", "Householder");
  locaStepperList.set("Number of Continuation Parameters", 2);
  locaStepperList.set("Epsilon", 0.1);
  locaStepperList.set("Max Charts", 10000);
  locaStepperList.set("Verbosity", 1);
  locaStepperList.set("Page Charts", 1);
  locaStepperList.set("Dump Polyhedra", true);
  locaStepperList.set("Dump Centers", false);
  locaStepperList.set("Filename", "MFresults");
  locaStepperList.set("Enable Arc Length Scaling", false);
  locaStepperList.set("Max Nonlinear Iterations", 15);
  locaStepperList.set("Aggressiveness", 0.0);
  locaStepperList.set("Max Solution Component", 6.0);

  // Create sublist for each continuation parameter
  Teuchos::ParameterList& paramList1 = 
    locaStepperList.sublist("Continuation Parameter 1");
  paramList1.set("Parameter Name", "Right BC");
  paramList1.set("Initial Value", 0.1);
  paramList1.set("Max Value", 4.0);
  paramList1.set("Min Value", 0.0);
  paramList1.set("Initial Step Size", 0.1);
  paramList1.set("Max Step Size", 0.2);
  paramList1.set("Min Step Size", 1.0e-3);

  Teuchos::ParameterList& paramList2 = 
    locaStepperList.sublist("Continuation Parameter 2");
  paramList2.set("Parameter Name", "Nonlinear Factor");
  paramList2.set("Initial Value", 1.0);
  paramList2.set("Max Value", 4.0);
  paramList2.set("Min Value", 0.0);
  paramList2.set("Initial Step Size", 0.1);
  paramList2.set("Max Step Size", 0.2);
  paramList2.set("Min Step Size", 1.0e-3);


  // Create bifurcation sublist
  Teuchos::ParameterList& bifurcationList = 
    locaParamsList.sublist("Bifurcation");
  bifurcationList.set("Type", "None");

  // Create predictor sublist
  Teuchos::ParameterList& predictorList = locaParamsList.sublist("Predictor");
  //predictorList.set("Method", "Constant");
  predictorList.set("Method", "Tangent");
  //predictorList.set("Method", "Secant");

  // Create the "Solver" parameters sublist to be used with NOX Solvers
  Teuchos::ParameterList& nlParams = paramList->sublist("NOX");
  nlParams.set("Nonlinear Solver", "Line Search Based");

  // Create the NOX printing parameter list
  Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
  nlPrintParams.set("MyPID", MyPID); 
  nlPrintParams.set("Output Information", 
			     NOX::Utils::OuterIteration + 
			     NOX::Utils::OuterIterationStatusTest + 
			     NOX::Utils::InnerIteration +
			     NOX::Utils::LinearSolverDetails +
			     NOX::Utils::Parameters + 
			     NOX::Utils::Details + 
			     NOX::Utils::Warning + 
			     NOX::Utils::StepperIteration +
			     NOX::Utils::StepperDetails);

  // Create the "Line Search" sublist for the "Line Search Based" solver
  Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
  searchParams.set("Method", "Full Step");
  searchParams.set("Max Iters", 7);
  searchParams.set("Default Step", 1.0000);
  searchParams.set("Recovery Step", 0.0001);
  searchParams.set("Minimum Step", 0.0001);

  // Create the "Direction" sublist for the "Line Search Based" solver
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newParams = dirParams.sublist("Newton");
  dirParams.set("Method", "Newton");
  newParams.set("Forcing Term Method", "Constant");
  //newParams.set("Forcing Term Method", "Type 1");
  //newParams.set("Forcing Term Method", "Type 2");
  //newParams.set("Forcing Term Minimum Tolerance", 1.0e-4);
  //newParams.set("Forcing Term Maximum Tolerance", 0.1);

  // Create the "Linear Solver" sublist for the "Direction" sublist
  Teuchos::ParameterList& lsParams = newParams.sublist("Linear Solver");
  lsParams.set("Aztec Solver", "GMRES");  
  lsParams.set("Max Iterations", 100);  
  lsParams.set("Tolerance", 1e-4);
  lsParams.set("Output Frequency", 50);    
  lsParams.set("Scaling", "None");             
  //lsParams.set("Scaling", "Row Sum");          
  //lsParams.set("Preconditioning", "None");   
  lsParams.set("Preconditioner", "Ifpack");
  //lsParams.set("Preconditioning", "AztecOO: Jacobian Matrix");   
  //lsParams.set("Preconditioning", "AztecOO: User RowMatrix"); 
  //lsParams.set("Preconditioning", "User Supplied Preconditioner");
  //lsParams.set("Aztec Preconditioner", "ilu"); 
  //lsParams.set("Overlap", 2);  
  //lsParams.set("Graph Fill", 2); 
  //lsParams.set("Aztec Preconditioner", "ilut"); 
  //lsParams.set("Overlap", 2);   
  //lsParams.set("Fill Factor", 2.0);   
  //lsParams.set("Drop Tolerance", 1.0e-12);   
  //lsParams.set("Aztec Preconditioner", "Polynomial"); 
  //lsParams.set("Polynomial Order", 6); 

  // Create and initialize the parameter vector
  LOCA::ParameterVector pVector;
  pVector.addParameter("Right BC", 0.1);
  pVector.addParameter("Nonlinear Factor",1.0);
  pVector.addParameter("Left BC", 0.0);

  // Create the interface between the test problem and the nonlinear solver
  // This is created by the user using inheritance of the abstract base class:
  Teuchos::RefCountPtr<Problem_Interface> interface = 
    Teuchos::rcp(new Problem_Interface(Problem));
  Teuchos::RefCountPtr<LOCA::Epetra::Interface::Required> iReq = interface;
  Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian> iJac = interface;
  
  // Create the Epetra_RowMatrixfor the Jacobian/Preconditioner
  Teuchos::RefCountPtr<Epetra_RowMatrix> Amat = 
    Teuchos::rcp(&Problem.getJacobian(),false);

  // Create the linear systems
  Teuchos::RefCountPtr<NOX::Epetra::LinearSystemAztecOO> linsys = 
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(nlPrintParams, lsParams,
						      iReq, iJac, Amat, soln));

  // Create the loca vector
  NOX::Epetra::Vector locaSoln(soln);

  // Create Epetra factory
  Teuchos::RefCountPtr<LOCA::Abstract::Factory> epetraFactory =
    Teuchos::rcp(new LOCA::Epetra::Factory);

  // Create global data object
  Teuchos::RefCountPtr<LOCA::GlobalData> globalData = 
    LOCA::createGlobalData(paramList, epetraFactory);

  // Create the Group
  Teuchos::RefCountPtr<LOCA::Epetra::Group> grp = 
    Teuchos::rcp(new LOCA::Epetra::Group(globalData, nlPrintParams, iReq, 
					 locaSoln, linsys,
					 pVector));
  grp->computeF();

  // Create the Solver convergence test
  //NOX::StatusTest::NormWRMS wrms(1.0e-2, 1.0e-8);
  Teuchos::RefCountPtr<NOX::StatusTest::NormF> wrms = 
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
  Teuchos::RefCountPtr<NOX::StatusTest::MaxIters> maxiters = 
    Teuchos::rcp(new NOX::StatusTest::MaxIters(searchParams.get("Max Iters", 10)));
  Teuchos::RefCountPtr<NOX::StatusTest::Combo> combo = 
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(wrms);
  combo->addStatusTest(maxiters);

  // Create the stepper  
  LOCA::MultiStepper stepper(globalData, grp, combo, paramList);
  LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

  if (status != LOCA::Abstract::Iterator::Finished) {
    if (globalData->locaUtils->isPrintType(NOX::Utils::Error))
      globalData->locaUtils->out() 
	<< "Stepper failed to converge!" << std::endl;
    }

  // Output the parameter list
  if (globalData->locaUtils->isPrintType(NOX::Utils::Parameters)) {
    globalData->locaUtils->out() 
      << std::endl << "Final Parameters" << std::endl
      << "****************" << std::endl;
    stepper.getList()->print(globalData->locaUtils->out());
    globalData->locaUtils->out() << std::endl;
  }

  destroyGlobalData(globalData);

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr ;
}
