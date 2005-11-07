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
#include "NOX_Parameter_Teuchos2NOX.H"

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
  Teuchos::RefCountPtr<NOX::Parameter::List> paramList = 
    Teuchos::rcp(new NOX::Parameter::List);

  // Create LOCA sublist
  NOX::Parameter::List& locaParamsList = paramList->sublist("LOCA");

  // Create the stepper sublist and set the stepper parameters
  NOX::Parameter::List& locaStepperList = locaParamsList.sublist("Stepper");
  //locaStepperList.setParameter("Continuation Method", "Natural");
  locaStepperList.setParameter("Continuation Method", "Arc Length");
  locaStepperList.setParameter("Bordered Solver Method", "Householder");
  locaStepperList.setParameter("Continuation Parameter", "Right BC");
  //locaStepperList.setParameter("Continuation Parameter", "Nonlinear Factor");
  locaStepperList.setParameter("Initial Value", 0.1/scale);
  locaStepperList.setParameter("Max Value", 100.0/scale);
  locaStepperList.setParameter("Min Value", 0.05/scale);
  locaStepperList.setParameter("Max Steps", 30);
  locaStepperList.setParameter("Max Nonlinear Iterations", 15);
  locaStepperList.setParameter("Enable Arc Length Scaling", true);
  locaStepperList.setParameter("Goal Arc Length Parameter Contribution", 0.5);
  locaStepperList.setParameter("Max Arc Length Parameter Contribution", 0.7);
  locaStepperList.setParameter("Initial Scale Factor", 1.0);
  locaStepperList.setParameter("Min Scale Factor", 1.0e-8);
  locaStepperList.setParameter("Enable Tangent Factor Step Size Scaling",false);
  locaStepperList.setParameter("Min Tangent Factor", 0.8);
  locaStepperList.setParameter("Tangent Factor Exponent",1.5);

  // Create bifurcation sublist
    NOX::Parameter::List& bifurcationList = 
      locaParamsList.sublist("Bifurcation");
    bifurcationList.setParameter("Method", "None");

  // Create Anasazi Eigensolver sublist (needs --with-loca-anasazi)
  locaStepperList.setParameter("Compute Eigenvalues",true);
  NOX::Parameter::List& aList = locaStepperList.sublist("Eigensolver");
  aList.setParameter("Method", "Anasazi");
  aList.setParameter("Block Size", 1);
  aList.setParameter("Arnoldi Size", 10);
  aList.setParameter("NEV", 3);
  aList.setParameter("Tol", 2.0e-7);
  aList.setParameter("Convergence Check", 1);
  aList.setParameter("Restarts",2);
  aList.setParameter("Frequency",1);
  aList.setParameter("Debug Level",0);
  
  // Create predictor sublist
  NOX::Parameter::List& predictorList = locaParamsList.sublist("Predictor");
  //predictorList.setParameter("Method", "Constant");
  predictorList.setParameter("Method", "Tangent");
  //predictorList.setParameter("Method", "Secant");

  // Create step size sublist
  NOX::Parameter::List& stepSizeList = locaParamsList.sublist("Step Size");
  //stepSizeList.setParameter("Method", "Constant");
  stepSizeList.setParameter("Method", "Adaptive");
  stepSizeList.setParameter("Initial Step Size", 0.1/scale);
  stepSizeList.setParameter("Min Step Size", 1.0e-3/scale);
  stepSizeList.setParameter("Max Step Size", 2000.0/scale);
  stepSizeList.setParameter("Aggressiveness", 0.1);
  stepSizeList.setParameter("Failed Step Reduction Factor", 0.5);
  stepSizeList.setParameter("Successful Step Increase Factor", 1.26); // for constant

  // Create the "Solver" parameters sublist to be used with NOX Solvers
  NOX::Parameter::List& nlParams = paramList->sublist("NOX");
  nlParams.setParameter("Nonlinear Solver", "Line Search Based");

  // Create the NOX printing parameter list
  NOX::Parameter::List& nlPrintParams = nlParams.sublist("Printing");
  nlPrintParams.setParameter("MyPID", MyPID); 
  nlPrintParams.setParameter("Output Information", 
			     NOX::Utils::OuterIteration + 
			     NOX::Utils::OuterIterationStatusTest + 
			     NOX::Utils::InnerIteration +
			     NOX::Utils::Parameters + 
			     NOX::Utils::Details + 
			     NOX::Utils::Warning + 
			     NOX::Utils::StepperIteration +
			     NOX::Utils::StepperDetails);

  // Create the "Line Search" sublist for the "Line Search Based" solver
  NOX::Parameter::List& searchParams = nlParams.sublist("Line Search");
  searchParams.setParameter("Method", "Full Step");
  searchParams.setParameter("Max Iters", 7);
  searchParams.setParameter("Default Step", 1.0000);
  searchParams.setParameter("Recovery Step", 0.0001);
  searchParams.setParameter("Minimum Step", 0.0001);

  // Create the "Direction" sublist for the "Line Search Based" solver
  NOX::Parameter::List& dirParams = nlParams.sublist("Direction");
  NOX::Parameter::List& newParams = dirParams.sublist("Newton");
  dirParams.setParameter("Method", "Newton");
  newParams.setParameter("Forcing Term Method", "Constant");

  // Create the "Linear Solver" sublist for the "Direction" sublist
  NOX::Parameter::List& lsParams = newParams.sublist("Linear Solver");
  lsParams.setParameter("Aztec Solver", "GMRES");  
  lsParams.setParameter("Max Iterations", 100);  
  lsParams.setParameter("Tolerance", 1e-4);
  lsParams.setParameter("Output Frequency", 1);    
  lsParams.setParameter("Scaling", "None");             
  lsParams.setParameter("Preconditioner", "Ifpack");
  //lsParams.setParameter("Preconditioner", "AztecOO");
  //lsParams.setParameter("Jacobian Operator", "Matrix-Free");
  //lsParams.setParameter("Preconditioner Operator", "Finite Difference");
  lsParams.setParameter("Aztec Preconditioner", "ilut"); 
  //lsParams.setParameter("Overlap", 2);   
  //lsParams.setParameter("Fill Factor", 2.0); 
  //lsParams.setParameter("Drop Tolerance", 1.0e-12);

  // Create and initialize the parameter vector
  LOCA::ParameterVector pVector;
  pVector.addParameter("Nonlinear Factor",1.0);
  pVector.addParameter("Left BC", 0.0);
  pVector.addParameter("Right BC", 0.1);

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
    Teuchos::rcp(new LOCA::Epetra::Group(globalData, nlPrintParams, 
					 iReq, locaSoln, linsys,
					 pVector));
  grp->computeF();

  // Create the Solver convergence test
  Teuchos::RefCountPtr<NOX::StatusTest::NormF> wrms = 
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
  Teuchos::RefCountPtr<NOX::StatusTest::MaxIters> maxiters = 
    Teuchos::rcp(new NOX::StatusTest::MaxIters(searchParams.getParameter("Max Iters", 10)));
  Teuchos::RefCountPtr<NOX::StatusTest::Combo> combo = 
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(wrms);
  combo->addStatusTest(maxiters);

#ifdef HAVE_TEUCHOS_EXPAT
  // Test writing of param list to XML file, and rereading it
  // into a new parameter list.
  NOX::Parameter::Teuchos2NOX pl_converter;
                                                                                                                                        
  cout << "Writing parameter list to \"input.xml\"" << cout;
  pl_converter.SaveToXMLFile("input.xml", *paramList);
                                                                                                                                        
  cout << "Reading parameter list from \"input.xml\"" << cout;
  Teuchos::RefCountPtr<NOX::Parameter::List> paramList2 = pl_converter.ReadFromXMLFile("input.xml");
#else
  Teuchos::RefCountPtr<NOX::Parameter::List> paramList2 = paramList;
#endif

  // Create the stepper  
  LOCA::NewStepper stepper(globalData, grp, combo, paramList);
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
    stepper.getParameterList()->print(globalData->locaUtils->out());
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
