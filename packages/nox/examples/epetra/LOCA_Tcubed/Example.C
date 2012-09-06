//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
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

// Required for reading and writing parameter lists from xml format
// Configure Trilinos with --enable-teuchos-extended
#ifdef HAVE_TEUCHOS_EXTENDED
#include "Teuchos_XMLParameterListHelpers.hpp"
#endif

using namespace std;

int main(int argc, char *argv[])
{
  int ierr = 0;
  
  // scale factor to test arc-length scaling
  double scale = 1.0;

  // Create output file to save solutions
  std::ofstream outFile("t3.dat");
  outFile.setf(ios::scientific, ios::floatfield);
  outFile.precision(14);

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
  int NumGlobalElements = 100 + 1;

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
  FiniteElementProblem Problem(NumGlobalElements, Comm, scale, &outFile);

  // Get the vector from the Problem
  Epetra_Vector& soln = Problem.getSolution();

  // Initialize Solution
  soln.PutScalar(1.0);
  
  // Begin LOCA Solver ************************************

  // Create parameter list
  Teuchos::RCP<Teuchos::ParameterList> paramList = 
    Teuchos::rcp(new Teuchos::ParameterList);

  // Create LOCA sublist
  Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

  // Create the stepper sublist and set the stepper parameters
  Teuchos::ParameterList& locaStepperList = locaParamsList.sublist("Stepper");
  locaStepperList.set("Continuation Method", "Arc Length");      // Default
  locaStepperList.set("Continuation Parameter", "Right BC");
  locaStepperList.set("Initial Value", 0.1/scale);
  locaStepperList.set("Max Value", 100.0/scale);
  locaStepperList.set("Min Value", 0.05/scale);
  locaStepperList.set("Max Steps", 30);
  locaStepperList.set("Max Nonlinear Iterations", 15);

  // Use Housholder method for solving bordered arclength equations
  locaStepperList.set("Bordered Solver Method", "Householder");

#ifdef HAVE_LOCA_ANASAZI
  // Create Anasazi Eigensolver sublist (needs --with-loca-anasazi)
  locaStepperList.set("Compute Eigenvalues",true);
  Teuchos::ParameterList& aList = locaStepperList.sublist("Eigensolver");
  aList.set("Method", "Anasazi");
  aList.set("Block Size", 1);        // Size of blocks
  aList.set("Num Blocks", 10);       // Size of Arnoldi factorization
  aList.set("Num Eigenvalues", 3);   // Number of eigenvalues
  aList.set("Convergence Tolerance", 2.0e-7);          // Tolerance
  aList.set("Step Size", 1);         // How often to check convergence
  aList.set("Maximum Restarts",2);   // Maximum number of restarts
  aList.set("Verbosity",  
	    Anasazi::Errors + 
	    Anasazi::Warnings +
	    Anasazi::FinalSummary);        // Verbosity
#else
    locaStepperList.set("Compute Eigenvalues",false);
#endif
  
  // Create predictor sublist
  Teuchos::ParameterList& predictorList = locaParamsList.sublist("Predictor");
  predictorList.set("Method", "Secant");      // Default

  // Create step size sublist
  Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
  stepSizeList.set("Method", "Adaptive");    // Dfault
  stepSizeList.set("Initial Step Size", 0.1/scale);
  stepSizeList.set("Min Step Size", 1.0e-3/scale);
  stepSizeList.set("Max Step Size", 2000.0/scale);
  stepSizeList.set("Aggressiveness", 0.1);

  // Create the "Solver" parameters sublist to be used with NOX Solvers
  Teuchos::ParameterList& nlParams = paramList->sublist("NOX");

  // Create the NOX printing parameter list
  Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
  nlPrintParams.set("MyPID", MyPID); 
  nlPrintParams.set("Output Information", 
		    NOX::Utils::OuterIteration + 
		    NOX::Utils::OuterIterationStatusTest + 
		    NOX::Utils::InnerIteration +
		    NOX::Utils::Details + 
		    NOX::Utils::Warning + 
		    NOX::Utils::StepperIteration +
		    NOX::Utils::StepperDetails +
		    NOX::Utils::StepperParameters);

  // Create the "Linear Solver" sublist for Newton's method
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newParams = dirParams.sublist("Newton");
  Teuchos::ParameterList& lsParams = newParams.sublist("Linear Solver");
  lsParams.set("Aztec Solver", "GMRES");  
  lsParams.set("Max Iterations", 100);  
  lsParams.set("Tolerance", 1e-4);
  lsParams.set("Output Frequency", 1);    
  lsParams.set("Scaling", "None");             
  lsParams.set("Preconditioner", "Ifpack");
  //lsParams.set("Preconditioner", "AztecOO");
  //lsParams.set("Jacobian Operator", "Matrix-Free");
  //lsParams.set("Preconditioner Operator", "Finite Difference");
  //lsParams.set("Aztec Preconditioner", "ilut"); 
  //lsParams.set("Overlap", 2);   
  //lsParams.set("Fill Factor", 2.0); 
  //lsParams.set("Drop Tolerance", 1.0e-12);

  // Create and initialize the parameter vector
  LOCA::ParameterVector pVector;
  pVector.addParameter("Nonlinear Factor",1.0);
  pVector.addParameter("Left BC", 0.0);
  pVector.addParameter("Right BC", 0.1);

  // Create the interface between the test problem and the nonlinear solver
  // This is created by the user using inheritance of the abstract base class:
  Teuchos::RCP<Problem_Interface> interface = 
    Teuchos::rcp(new Problem_Interface(Problem));
  Teuchos::RCP<LOCA::Epetra::Interface::Required> iReq = interface;
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;
  
  // Create the Epetra_RowMatrixfor the Jacobian/Preconditioner
  Teuchos::RCP<Epetra_RowMatrix> Amat = 
    Teuchos::rcp(&Problem.getJacobian(),false);

  // Create the linear systems
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linsys = 
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(nlPrintParams, lsParams,
						      iReq, iJac, Amat, soln));

  // Create the loca vector
  NOX::Epetra::Vector locaSoln(soln);

  // Create Epetra factory
  Teuchos::RCP<LOCA::Abstract::Factory> epetraFactory =
    Teuchos::rcp(new LOCA::Epetra::Factory);

  // Create global data object
  Teuchos::RCP<LOCA::GlobalData> globalData = 
    LOCA::createGlobalData(paramList, epetraFactory);

  // Create the Group
  Teuchos::RCP<LOCA::Epetra::Group> grp = 
    Teuchos::rcp(new LOCA::Epetra::Group(globalData, nlPrintParams, 
					 iReq, locaSoln, linsys,
					 pVector));
  grp->computeF();

  // Create the Solver convergence test
  Teuchos::RCP<NOX::StatusTest::NormF> wrms = 
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = 
    Teuchos::rcp(new NOX::StatusTest::MaxIters(locaStepperList.get("Max Nonlinear Iterations", 10)));
  Teuchos::RCP<NOX::StatusTest::Combo> combo = 
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(wrms);
  combo->addStatusTest(maxiters);

#ifdef HAVE_TEUCHOS_EXTENDED
  // Write the parameter list to a file
  cout << "Writing parameter list to \"input.xml\"" << endl;
  Teuchos::writeParameterListToXmlFile(*paramList, "input.xml");

  // Read in the parameter list from a file
  cout << "Reading parameter list from \"input.xml\"" << endl;
  Teuchos::RCP<Teuchos::ParameterList> paramList2 = 
    Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::updateParametersFromXmlFile("input.xml", paramList2.ptr());
  paramList = paramList2;
#endif

  // Create the stepper  
  LOCA::Stepper stepper(globalData, grp, combo, paramList);
  LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

  if (status == LOCA::Abstract::Iterator::Finished) 
    globalData->locaUtils->out() << "\nAll tests passed" << endl;
  else {
    if (globalData->locaUtils->isPrintType(NOX::Utils::Error))
      globalData->locaUtils->out() 
	<< "Stepper failed to converge!" << std::endl;
    }

  // Output the parameter list
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperParameters)) {
    globalData->locaUtils->out() 
      << std::endl << "Final Parameters" << std::endl
      << "****************" << std::endl;
    stepper.getList()->print(globalData->locaUtils->out());
    globalData->locaUtils->out() << std::endl;
  }

  LOCA::destroyGlobalData(globalData);

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr ;
}
