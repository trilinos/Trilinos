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

using namespace std;

int main(int argc, char *argv[])
{
  int ierr = 0;

  double nonlinear_factor = 1.0;
  double left_bc = 0.0;
  double right_bc = 0.40;

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
    std::cout << "numGlobalBlocks = " << NumGlobalElements 
	 << " cannot be < number of processors = " << NumProc << std::endl;
    exit(1);
  }

  // Create the FiniteElementProblem class.  This creates all required
  // Epetra objects for the problem and allows calls to the 
  // function (RHS) and Jacobian evaluation routines.
  FiniteElementProblem Problem(NumGlobalElements, Comm);

  // Get the vector from the Problem
  Epetra_Vector& soln = Problem.getSolution();

  // Initialize Solution
  soln.PutScalar(0.1);

  // Create initial guess for the null vector of jacobian
  Teuchos::RCP<NOX::Abstract::Vector> solnTwo = 
    Teuchos::rcp(new NOX::Epetra::Vector(soln));  
  solnTwo->init(2.5);             // initial value 1.0
  
  // Begin LOCA Solver ************************************

  // Create parameter list
  Teuchos::RCP<Teuchos::ParameterList> paramList = 
    Teuchos::rcp(new Teuchos::ParameterList);

  // Create LOCA sublist
  Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

  // Create the stepper sublist and set the stepper parameters
  Teuchos::ParameterList& locaStepperList = locaParamsList.sublist("Stepper");
  locaStepperList.set("Continuation Method", "Natural");
  //locaStepperList.set("Bordered Solver Method", "Nested");
  //locaStepperList.set("Bordered Solver Method", "Householder");
  locaStepperList.set("Continuation Parameter", "Nonlinear Factor");
  locaStepperList.set("Initial Value", nonlinear_factor);
  locaStepperList.set("Max Value", 1.6);
  locaStepperList.set("Min Value", 0.00);
  locaStepperList.set("Max Steps", 20);
  locaStepperList.set("Max Nonlinear Iterations", 15);

  // Create bifurcation sublist
  Teuchos::ParameterList& bifurcationList = 
    locaParamsList.sublist("Bifurcation");
  bifurcationList.set("Type", "Phase Transition");
  bifurcationList.set("Bifurcation Parameter", "Right BC");

  bifurcationList.set("Second Solution Vector", solnTwo);
  
  // Create predictor sublist
  Teuchos::ParameterList& predictorList = locaParamsList.sublist("Predictor");
  predictorList.set("Method", "Secant");

  // Create step size sublist
  Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
  stepSizeList.set("Method", "Constant");
  stepSizeList.set("Initial Step Size", 0.1);
  stepSizeList.set("Min Step Size", 1.0e-3);
  stepSizeList.set("Max Step Size", 2000.0);
  stepSizeList.set("Aggressiveness", 0.1);

  // Create the "Solver" parameters sublist to be used with NOX Solvers
  Teuchos::ParameterList& nlParams = paramList->sublist("NOX");

  // Create the NOX printing parameter list
  Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
  nlPrintParams.set("MyPID", MyPID); 
  nlPrintParams.set("Output Precision", 6); 
  nlPrintParams.set("Output Information", 
		    NOX::Utils::OuterIteration + 
		    NOX::Utils::OuterIterationStatusTest + 
		    NOX::Utils::InnerIteration +
		    NOX::Utils::Details + 
		    NOX::Utils::LinearSolverDetails +
		    NOX::Utils::Warning + 
		    NOX::Utils::StepperIteration +
		    NOX::Utils::StepperDetails +
		    NOX::Utils::StepperParameters);

  // Create the "Linear Solver" sublist for Newton's method
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newParams = dirParams.sublist("Newton");
  Teuchos::ParameterList& lsParams = newParams.sublist("Linear Solver");
  lsParams.set("Aztec Solver", "GMRES");  
  lsParams.set("Max Iterations", 200);  
  lsParams.set("Tolerance", 1e-6);
  lsParams.set("Output Frequency", 50);    
  //lsParams.set("Scaling", "None");             
  //lsParams.set("Scaling", "Row Sum");  
  lsParams.set("Compute Scaling Manually", false);
  lsParams.set("Preconditioner", "Ifpack");
  lsParams.set("Ifpack Preconditioner", "ILU");

  //lsParams.set("Preconditioner", "New Ifpack");
  //Teuchos::ParameterList& ifpackParams = lsParams.sublist("Ifpack");
  //ifpackParams.set("fact: level-of-fill", 1);

  // Create and initialize the parameter vector
  LOCA::ParameterVector pVector;
  pVector.addParameter("Nonlinear Factor",nonlinear_factor);
  pVector.addParameter("Left BC", left_bc);
  pVector.addParameter("Right BC", right_bc);

  // Create the interface between the test problem and the nonlinear solver
  // This is created by the user using inheritance of the abstract base class:
  Teuchos::RCP<Problem_Interface> interface = 
    Teuchos::rcp(new Problem_Interface(Problem));
  Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent> iReq = interface;
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;
  
  // Create the Epetra_RowMatrixfor the Jacobian/Preconditioner
  Teuchos::RCP<Epetra_RowMatrix> Amat = 
    Teuchos::rcp(&Problem.getJacobian(),false);

  // Create scaling object
  Teuchos::RCP<NOX::Epetra::Scaling> scaling = Teuchos::null;
//   scaling = Teuchos::rcp(new NOX::Epetra::Scaling);
//   Teuchos::RCP<Epetra_Vector> scalingVector = 
//     Teuchos::rcp(new Epetra_Vector(soln.Map()));
//   //scaling->addRowSumScaling(NOX::Epetra::Scaling::Left, scalingVector);
//   scaling->addColSumScaling(NOX::Epetra::Scaling::Right, scalingVector);

  // Create transpose scaling object
  Teuchos::RCP<NOX::Epetra::Scaling> trans_scaling = Teuchos::null;
//   trans_scaling = Teuchos::rcp(new NOX::Epetra::Scaling);
//   Teuchos::RCP<Epetra_Vector> transScalingVector = 
//     Teuchos::rcp(new Epetra_Vector(soln.Map()));
//   trans_scaling->addRowSumScaling(NOX::Epetra::Scaling::Right, 
// 				  transScalingVector);
//   trans_scaling->addColSumScaling(NOX::Epetra::Scaling::Left, 
// 				  transScalingVector);
  //bifurcationList.set("Transpose Scaling", trans_scaling);

  // Create the linear systems
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linsys = 
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(nlPrintParams, lsParams,
						      iReq, iJac, Amat, soln,
						      scaling));

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
    Teuchos::rcp(new LOCA::Epetra::Group(globalData, nlPrintParams, iReq, 
					 locaSoln, linsys, linsys,
					 pVector));

  // Inject FreeEnergy interface into the group
  Teuchos::RCP<LOCA::Epetra::Interface::FreeEnergy> iFE = interface;
  grp->setFreeEnergyInterface(iFE);

  grp->computeF();

  // Create the Solver convergence test
  //NOX::StatusTest::NormWRMS wrms(1.0e-2, 1.0e-8);
  Teuchos::RCP<NOX::StatusTest::NormF> wrms = 
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-12));
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = 
    Teuchos::rcp(new NOX::StatusTest::MaxIters(locaStepperList.get("Max Nonlinear Iterations", 10)));
  Teuchos::RCP<NOX::StatusTest::Combo> combo = 
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(wrms);
  combo->addStatusTest(maxiters);
  
  // Create the stepper  
  LOCA::Stepper stepper(globalData, grp, combo, paramList);
  LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();
  
  if (status == LOCA::Abstract::Iterator::Finished) 
    globalData->locaUtils->out() << "All tests passed" << std::endl;
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
