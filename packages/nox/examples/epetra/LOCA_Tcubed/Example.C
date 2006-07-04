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
#ifdef HAVE_TEUCHOS_EXPAT
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"
#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#endif

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
  //locaStepperList.set("Continuation Method", "Natural");
  locaStepperList.set("Continuation Method", "Arc Length");
  locaStepperList.set("Bordered Solver Method", "Householder");
  locaStepperList.set("Continuation Parameter", "Right BC");
  //locaStepperList.set("Continuation Parameter", "Nonlinear Factor");
  locaStepperList.set("Initial Value", 0.1/scale);
  locaStepperList.set("Max Value", 100.0/scale);
  locaStepperList.set("Min Value", 0.05/scale);
  locaStepperList.set("Max Steps", 30);
  locaStepperList.set("Max Nonlinear Iterations", 15);
  locaStepperList.set("Enable Arc Length Scaling", true);
  locaStepperList.set("Goal Arc Length Parameter Contribution", 0.5);
  locaStepperList.set("Max Arc Length Parameter Contribution", 0.7);
  locaStepperList.set("Initial Scale Factor", 1.0);
  locaStepperList.set("Min Scale Factor", 1.0e-8);
  locaStepperList.set("Enable Tangent Factor Step Size Scaling",false);
  locaStepperList.set("Min Tangent Factor", 0.8);
  locaStepperList.set("Tangent Factor Exponent",1.5);

  // Create bifurcation sublist
    Teuchos::ParameterList& bifurcationList = 
      locaParamsList.sublist("Bifurcation");
    bifurcationList.set("Type", "None");

  // Create Anasazi Eigensolver sublist (needs --with-loca-anasazi)
  locaStepperList.set("Compute Eigenvalues",true);
  Teuchos::ParameterList& aList = locaStepperList.sublist("Eigensolver");
  aList.set("Method", "Anasazi");
  aList.set("Block Size", 1);
  aList.set("Arnoldi Size", 10);
  aList.set("NEV", 3);
  aList.set("Tol", 2.0e-7);
  aList.set("Convergence Check", 1);
  aList.set("Restarts",2);
  aList.set("Frequency",1);
  aList.set("Debug Level",0);
  
  // Create predictor sublist
  Teuchos::ParameterList& predictorList = locaParamsList.sublist("Predictor");
  //predictorList.set("Method", "Constant");
  predictorList.set("Method", "Tangent");
  //predictorList.set("Method", "Secant");

  // Create step size sublist
  Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
  //stepSizeList.set("Method", "Constant");
  stepSizeList.set("Method", "Adaptive");
  stepSizeList.set("Initial Step Size", 0.1/scale);
  stepSizeList.set("Min Step Size", 1.0e-3/scale);
  stepSizeList.set("Max Step Size", 2000.0/scale);
  stepSizeList.set("Aggressiveness", 0.1);
  stepSizeList.set("Failed Step Reduction Factor", 0.5);
  stepSizeList.set("Successful Step Increase Factor", 1.26); // for constant

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

  // Create the "Linear Solver" sublist for the "Direction" sublist
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
  lsParams.set("Aztec Preconditioner", "ilut"); 
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
    Teuchos::rcp(new NOX::StatusTest::MaxIters(searchParams.get("Max Iters", 10)));
  Teuchos::RefCountPtr<NOX::StatusTest::Combo> combo = 
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(wrms);
  combo->addStatusTest(maxiters);

#ifdef HAVE_TEUCHOS_EXPAT

  // Write the parameter list to a file
  cout << "Writing parameter list to \"input.xml\"" << endl;
  std::string output_filename = "input.xml";
  Teuchos::XMLParameterListWriter xml_converter;
  Teuchos::XMLObject xml_pl = xml_converter.toXML(*paramList);
  ofstream of(output_filename.c_str()); 
  of << xml_pl << endl;
  of.close();
  
  // Read in the parameter list from a file
  cout << "Reading parameter list from \"input.xml\"" << endl;
  Teuchos::FileInputSource fileSrc(output_filename);
  // Convert the file data into an xml object
  Teuchos::XMLObject xml_obj = fileSrc.getObject();
  // Create a xml to teuchos::parameterlist converter
  Teuchos::XMLParameterListReader xml_to_pl_converter;
  // Create a teuchos parameter list
  Teuchos::RefCountPtr<Teuchos::ParameterList> paramList2
    = Teuchos::rcp(new Teuchos::ParameterList);
  // convert the teuchos xml object to a parameter list
  (*paramList2) = xml_to_pl_converter.toParameterList(xml_obj);
  
#else
  Teuchos::RefCountPtr<Teuchos::ParameterList> paramList2 = paramList;
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
