// $Id$ 
// $Source$ 

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
// Questions? Contact Andy Salinger (agsalin@sandia.gov) or Eric Phipps
// (etphipp@sandia.gov), Sandia National Laboratories.
//
// ************************************************************************
//@HEADER

#include "LOCA.H"
#include "LOCA_Epetra.H"

#include "LOCA_BorderedSolver_AbstractStrategy.H"
#include "LOCA_Parameter_SublistParser.H"
#include "NOX_TestCompare.H"

// Trilinos Objects
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"

// User's application specific files 
#include "Problem_Interface.H" 
#include "Tcubed_FiniteElementProblem.H"
#include "LinearConstraint.H"

// Global variables used in main() and testSolve()
Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup> grp;
Teuchos::RefCountPtr<LinearConstraint> constraints;
Teuchos::RefCountPtr<LOCA::Parameter::SublistParser> parsedParams;
Teuchos::RefCountPtr<LOCA::BorderedSolver::AbstractStrategy> bordering;
Teuchos::RefCountPtr<LOCA::BorderedSolver::AbstractStrategy> householder;
Teuchos::RefCountPtr<LOCA::GlobalData> globalData;
Teuchos::RefCountPtr<NOX::TestCompare> testCompare;
Teuchos::RefCountPtr<NOX::Abstract::MultiVector> A;
Teuchos::RefCountPtr<NOX::Abstract::MultiVector> B;
Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> C;
Teuchos::RefCountPtr<NOX::Abstract::MultiVector> F;
Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> G;
Teuchos::RefCountPtr<NOX::Abstract::MultiVector> X_bordering;
Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> Y_bordering;
Teuchos::RefCountPtr<NOX::Abstract::MultiVector> X_householder;
Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> Y_householder;

int  
testSolve(bool flagA, bool flagB, bool flagC, bool flagF, bool flagG,
	  double reltol, double abstol, 
	  const string& testName) {
  int ierr = 0;

  if (globalData->locaUtils->isPrintType(NOX::Utils::TestDetails))
    globalData->locaUtils->out() 
      << std::endl << "***** " << testName << " *****" << std::endl;

  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> a = 
    Teuchos::null;
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> c = 
    Teuchos::null;
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> f = 
    Teuchos::null;
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> g = 
    Teuchos::null;

  if (!flagA)
    a = A;
  if (!flagC)
    c = C;
  if (!flagF)
    f = F;
  if (!flagG)
    g = G;
  constraints->setIsZeroDX(flagB);

  // Set up bordered problem
  bordering->setMatrixBlocks(grp, a, constraints, c);
  bordering->initForSolve();

  householder->setMatrixBlocks(grp, a, constraints, c);
  householder->initForSolve();

  X_bordering->init(0.0);
  Y_bordering->putScalar(0.0);
  X_householder->init(0.0);
  Y_householder->putScalar(0.0);

  // Solve using bordering
  NOX::Abstract::Group::ReturnType borderingStatus = 
    bordering->applyInverse(*(parsedParams->getSublist("Linear Solver")),
			    f.get(), g.get(), *X_bordering, *Y_bordering);
  if (borderingStatus == NOX::Abstract::Group::Failed)
    ++ierr;

  // Solve using householder
  NOX::Abstract::Group::ReturnType householderStatus = 
    householder->applyInverse(*(parsedParams->getSublist("Linear Solver")),
			      f.get(), g.get(), *X_householder, 
			      *Y_householder);
  if (householderStatus == NOX::Abstract::Group::Failed)
    ++ierr;

  ierr += testCompare->testMultiVector(*X_householder, *X_bordering,
				       reltol, abstol, "Solution Component");

  ierr += testCompare->testMatrix(*Y_householder, *Y_bordering, 
				  reltol, abstol, "Parameter Component");

  return ierr;
}
	  

int main(int argc, char *argv[])
{
  int nConstraints = 10;
  int nRHS = 7;

  double left_bc = 0.0;
  double right_bc = 1.0;
  double nonlinear_factor = 1.0;
  int ierr = 0;
  double reltol = 1.0e-8;
  double abstol = 1.0e-8;
  double lstol = 1.0e-11;

  int MyPID;

  try {

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

    // Get the total number of processors
    MyPID = Comm.MyPID();
    int NumProc = Comm.NumProc();
    
    bool verbose = false;
    // Check for verbose output
    if (argc>1)
      if (argv[1][0]=='-' && argv[1][1]=='v') 
	verbose = true;
    
    // Get the number of elements from the command line
    int NumGlobalElements = 0;
    if ((argc > 2) && (verbose))
      NumGlobalElements = atoi(argv[2]) + 1;
    else if ((argc > 1) && (!verbose))
      NumGlobalElements = atoi(argv[1]) + 1;
    else 
      NumGlobalElements = 101;

    // The number of unknowns must be at least equal to the 
    // number of processors.
    if (NumGlobalElements < NumProc) {
      cout << "numGlobalBlocks = " << NumGlobalElements 
	   << " cannot be < number of processors = " << NumProc << endl;
      exit(1);
    }

    // Create parameter list
    Teuchos::RefCountPtr<Teuchos::ParameterList> paramList = 
      Teuchos::rcp(new Teuchos::ParameterList);

    // Create LOCA sublist
    Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

    // Create the constraints list
    Teuchos::ParameterList& constraintsList = 
      locaParamsList.sublist("Constraints");
    constraintsList.set("Bordered Solver Method", "Bordering");

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    Teuchos::ParameterList& nlParams = paramList->sublist("NOX");

    Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.set("MyPID", MyPID);
    if (verbose)
       nlPrintParams.set("Output Information", 
				  NOX::Utils::Error +
				  NOX::Utils::Details +
				  NOX::Utils::OuterIteration + 
				  NOX::Utils::InnerIteration + 
				  NOX::Utils::Warning +
				  NOX::Utils::TestDetails + 
				  NOX::Utils::StepperIteration +
				  NOX::Utils::StepperDetails);
     else
       nlPrintParams.set("Output Information", NOX::Utils::Error);

    // Create the "Direction" sublist for the "Line Search Based" solver
    Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
    Teuchos::ParameterList& newParams = dirParams.sublist("Newton");
    Teuchos::ParameterList& lsParams = newParams.sublist("Linear Solver");
    lsParams.set("Aztec Solver", "GMRES");  
    lsParams.set("Max Iterations", 100);  
    lsParams.set("Tolerance", lstol);
    if (verbose)
      lsParams.set("Output Frequency", 1);
    else
      lsParams.set("Output Frequency", 0);
    lsParams.set("Scaling", "None");             
    lsParams.set("Preconditioner", "Ifpack");
    //lsParams.set("Preconditioner", "AztecOO");
    //lsParams.set("Jacobian Operator", "Matrix-Free");
    //lsParams.set("Preconditioner Operator", "Finite Difference");
    lsParams.set("Aztec Preconditioner", "ilut"); 
    //lsParams.set("Overlap", 2);   
    //lsParams.set("Fill Factor", 2.0); 
    //lsParams.set("Drop Tolerance", 1.0e-12);
    lsParams.set("Max Age Of Prec", -2);

    // Create the FiniteElementProblem class.  This creates all required
    // Epetra objects for the problem and allows calls to the 
    // function (RHS) and Jacobian evaluation routines.
    Tcubed_FiniteElementProblem Problem(NumGlobalElements, Comm);

    // Get the vector from the Problem
    Epetra_Vector& soln = Problem.getSolution();

    // Initialize Solution
    soln.PutScalar(0.0);

    // Create and initialize the parameter vector
    LOCA::ParameterVector pVector;
    pVector.addParameter("Nonlinear Factor",nonlinear_factor);
    pVector.addParameter("Left BC", left_bc);
    pVector.addParameter("Right BC", right_bc);

    // Create the interface between the test problem and the nonlinear solver
    // This is created by the user using inheritance of the abstract base 
    // class:
    Teuchos::RefCountPtr<Problem_Interface> interface = 
      Teuchos::rcp(new Problem_Interface(Problem));
    Teuchos::RefCountPtr<LOCA::Epetra::Interface::Required> iReq = interface;
    Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian> iJac = interface;
    
    // Create the Epetra_RowMatrixfor the Jacobian/Preconditioner
    Teuchos::RefCountPtr<Epetra_RowMatrix> Amat = 
      Teuchos::rcp(&Problem.getJacobian(),false);
    
    // Create the linear systems
    Teuchos::RefCountPtr<NOX::Epetra::LinearSystemAztecOO> linsys = 
      Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(nlPrintParams, 
							lsParams, iReq, iJac, 
							Amat, soln));

    // Create the loca vector
    NOX::Epetra::Vector locaSoln(soln);

    // Create Epetra factory
    Teuchos::RefCountPtr<LOCA::Abstract::Factory> epetraFactory =
      Teuchos::rcp(new LOCA::Epetra::Factory);

     // Create global data object
    globalData = LOCA::createGlobalData(paramList, epetraFactory);

    // Create parsed parameter list
    parsedParams = 
      Teuchos::rcp(new LOCA::Parameter::SublistParser(globalData));
    parsedParams->parseSublists(paramList);

    // Create the Group
    grp = Teuchos::rcp(new LOCA::Epetra::Group(globalData, nlPrintParams, 
					       iReq, locaSoln, 
					       linsys, pVector));
    
    // Change initial guess to a random vector
    Teuchos::RefCountPtr<NOX::Abstract::Vector> xnew = 
      grp->getX().clone();
    xnew->random();
    grp->setX(*xnew);

    // Create the constraints object & constraint param IDs list
    constraints = 
      Teuchos::rcp(new LinearConstraint(nConstraints, LOCA::ParameterVector(),
					locaSoln));

    // Create bordering solver
    bordering
      = globalData->locaFactory->createBorderedSolverStrategy(
				     parsedParams, 
				     parsedParams->getSublist("Constraints"));

    // Change strategy to Householder
    constraintsList.set("Bordered Solver Method", 
				 "Householder");

    // Create householder solver
    householder
      = globalData->locaFactory->createBorderedSolverStrategy(
				     parsedParams, 
				     parsedParams->getSublist("Constraints"));

    // Check some statistics on the solution
    testCompare = Teuchos::rcp(new NOX::TestCompare(
				                 globalData->locaUtils->out(), 
						 *(globalData->locaUtils)));

    // Evaluate blocks
    grp->computeF();
    grp->computeJacobian();

    // A
    A = grp->getX().createMultiVector(nConstraints);
    A->random();
    
    // B
    constraints->setX(grp->getX());
    B = grp->getX().createMultiVector(nConstraints);
    B->random();
    constraints->setDgDx(*B);
    constraints->computeConstraints();    
    constraints->computeDX();

    // C
    C = Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(nConstraints,
								 nConstraints));
    C->random();

    // Set up left- and right-hand sides
    F = grp->getX().createMultiVector(nRHS);
    F->random();
    G = Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(nConstraints,
								 nRHS));
    G->random();
    X_bordering = F->clone(NOX::ShapeCopy);
    Y_bordering = 
      Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(nConstraints,
							       nRHS));
    X_householder = F->clone(NOX::ShapeCopy);
    Y_householder = 
      Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(nConstraints,
							       nRHS));

    string testName;

    // Test all nonzero
    testName = "Testing all nonzero";
    ierr += testSolve(false, false, false, false, false,
		      reltol, abstol, testName);

    // Test A = 0
    testName = "Testing A=0";
    ierr += testSolve(true, false, false, false, false,
		      reltol, abstol, testName);

    // Test B = 0
    testName = "Testing B=0";
    ierr += testSolve(false, true, false, false, false,
		      reltol, abstol, testName);

    // Test C = 0
    testName = "Testing C=0";
    ierr += testSolve(false, false, true, false, false,
		      reltol, abstol, testName);

    // Test F = 0
    testName = "Testing F=0";
    ierr += testSolve(false, false, false, true, false,
		      reltol, abstol, testName);

    // Test G = 0
    testName = "Testing G=0";
    ierr += testSolve(false, false, false, false, true,
		      reltol, abstol, testName);

    // Test A,B = 0
    testName = "Testing A,B=0";
    ierr += testSolve(true, true, false, false, false,
		      reltol, abstol, testName);

    // Test A,F = 0
    testName = "Testing A,F=0";
    ierr += testSolve(true, false, false, true, false,
		      reltol, abstol, testName);

    // Test A,G = 0
    testName = "Testing A,G=0";
    ierr += testSolve(true, false, false, false, true,
		      reltol, abstol, testName);

    // Test B,F = 0
    testName = "Testing B,F=0";
    ierr += testSolve(false, true, false, true, false,
		      reltol, abstol, testName);

    // Test B,G = 0
    testName = "Testing B,G=0";
    ierr += testSolve(false, true, false, false, true,
		      reltol, abstol, testName);

    // Test C,F = 0
    testName = "Testing C,F=0";
    ierr += testSolve(false, false, true, true, false,
		      reltol, abstol, testName);

    // Test C,G = 0
    testName = "Testing C,G=0";
    ierr += testSolve(false, false, true, false, true,
		      reltol, abstol, testName);

    // Test F,G = 0
    testName = "Testing F,G=0";
    ierr += testSolve(false, false, false, true, true,
		      reltol, abstol, testName);

    // Test A,B,F = 0
    testName = "Testing A,B,F=0";
    ierr += testSolve(true, true, false, true, false,
		      reltol, abstol, testName);

    // Test A,B,G = 0
    testName = "Testing A,B,G=0";
    ierr += testSolve(true, true, false, false, true,
		      reltol, abstol, testName);

    // Test A,F,G = 0
    testName = "Testing A,F,G=0";
    ierr += testSolve(true, false, false, true, true,
		      reltol, abstol, testName);

    // Test B,F,G = 0
    testName = "Testing B,F,G=0";
    ierr += testSolve(false, true, false, true, true,
		      reltol, abstol, testName);

    // Test C,F,G = 0
    testName = "Testing C,F,G=0";
    ierr += testSolve(false, false, true, true, true,
		      reltol, abstol, testName);

    // Test A,B,F,G = 0
    testName = "Testing A,B,F,G=0";
    ierr += testSolve(true, true, false, true, true,
		      reltol, abstol, testName);

    destroyGlobalData(globalData);

  }

  catch (std::exception& e) {
    cout << e.what() << endl;
    ierr = 1;
  }
  catch (const char *s) {
    cout << s << endl;
    ierr = 1;
  }
  catch (...) {
    cout << "Caught unknown exception!" << endl;
    ierr = 1;
  }

  if (MyPID == 0) {
    if (ierr == 0)
      cout << "All tests passed!" << endl;
    else
      cout << ierr << " test(s) failed!" << endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return ierr;
}
