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
#include "LOCA_LAPACK.H"

#include "LOCA_BorderedSystem_AbstractStrategy.H"
#include "LOCA_Parameter_SublistParser.H"

#include "ChanProblemInterface.H"
#include "NormConstraint.H"
#include "NOX_TestCompare.H"

// Global variables used in main() and testSolve()
Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup> grp;
Teuchos::RefCountPtr<NormConstraint> constraints;
Teuchos::RefCountPtr<LOCA::Parameter::SublistParser> parsedParams;
Teuchos::RefCountPtr<LOCA::BorderedSystem::AbstractStrategy> bordering;
Teuchos::RefCountPtr<LOCA::BorderedSystem::AbstractStrategy> direct;
Teuchos::RefCountPtr<LOCA::GlobalData> globalData;
Teuchos::RefCountPtr<NOX::TestCompare> testCompare;
Teuchos::RefCountPtr<NOX::Abstract::MultiVector> A;
Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> C;
Teuchos::RefCountPtr<NOX::Abstract::MultiVector> F;
Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> G;
Teuchos::RefCountPtr<NOX::Abstract::MultiVector> X_bordering;
Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> Y_bordering;
Teuchos::RefCountPtr<NOX::Abstract::MultiVector> X_direct;
Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> Y_direct;

int  
testSolve(bool flagA, bool flagB, bool flagC, bool flagF, bool flagG,
	  bool contiguous, double reltol, double abstol, 
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
  bordering->setIsContiguous(contiguous);
  bordering->setMatrixBlocks(grp, a, constraints, c);

  direct->setIsContiguous(contiguous);
  direct->setMatrixBlocks(grp, a, constraints, c);

  // Solve using bordering
  NOX::Abstract::Group::ReturnType borderingStatus = 
    bordering->applyInverse(*(parsedParams->getSublist("Linear Solver")),
			    f.get(), g.get(), *X_bordering, *Y_bordering);
  if (borderingStatus != NOX::Abstract::Group::Ok)
    ++ierr;

  // Solve using direct
  NOX::Abstract::Group::ReturnType directStatus = 
    direct->applyInverse(*(parsedParams->getSublist("Linear Solver")),
			 f.get(), g.get(), *X_direct, *Y_direct);
  if (directStatus != NOX::Abstract::Group::Ok)
    ++ierr;
  
  for (int i=0; i<Y_bordering->numCols(); i++) {
    stringstream sstr;
    sstr << "Column " << i;
    ierr += testCompare->testVector((*X_bordering)[i],
				    (*X_direct)[i], reltol, abstol,
				    sstr.str());
  }

  for (int i=0; i<Y_bordering->numRows(); i++) 
    for (int j=0; j<Y_bordering->numCols(); j++) {
      stringstream sstr;
      sstr << "Scalars entry (" << i << "," << j << ")";
      ierr += testCompare->testValue((*Y_bordering)(i,j),
				     (*Y_direct)(i,j), reltol,
				     sstr.str(),
				     NOX::TestCompare::Relative);
    }

  return ierr;
}
	  

int main(int argc, char *argv[])
{
  int n = 10;
  double alpha = 1.0;
  double beta = 0.0;
  double gamma = 2.0;
  double scale = 1.0;
  int ierr = 0;
  double reltol = 1.0e-13;
  double abstol = 1.0e-13;

  alpha = alpha / scale;

  try {

    bool verbose = false;
    // Check for verbose output
    if (argc>1)
      if (argv[1][0]=='-' && argv[1][1]=='v') 
	verbose = true;

    // Create parameter list
    Teuchos::RefCountPtr<NOX::Parameter::List> paramList = 
      Teuchos::rcp(new NOX::Parameter::List);

    // Create LOCA sublist
    NOX::Parameter::List& locaParamsList = paramList->sublist("LOCA");

    // Create the constraints list
    NOX::Parameter::List& constraintsList = 
      locaParamsList.sublist("Constraints");
    constraintsList.setParameter("Bordered Solver Method", "Bordering");

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    NOX::Parameter::List& nlParams = paramList->sublist("NOX");

    NOX::Parameter::List& nlPrintParams = nlParams.sublist("Printing");
    if (verbose)
       nlPrintParams.setParameter("Output Information", 
				  NOX::Utils::Error +
				  NOX::Utils::Details +
				  NOX::Utils::OuterIteration + 
				  NOX::Utils::InnerIteration + 
				  NOX::Utils::Warning +
				  NOX::Utils::TestDetails + 
				  NOX::Utils::StepperIteration +
				  NOX::Utils::StepperDetails);
     else
       nlPrintParams.setParameter("Output Information", NOX::Utils::Error);

    // Create LAPACK factory
    Teuchos::RefCountPtr<LOCA::Abstract::Factory> lapackFactory =
      Teuchos::rcp(new LOCA::LAPACK::Factory);

    // Create global data object
    globalData = LOCA::createGlobalData(paramList, lapackFactory);

    // Create parsed parameter list
    parsedParams = 
      Teuchos::rcp(new LOCA::Parameter::SublistParser(globalData));
    parsedParams->parseSublists(paramList);

    // Set up the problem interface
    ChanProblemInterface chan(globalData, n, alpha, beta, scale);
    LOCA::ParameterVector p;
    p.addParameter("alpha",alpha);
    p.addParameter("beta",beta);
    p.addParameter("gamma",gamma);
    p.addParameter("scale",scale);
  
    // Create a group which uses that problem interface. The group will
    // be initialized to contain the default initial guess for the
    // specified problem.
    grp = Teuchos::rcp(new LOCA::LAPACK::Group(globalData, chan));
    grp->setParams(p);

    // Change initial guess to a random vector
    Teuchos::RefCountPtr<NOX::Abstract::Vector> xnew = grp->getX().clone();
    xnew->random();
    grp->setX(*xnew);

    // Create the constraints object & constraint param IDs list
    constraints = Teuchos::rcp(new NormConstraint(n, p));
    Teuchos::RefCountPtr< vector<int> > constraintParamIDs = 
      Teuchos::rcp(new vector<int>(1));
    (*constraintParamIDs)[0] = p.getIndex("alpha");

    // Create bordering solver
    bordering
      = globalData->locaFactory->createBorderedSystemStrategy(
				     parsedParams, 
				     parsedParams->getSublist("Constraints"));

    // Change strategy to LAPACK Direct Solve
    constraintsList.setParameter("Bordered Solver Method", 
				 "LAPACK Direct Solve");

    // Create direct solver
    direct
      = globalData->locaFactory->createBorderedSystemStrategy(
				     parsedParams, 
				     parsedParams->getSublist("Constraints"));

    // Check some statistics on the solution
    testCompare = Teuchos::rcp(new NOX::TestCompare(
				                 globalData->locaUtils->out(), 
						 *(globalData->locaUtils)));

    // Evaluate blocks
    grp->computeF();
    grp->computeJacobian();

    // B
    constraints->setX(grp->getX());
    constraints->computeConstraints();
    constraints->computeDX();

    // C
    C = Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(1,1));
    C->random();

    // Set up left- and right-hand sides
    F = grp->getX().createMultiVector(3);
    F->random();
    G = Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(1,2));
    G->random();
    X_bordering = F->clone(3);
    Y_bordering = 
      Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(1,2));
    X_direct = F->clone(3);
    Y_direct = Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(1,2));

    // Setup A block as a view of the last column of F
    vector<int> indexA(1); 
    indexA[0] = 2;
    A = F->subView(indexA);

    string testName;

    // Test all nonzero, contiguous
    testName = "Testing all nonzero, contiguous";
    ierr += testSolve(false, false, false, false, false, true,
		      reltol, abstol, testName);

    // Test B = 0, contiguous
    testName = "Testing B=0, contiguous";
    ierr += testSolve(false, true, false, false, false, true,
		      reltol, abstol, testName);

    // Test C = 0, contiguous
    testName = "Testing C=0, contiguous";
    ierr += testSolve(false, false, true, false, false, true,
		      reltol, abstol, testName);

    // Test G = 0, contiguous
    testName = "Testing G=0, contiguous";
    ierr += testSolve(false, false, false, false, true, true,
		      reltol, abstol, testName);

    // Test B,G = 0, contiguous
    testName = "Testing B,G=0, contiguous";
    ierr += testSolve(false, true, false, false, true, true,
		      reltol, abstol, testName);

    // Test C,G = 0, contiguous
    testName = "Testing C,G=0, contiguous";
    ierr += testSolve(false, false, true, false, true, true,
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

  if (ierr == 0)
    cout << "All tests passed!" << endl;
  else
    cout << ierr << " test(s) failed!" << endl;

  return ierr;
}
