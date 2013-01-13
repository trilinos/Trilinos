// $Id$ 
// $Source$ 

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

#include "LOCA.H"
#include "LOCA_LAPACK.H"

#include "LOCA_BorderedSolver_AbstractStrategy.H"
#include "LOCA_Parameter_SublistParser.H"
#include "LOCA_BorderedSolver_JacobianOperator.H"

#include "ChanProblemInterface.H"
#include "LOCALinearConstraint.H"
#include "NOX_TestCompare.H"

// Global variables used in main() and testSolve()
Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> grp;
Teuchos::RCP<LOCA::BorderedSolver::JacobianOperator> op;
Teuchos::RCP<LinearConstraint> constraints;
Teuchos::RCP<LOCA::Parameter::SublistParser> parsedParams;
Teuchos::RCP<LOCA::BorderedSolver::AbstractStrategy> bordering;
Teuchos::RCP<LOCA::BorderedSolver::AbstractStrategy> direct;
Teuchos::RCP<LOCA::GlobalData> globalData;
Teuchos::RCP<NOX::TestCompare> testCompare;
Teuchos::RCP<NOX::Abstract::MultiVector> A;
Teuchos::RCP<NOX::Abstract::MultiVector> B;
Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> C;
Teuchos::RCP<NOX::Abstract::MultiVector> F;
Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> G;
Teuchos::RCP<NOX::Abstract::MultiVector> X_bordering;
Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> Y_bordering;
Teuchos::RCP<NOX::Abstract::MultiVector> X_direct;
Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> Y_direct;

int  
testSolve(bool flagA, bool flagB, bool flagC, bool flagF, bool flagG,
	  double reltol, double abstol, const std::string& testName) {
  int ierr = 0;

  if (globalData->locaUtils->isPrintType(NOX::Utils::TestDetails))
    globalData->locaUtils->out() 
      << std::endl << "***** " << testName << " *****" << std::endl;

  Teuchos::RCP<NOX::Abstract::MultiVector> a = 
    Teuchos::null;
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> c = 
    Teuchos::null;
  Teuchos::RCP<NOX::Abstract::MultiVector> f = 
    Teuchos::null;
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> g = 
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
  bordering->setMatrixBlocks(op, a, constraints, c);
  bordering->initForSolve();

  direct->setMatrixBlocks(op, a, constraints, c);
  direct->initForSolve();

  X_bordering->init(0.0);
  Y_bordering->putScalar(0.0);
  X_direct->init(0.0);
  Y_direct->putScalar(0.0);

  // Solve using bordering
  NOX::Abstract::Group::ReturnType borderingStatus = 
    bordering->applyInverseTranspose(*(parsedParams->getSublist("Linear Solver")),
			    f.get(), g.get(), *X_bordering, *Y_bordering);
  if (borderingStatus != NOX::Abstract::Group::Ok)
    ++ierr;

  // Solve using direct
  NOX::Abstract::Group::ReturnType directStatus = 
    direct->applyInverseTranspose(*(parsedParams->getSublist("Linear Solver")),
			 f.get(), g.get(), *X_direct, *Y_direct);
  if (directStatus != NOX::Abstract::Group::Ok)
    ++ierr;
  
  for (int i=0; i<Y_bordering->numCols(); i++) {
    std::stringstream sstr;
    sstr << "Column " << i;
    ierr += testCompare->testVector((*X_bordering)[i],
				    (*X_direct)[i], reltol, abstol,
				    sstr.str());
  }

  for (int i=0; i<Y_bordering->numRows(); i++) 
    for (int j=0; j<Y_bordering->numCols(); j++) {
      std::stringstream sstr;
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
  int nConstraints = 10;
  int nRHS = 7;

  int n = 100;
  double alpha = 1.0;
  double beta = 0.0;
  double gamma = 2.0;
  double scale = 1.0;
  int ierr = 0;
  double reltol = 1.0e-9;
  double abstol = 1.0e-9;

  alpha = alpha / scale;

  try {

    bool verbose = false;
    // Check for verbose output
    if (argc>1)
      if (argv[1][0]=='-' && argv[1][1]=='v') 
	verbose = true;

    // Create parameter list
    Teuchos::RCP<Teuchos::ParameterList> paramList = 
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

    // Create LAPACK factory
    Teuchos::RCP<LOCA::Abstract::Factory> lapackFactory =
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

    // Create Jacobian operator
    op = Teuchos::rcp(new LOCA::BorderedSolver::JacobianOperator(grp));

    // Change initial guess to a random vector
    Teuchos::RCP<NOX::Abstract::Vector> xnew = grp->getX().clone();
    xnew->random();
    grp->setX(*xnew);

    // Create the constraints object & constraint param IDs list
    constraints = Teuchos::rcp(new LinearConstraint(nConstraints, p, *xnew));
    Teuchos::RCP< std::vector<int> > constraintParamIDs = 
      Teuchos::rcp(new std::vector<int>(1));
    (*constraintParamIDs)[0] = p.getIndex("alpha");

    // Create bordering solver
    bordering
      = globalData->locaFactory->createBorderedSolverStrategy(
				     parsedParams, 
				     parsedParams->getSublist("Constraints"));

    // Change strategy to LAPACK Direct Solve
    constraintsList.set("Bordered Solver Method", 
				 "LAPACK Direct Solve");

    // Create direct solver
    direct
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
    C = 
      Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(nConstraints,
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
    X_direct = F->clone(NOX::ShapeCopy);
    Y_direct = 
      Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(nConstraints,
							       nRHS));

    std::string testName;

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

    LOCA::destroyGlobalData(globalData);
  }

  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
    ierr = 1;
  }
  catch (const char *s) {
    std::cout << s << std::endl;
    ierr = 1;
  }
  catch (...) {
    std::cout << "Caught unknown exception!" << std::endl;
    ierr = 1;
  }

  if (ierr == 0)
    std::cout << "All tests passed!" << std::endl;
  else
    std::cout << ierr << " test(s) failed!" << std::endl;

  return ierr;
}
