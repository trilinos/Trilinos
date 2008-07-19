// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
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

/*!
  \file StatusTests.C
  
  This is a test for the nox status tests and the
  NOX::StatusTest::Factory object.  It relies on the Broyden problem
  written by Brett Bader described below.

  This test problem is a modified extension of the "Broyden
  Tridiagonal Problem" from Jorge J. More', Burton S. Garbow, and
  Kenneth E. Hillstrom, Testing Unconstrained Optimization Software,
  ACM TOMS, Vol. 7, No. 1, March 1981, pp. 14-41.  The modification
  involves squaring the last equation fn(x) and using it in a
  homotopy-type equation.

  The parameter "lambda" is a homotopy-type parameter that may be
  varied from 0 to 1 to adjust the ill-conditioning of the problem.  A
  value of 0 is the original, unmodified problem, while a value of 1
  is that problem with the last equation squared.  Typical values for
  increasingly ill-conditioned problems might be 0.9, 0.99, 0.999,
  etc.

  The standard starting point is x(i) = -1, but setting x(i) = 0 tests
  the selected global strategy.

  \authors Roger Pawlowski (SNL 1416)
*/


#include "NOX.H"
#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "NOX_LAPACK_Group.H"
#include "NOX_TestUtils.H"
#include "Teuchos_ScalarTraits.hpp"
#ifdef HAVE_TEUCHOS_EXTENDED
#include "Teuchos_XMLParameterListHelpers.hpp"
#endif

//! Interface to modified Broyden problem defined in Broyden.C
class Broyden : public NOX::LAPACK::Interface {

public:
 
  //! Constructor
  Broyden(int m, double lambdaVal=0) : 
    initialGuess(m),
    solution(m)
  {
    n = m;
    lambda = lambdaVal;

    cout << "Broyden ill-conditioning: lambda = " << lambda << "\n"; 
    
    for (int i=0; i<n; i++) {
      // initialGuess(i) = -100;   // Case for lambdaBar != 1.0
      initialGuess(i) = 0;      // General testing
      // initialGuess(i) = -1;     // Standard
      solution(i) = 1;
    }
    fevals = 0;
  };

  //! Destructor
  ~Broyden() { cout << "Function evaluations: " << fevals << "\n"; };

  const NOX::LAPACK::Vector& getInitialGuess()
  {
    return initialGuess;
  };

  //! Return true solution vector
  const NOX::LAPACK::Vector& getSolution()
  {
    return solution;
  };

  //! Sets the initial guess to NaN so that we can test the FiniteValues test.
  void setNaNSolution()
  {
    for (int i=0; i < 100; ++i)
      initialGuess(i) = 1.0 / Teuchos::ScalarTraits<double>::zero();
  };

  bool computeF(NOX::LAPACK::Vector& f, const NOX::LAPACK::Vector &x)
  {
    double fn;
    
    f(0) = (3-2*x(0))*x(0) - 2*x(1) + 1;
    for (int i=1; i<n-1; i++)
      f(i) = (3-2*x(i))*x(i) - x(i-1) - 2*x(i+1) + 1;

    fn = (3-2*x(n-1))*x(n-1) - x(n-2) + 1;
    f(n-1) = (1-lambda)*fn + lambda*(fn*fn);

    fevals++;
    return true;
  };
  
  bool computeJacobian(NOX::LAPACK::Matrix<double>& J, 
		       const NOX::LAPACK::Vector & x)
  {
    double fn;
    double dfndxn;
    
    // F(0) = (3-2*x(0))*x(0) - 2*x(1) + 1;
    J(0,0) = 3 - 4*x(0);
    J(0,1) = -2;

    // F(i) = (3-2*x(i))*x(i) - x(i-1) - 2*x(i+1) + 1;
    for (int i=1; i<n-1; i++) {
      J(i,i-1) = -1;
      J(i,i)   = 3 - 4*x(i);
      J(i,i+1) = -2;
    }

    // F(n-1) = ((3-2*x(n-1))*x(n-1) - x(n-2) + 1)^2;
    fn = (3-2*x(n-1))*x(n-1) - x(n-2) + 1;
    dfndxn = 3-4*x(n-1);
    J(n-1,n-1) = (1-lambda)*(dfndxn) + lambda*(2*dfndxn*fn);
    J(n-1,n-2) = (1-lambda)*(-1) + lambda*(-2*fn);

    return true;
  };

private:
  
  //! Problem size
  int n;
  //! Number of calls to computeF
  int fevals;
  //! Ill-conditioning parameters
  double lambda;
  //! Initial guess
  NOX::LAPACK::Vector initialGuess;
  //! Correct answer
  NOX::LAPACK::Vector solution;

};

class MyTest : public NOX::StatusTest::Generic {

public:

  MyTest(double tol) : status(NOX::StatusTest::Unconverged), tolerance(tol) {}

  ~MyTest() {}

  NOX::StatusTest::StatusType 
  checkStatus(const NOX::Solver::Generic &problem, 
	      NOX::StatusTest::CheckType checkType) 
  {
    
    if (!problem.getSolutionGroup().isF()) {
      std::string msg = "Error - MyTest::checkStatus() - The residual F is not up-to-date with the solution.  Please fix your algorithm!";
      TEST_FOR_EXCEPTION(true, std::logic_error, msg);
    }
    

    switch (checkType) {
    case NOX::StatusTest::Complete:
    case NOX::StatusTest::Minimal:
      norm_f = problem.getSolutionGroup().getNormF();
      status = (norm_f < tolerance) ? NOX::StatusTest::Converged : NOX::StatusTest::Unconverged;
      break;      
    case NOX::StatusTest::None:
    default:
      norm_f = -1.0;
      status = NOX::StatusTest::Unevaluated;
      break;
    }
    
    return status;
  }
  
  NOX::StatusTest::StatusType getStatus() const 
  { return status; }

  std::ostream& print(std::ostream &stream, int indent) const 
  {
    for (int j = 0; j < indent; j ++)
      stream << ' ';
    stream << status;
    stream << "MyTest = " << norm_f << " < " << tolerance;
    stream << endl;
    return stream;
  }
    
private:
  NOX::StatusTest::StatusType status;
  double tolerance;
  double norm_f;
};


//! Main subroutine of Broyden.C
int main(int argc, char *argv[])
{
  // Basic declarations to clean up the code
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using NOX::StatusTest::Combo;
  using NOX::StatusTest::FiniteValue;
  using NOX::StatusTest::MaxIters;
  using NOX::StatusTest::Divergence;
  using NOX::StatusTest::Stagnation;
  using namespace NOX::StatusTest;

  cout << "Started" << endl;

  int final_status_value = 0; // zero = success, !0 = failed

  // Create the top level parameter list
  RCP<Teuchos::ParameterList> solverParametersPtr =
    rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& solverParameters = *solverParametersPtr;

  // Set the nonlinear solver method
  //solverParameters.set("Nonlinear Solver", "Tensor-Krylov Based");
  //solverParameters.set("Nonlinear Solver", "Tensor Based");
  solverParameters.set("Nonlinear Solver", "Line Search Based");
  
  // Sublist for printing parameters
  Teuchos::ParameterList& printParams = solverParameters.sublist("Printing");
  //printParams.set("MyPID", 0); 
  printParams.set("Output Precision", 3);
  printParams.set("Output Processor", 0);
  printParams.set("Output Information", 
			NOX::Utils::OuterIteration + 
			NOX::Utils::OuterIterationStatusTest + 
			NOX::Utils::InnerIteration +
			NOX::Utils::Details + 
			NOX::Utils::Warning);
  NOX::Utils utils(printParams);


  // Convergence tests and factory

  cout << "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
  cout << "Testing Convergence tests (NormF, NormUpdate, NormWRMS) ..."
       << endl;
  cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;

  // Set up the problem interface
  Broyden broyden(100,0.99);
  
  // Create a group which uses that problem interface. The group will
  // be initialized to contain the default initial guess for the
  // specified problem.
  RCP<NOX::LAPACK::Group> grp = rcp(new NOX::LAPACK::Group(broyden));

  // **************************************
  // Create the convergence tests
  // **************************************
  Teuchos::ParameterList stl;
  stl.set("Test Type", "Combo");
  stl.set("Combo Type", "OR");
  stl.set("Number of Tests", 5);
  Teuchos::ParameterList& conv = stl.sublist("Test 0");
  Teuchos::ParameterList& fv = stl.sublist("Test 1");
  Teuchos::ParameterList& divergence = stl.sublist("Test 2");
  Teuchos::ParameterList& stagnation = stl.sublist("Test 3");
  Teuchos::ParameterList& maxiters = stl.sublist("Test 4");

  conv.set("Test Type", "Combo");
  conv.set("Combo Type", "AND");
  conv.set("Number of Tests", 3);
  Teuchos::ParameterList& normF = conv.sublist("Test 0");
  Teuchos::ParameterList& normWRMS = conv.sublist("Test 1");
  Teuchos::ParameterList& normUpdate = conv.sublist("Test 2");
  normF.set("Test Type", "NormF");
  normF.set("Tolerance", 1.0e-12);
  normF.set("Norm Type", "Two Norm");
  normF.set("Scale Type", "Unscaled");
  normWRMS.set("Test Type", "NormWRMS");
  normWRMS.set("Absolute Tolerance", 1.0e-8);
  normWRMS.set("Relative Tolerance", 1.0e-5);
  normWRMS.set("Tolerance", 1.0);
  normWRMS.set("BDF Multiplier", 1.0);
  normWRMS.set("Alpha", 1.0);
  normWRMS.set("Beta", 0.5);
  normUpdate.set("Test Type", "NormUpdate");
  normUpdate.set("Norm Type", "One Norm");
  normUpdate.set("Scale Type", "Scaled");

  fv.set("Test Type", "FiniteValue");
  fv.set("Vector Type", "F Vector");
  fv.set("Norm Type", "Two Norm");

  divergence.set("Test Type", "Divergence");
  divergence.set("Tolerance", 1.0e+20);
  divergence.set("Consecutive Iterations", 3);
  
  stagnation.set("Test Type", "Stagnation");
  stagnation.set("Tolerance", 1.0);
  stagnation.set("Consecutive Iterations", 5);
  
  maxiters.set("Test Type", "MaxIters");
  maxiters.set("Maximum Iterations", 20);

  Teuchos::RCP<NOX::StatusTest::Generic> statusTestsCombo;
  Teuchos::RCP<Teuchos::ParameterList> st_params;

#ifdef HAVE_TEUCHOS_EXTENDED
  cout << "Writing parameter list to \"input.xml\"" << endl;
  Teuchos::writeParameterListToXmlFile(stl, "input.xml");
  cout << "Reading parameter list from \"input.xml\"" << endl;
  statusTestsCombo = NOX::StatusTest::buildStatusTests("input.xml", utils);
#else
  statusTestsCombo = NOX::StatusTest::buildStatusTests(stl, utils);
#endif

  // **************************************
  // Finished: Create the convergence tests
  // **************************************
  
  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver = 
    NOX::Solver::buildSolver(grp, statusTestsCombo, solverParametersPtr);

  // Solve the nonlinear system
  NOX::StatusTest::StatusType status = solver->solve();

  // Print final status
  if (status == NOX::StatusTest::Converged && 
      solver->getNumIterations() == 12) {
    final_status_value += 0;
    cout << "\nConvergence tests passed!" << endl;
  }
  else {
    final_status_value += 1;
    cout << "\nConvergence tests Failed!" << endl;
  }

  // Re-run test with complete checks of status tests
  {
    cout << "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "re-runningConvergence tests (NormF, NormUpdate, NormWRMS)\n"
	 << "with NOX::Solver::CheckType = Complete"
	 << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    
    Broyden interface(100,0.99);
    RCP<NOX::LAPACK::Group> group = rcp(new NOX::LAPACK::Group(interface));
 
    Teuchos::RCP<Teuchos::ParameterList> tmpParams =
      Teuchos::rcp(new Teuchos::ParameterList);
    
    tmpParams = solverParametersPtr;

    tmpParams->sublist("Solver Options").
      set("Status Test Check Type", "Complete");
    
    Teuchos::RCP<NOX::Solver::Generic> solver = 
      NOX::Solver::buildSolver(group, statusTestsCombo, tmpParams);
  
    status = solver->solve();

    cout << *tmpParams << endl;

    if (status == NOX::StatusTest::Converged && 
	solver->getNumIterations() == 12) {
      final_status_value += 0;
      cout << "Convergence (Complete) tests passed!\n" << endl;
    }
    else {
      final_status_value += 1;
      cout << "Convergence (Complete)tests Failed!\n" << endl;
    }

  }

  // Test the user defined status test
  {
    cout << "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "Testing User Defined Status Test" << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" << endl;
    Broyden interface(100,0.99);
    RCP<NOX::LAPACK::Group> group = rcp(new NOX::LAPACK::Group(interface));
    
    ParameterList tmp_list;
    tmp_list.set("Test Type", "Combo");
    tmp_list.set("Combo Type", "OR");
    tmp_list.set("Number of Tests", 2);
    ParameterList& max_iters = tmp_list.sublist("Test 0");
    ParameterList& user_defined = tmp_list.sublist("Test 1");

    max_iters.set("Test Type", "MaxIters");
    max_iters.set("Maximum Iterations", 20);

    user_defined.set("Test Type", "User Defined");
    Teuchos::RCP<NOX::StatusTest::Generic> myTest = 
      Teuchos::rcp(new MyTest(1.0e-3));
    user_defined.set("User Status Test", myTest);

    std::cout << tmp_list << endl;

    RCP<NOX::StatusTest::Generic> combo = buildStatusTests(tmp_list, utils);

    Teuchos::RCP<NOX::Solver::Generic> solver = 
      NOX::Solver::buildSolver(group, combo, solverParametersPtr);
  
    status = solver->solve();

    if (status == NOX::StatusTest::Converged && 
	solver->getNumIterations() == 11) {
      final_status_value += 0;
      cout << "\nUser Defined test passed!\n" << endl;
    }
    else {
      final_status_value += 1;
      cout << "\nUser Defined test failed!\n" << endl;
    }
  }


  // Tagging
  {
    cout << "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "Testing NOX::StatusTest::Factory Tagging Option" << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" << endl;

    std::map< std::string, Teuchos::RCP<NOX::StatusTest::Generic> > my_map;
    normF.set("Tag", "Norm F Status Test");

    Teuchos::RCP<NOX::StatusTest::Generic> all_tests = 
      NOX::StatusTest::buildStatusTests(stl, utils, &my_map);
    
    Teuchos::RCP<NOX::StatusTest::Generic> my_test = 
      my_map["Norm F Status Test"];

    Teuchos::RCP<NOX::StatusTest::NormF> my_normF_test= 
      Teuchos::rcp_dynamic_cast<NOX::StatusTest::NormF>(my_test);

    my_normF_test->print(std::cout);

    if (!Teuchos::is_null(my_normF_test)) {
      final_status_value += 0;
      cout << "\nTagging test passed!" << endl;
    }
    else {
      final_status_value += 1;
      cout << "\nTagging test failed!" << endl;
    }
  }

  // NormF RELATIVE tolerance
  {
    
    cout << "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "Testing NOX::StatusTest::NormF RELATIVE Test" << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" << endl;

    Teuchos::ParameterList p;

    p.set("Test Type", "Combo");
    p.set("Combo Type", "OR");
    p.set("Number of Tests", 2);

    Teuchos::ParameterList& conv = p.sublist("Test 0");
    conv.set("Test Type", "Combo");
    conv.set("Combo Type", "AND");
    conv.set("Number of Tests", 3);    

    Teuchos::ParameterList& normF_rel = conv.sublist("Test 0");
    normF_rel.set("Test Type", "NormF");
    normF_rel.set("Tolerance", 1.0e-2);
    normF_rel.set("Norm Type", "Two Norm");
    normF_rel.set("Scale Type", "Unscaled");
    normF_rel.set("Tag", "Rel Test");
    
    Broyden interface(100,0.99);
    RCP<NOX::LAPACK::Group> group = rcp(new NOX::LAPACK::Group(interface));
    Teuchos::RCP<NOX::Abstract::Group> ig = group;

    normF_rel.set("Initial Guess", ig);

    Teuchos::ParameterList& normF_rel2 = conv.sublist("Test 1");
    normF_rel2.set("Test Type", "NormF");
    normF_rel2.set("Tolerance", 1.0e-2);
    normF_rel2.set("Norm Type", "Two Norm");
    normF_rel2.set("Scale Type", "Unscaled");
    normF_rel2.set("Tag", "Rel Test 2");
    normF_rel2.set("Initial Guess", ig);

    Teuchos::ParameterList& normF_abs = conv.sublist("Test 2");
    normF_abs.set("Test Type", "NormF");
    normF_abs.set("Tolerance", 1.0e-2);
    normF_abs.set("Norm Type", "Two Norm");
    normF_abs.set("Scale Type", "Unscaled");
    normF_abs.set("Tag", "Abs Test");

    Teuchos::ParameterList& max_iters = p.sublist("Test 1");
    max_iters.set("Test Type", "MaxIters");
    max_iters.set("Maximum Iterations", 15);

    std::map< std::string, Teuchos::RCP<NOX::StatusTest::Generic> > tag_map;
    statusTestsCombo = NOX::StatusTest::buildStatusTests(p, utils, &tag_map);
    
    RCP<Teuchos::ParameterList> sp = rcp(new Teuchos::ParameterList);
    sp->set("Nonlinear Solver", "Line Search Based");
    sp->sublist("Printing") = printParams;
    sp->sublist("Solver Options").set("Status Test Check Type", "Complete");

    Teuchos::RCP<NOX::Solver::Generic> solver = 
      NOX::Solver::buildSolver(ig, statusTestsCombo, sp);

    NOX::StatusTest::StatusType status = solver->solve();
    
    if (status == NOX::StatusTest::Converged && 
	solver->getNumIterations() == 10) {
      final_status_value += 0;
      cout << "\nNormF RELATIVE test 1/2 passed!\n" << endl;
    }
    else {
      final_status_value += 1;
      cout << "\nNormF RELATIVE test 1/2 failed!\n" << endl;
    }

    // reset absolute and relative tolerance and keep going
    Teuchos::RCP<NOX::StatusTest::NormF> abs_test = 
      Teuchos::rcp_dynamic_cast<NOX::StatusTest::NormF>(tag_map["Abs Test"]);
    Teuchos::RCP<NOX::StatusTest::NormF> rel_test = 
      Teuchos::rcp_dynamic_cast<NOX::StatusTest::NormF>(tag_map["Rel Test"]);
    Teuchos::RCP<NOX::StatusTest::NormF> rel_test2 = 
      Teuchos::rcp_dynamic_cast<NOX::StatusTest::NormF>(tag_map["Rel Test 2"]);
    NOX::Abstract::Group& nonconst_solution = 
      const_cast<NOX::Abstract::Group&>(solver->getSolutionGroup());

    abs_test->reset(1.0e-8);
    rel_test->reset(nonconst_solution, 1.0e-4);
    rel_test2->reset(1.0e-4);

    solver->reset(solver->getSolutionGroup().getX());
    status = solver->solve();
    
    if (status == NOX::StatusTest::Converged && 
	solver->getNumIterations() == 2) {
      final_status_value += 0;
      cout << "\nNormF RELATIVE test 2/2 passed!\n" << endl;
    }
    else {
      final_status_value += 1;
      cout << "\nNormF RELATIVE test 2/2 failed!\n" << endl;
    }

  }

  // **********************
  // Individual status test options - for failure tests
  // **********************

  // Max Iters
  {
    cout << "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "Testing NOX::StatusTest::MaxIters" << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" << endl;
    Broyden interface(100,0.99);
    RCP<NOX::LAPACK::Group> group = rcp(new NOX::LAPACK::Group(interface));
    ParameterList p;
    p.set("Test Type", "MaxIters");
    p.set("Maximum Iterations", 5);
    Teuchos::RCP<NOX::Solver::Generic> solver = 
      NOX::Solver::buildSolver(group, 
			       NOX::StatusTest::buildStatusTests(p, utils), 
			       solverParametersPtr);
    status = solver->solve();
    
    // A failure reported by max iters is a passing test
    if (status == NOX::StatusTest::Failed) {
      final_status_value += 0;
      cout << "\nMaxIters test passed!\n" << endl;
    }
    else {
      final_status_value += 1;
      cout << "\nMaxIters test failed!\n" << endl;
    }
  }

  // Finite Value
  {
    cout << "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "Testing NOX::StatusTest::FiniteValue" << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" << endl;
    Broyden interface(100,1000.0);
    interface.setNaNSolution();
    RCP<NOX::LAPACK::Group> group = rcp(new NOX::LAPACK::Group(interface));
    
    RCP<Combo> combo = rcp(new Combo(Combo::OR));
    RCP<FiniteValue> fvst = rcp(new FiniteValue);
    RCP<MaxIters> mist = rcp(new MaxIters(20));
    combo->addStatusTest(fvst);
    combo->addStatusTest(mist);

    Teuchos::RCP<NOX::Solver::Generic> solver = 
      NOX::Solver::buildSolver(group, combo, solverParametersPtr);
    status = solver->solve();

    // A failure reported by finite value is a passing test
    if (status == NOX::StatusTest::Failed && 
	fvst->getStatus() == NOX::StatusTest::Failed) {
      final_status_value += 0;
      cout << "\nFiniteValue test passed!\n" << endl;
    }
    else {
      final_status_value += 1;
      cout << "\nFiniteValue test failed!\n" << endl;
    }
  }

  // Divergence
  {
    cout << "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "Testing NOX::StatusTest::Divergence" << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" << endl;
    Broyden interface(100,0.99);
    Teuchos::RCP<NOX::LAPACK::Group> group = 
      Teuchos::rcp(new NOX::LAPACK::Group(interface));

    RCP<Combo> combo = rcp(new Combo(Combo::OR));
    RCP<Divergence> divst = rcp(new Divergence(1.0e-5, 2));
    RCP<MaxIters> mist = rcp(new MaxIters(20));
    combo->addStatusTest(divst);
    combo->addStatusTest(mist);

    Teuchos::RCP<NOX::Solver::Generic> solver = 
      NOX::Solver::buildSolver(group, combo, solverParametersPtr);
    status = solver->solve();
    
    // A failure reported by divergence is a passing test
    if (status == NOX::StatusTest::Failed &&
	divst->getStatus() == NOX::StatusTest::Failed) {
      final_status_value += 0;
      cout << "\nDivergence test passed!\n" << endl;
    }
    else {
      final_status_value += 1;
      cout << "\nDivergence test failed!\n" << endl;
    }

  }

  // Stagnation
  {
    cout << "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "Testing NOX::StatusTest::Stagnation" << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" << endl;
    Broyden interface(100,0.99);
    Teuchos::RCP<NOX::LAPACK::Group> group = 
      Teuchos::rcp(new NOX::LAPACK::Group(interface));

    RCP<Combo> combo = rcp(new Combo(Combo::OR));
    RCP<Stagnation> stagst = rcp(new Stagnation(2, 1.0e-2));
    RCP<MaxIters> mist = rcp(new MaxIters(20));
    combo->addStatusTest(stagst);
    combo->addStatusTest(mist);
    
    Teuchos::RCP<NOX::Solver::Generic> solver = 
      NOX::Solver::buildSolver(group, combo, solverParametersPtr);
    status = solver->solve();
    
    // A failure reported by stagnation is a passing test
    if (status == NOX::StatusTest::Failed &&
	stagst->getStatus() == NOX::StatusTest::Failed) {
      final_status_value += 0;
      cout << "\nStagnation test passed!\n" << endl;
    }
    else {
      final_status_value += 1;
      cout << "\nStagnation test failed!\n" << endl;
    }
  }

  // **********************
  // Finished: Individual status test options
  // **********************

  cout << endl;
  if (final_status_value == 0)
    cout << "Test passed!" << endl;    
  else 
    cout << "Test failed!" << endl;
  cout << endl;

  return final_status_value;
}
