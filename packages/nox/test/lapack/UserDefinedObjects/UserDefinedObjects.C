// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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

/*!
  \file UserDefinedObjects.C
  
  This is a test for the "User Defined" direction and "User Defined"
  line searches.  It relies on the Broyden problem written by Brett
  Bader described below.

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
#ifdef HAVE_TEUCHOS_EXTENDED
#include "Teuchos_XMLParameterListHelpers.hpp"
#endif

//! These wil represent the user defined objects
#include "NOX_Direction_Newton.H"
#include "NOX_LineSearch_Polynomial.H"

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

    std::cout << "Broyden ill-conditioning: lambda = " << lambda << "\n"; 
    
    for (int i=0; i<n; i++) {
      // initialGuess(i) = -100;   // Case for lambdaBar != 1.0
      initialGuess(i) = 0;      // General testing
      // initialGuess(i) = -1;     // Standard
      solution(i) = 1;
    }
    fevals = 0;
  }

  //! Destructor
  ~Broyden() { std::cout << "Function evaluations: " << fevals << "\n"; }

  const NOX::LAPACK::Vector& getInitialGuess()
  {
    return initialGuess;
  }

  //! Return true solution vector
  const NOX::LAPACK::Vector& getSolution()
  {
    return solution;
  }

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
  }
  
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
  }

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

//! Main subroutine of Broyden.C
int main(int argc, char *argv[])
{
  // Basic declarations to clean up the code
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using namespace NOX::StatusTest;

  std::cout << "Started" << std::endl;

  int final_status_value = 0; // zero = success, !0 = failed

  RCP<ParameterList> solverParametersPtr = rcp(new ParameterList);
  ParameterList& solverParameters = *solverParametersPtr;
  solverParameters.set("Nonlinear Solver", "Line Search Based");
  
  // Create a user defined direction using the template builder
  ParameterList& dl = solverParametersPtr->sublist("Direction");
  dl.set("Method", "User Defined");
  RCP<NOX::Direction::UserDefinedFactory> uddf =
    rcp(new NOX::Direction::UserDefinedFactoryT<NOX::Direction::Newton>);
  dl.set("User Defined Direction Factory", uddf);
  
  // Create a user defined line search using the template builder
  ParameterList& lsl = solverParametersPtr->sublist("Line Search");
  lsl.set("Method", "User Defined");
  RCP<NOX::LineSearch::UserDefinedFactory> udlsf =
    rcp(new NOX::LineSearch::UserDefinedFactoryT<NOX::LineSearch::Polynomial>);
  lsl.set("User Defined Line Search Factory", udlsf);

  ParameterList& printParams = solverParameters.sublist("Printing");
  printParams.set("Output Precision", 3);
  printParams.set("Output Processor", 0);
  printParams.set("Output Information", 
		  NOX::Utils::OuterIteration + 
		  NOX::Utils::OuterIterationStatusTest + 
		  NOX::Utils::InnerIteration +
		  NOX::Utils::Details + 
		  NOX::Utils::Warning);

  Teuchos::ParameterList stl;
  stl.set("Test Type", "Combo");
  stl.set("Combo Type", "OR");
  stl.set("Number of Tests", 3);
  Teuchos::ParameterList& conv = stl.sublist("Test 0");
  Teuchos::ParameterList& fv = stl.sublist("Test 1");
  Teuchos::ParameterList& maxiters = stl.sublist("Test 2");

  conv.set("Test Type", "Combo");
  conv.set("Combo Type", "AND");
  conv.set("Number of Tests", 2);
  Teuchos::ParameterList& normF = conv.sublist("Test 0");
  Teuchos::ParameterList& normWRMS = conv.sublist("Test 1");
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

  fv.set("Test Type", "FiniteValue");
  fv.set("Vector Type", "F Vector");
  fv.set("Norm Type", "Two Norm");

  maxiters.set("Test Type", "MaxIters");
  maxiters.set("Maximum Iterations", 20);

  NOX::Utils utils(printParams);
  RCP<NOX::StatusTest::Generic> status_tests = buildStatusTests(stl, utils);

  Broyden broyden(100, 0.99);
  RCP<NOX::LAPACK::Group> grp = rcp(new NOX::LAPACK::Group(broyden));

  RCP<NOX::Solver::Generic> solver = 
    NOX::Solver::buildSolver(grp, status_tests, solverParametersPtr);

  NOX::StatusTest::StatusType status = solver->solve();

  std::cout << *solverParametersPtr << std::endl;

  if (status != NOX::StatusTest::Converged || 
      solver->getNumIterations() != 12) {
    final_status_value += 1;
    std::cout << "\nTest failed!\n" << std::endl;
  }
  else 
    std::cout << "\nTest passed!\n" << std::endl;    
    
  return final_status_value;
}
