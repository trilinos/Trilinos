#ifdef WITH_PRERELEASE
//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

/*!
  \file Broyden.C
  
   This is an example of using NOX with the NOX::Solver::TensorBased 
   tensor-Krylov method.

   This test problem is a modified extension of the "Broyden
   Tridiagonal Problem" from Jorge J. More', Burton S. Garbow, and
   Kenneth E. Hillstrom, Testing Unconstrained Optimization Software,
   ACM TOMS, Vol. 7, No. 1, March 1981, pp. 14-41.  The modification
   involves squaring the last equation fn(x) and using it in a
   homotopy-type equation.

   The parameter "lambda" is a homotopy-type parameter that may be
   varied from 0 to 1 to adjust the ill-conditioning of the problem.
   A value of 0 is the original, unmodified problem, while a value of
   1 is that problem with the last equation squared.  Typical values
   for increasingly ill-conditioned problems might be 0.9, 0.99,
   0.999, etc.

   The standard starting point is x(i) = -1, but x(i) = 0 requires a 
   few linesearches to test that code functionality.

   \author Brett Bader, UC Boulder, 2002
*/


#include "NOX.H"
#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "NOX_Example_Group.H"
#include "NOX_Solver_TensorBased.H"  

//! Interface to modified Broyden problem defined in Broyden.C
class Broyden : public NOX::Example::Interface {

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
      initialGuess(i) = 0;  // -1;
      solution(i) = 1;
    }
    fevals = 0;
  };

  //! Destructor
  ~Broyden() { cout << "Function evaluations: " << fevals << "\n"; };

  const NOX::Example::Vector& getInitialGuess()
  {
    return initialGuess;
  };

  const NOX::Example::Vector& getSolution()
  {
    return solution;
  };

  bool computeF(NOX::Example::Vector& f, const NOX::Example::Vector &x)
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
  
  bool computeJacobian(NOX::Example::Matrix& J, NOX::Example::Vector & x)
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
  NOX::Example::Vector initialGuess;
  //! Correct answer
  NOX::Example::Vector solution;

};

//! Main subroutine of Broyden.C
int main()
{
  // Set up the problem interface
  Broyden broyden(100,0.99);
  
  // Create a group which uses that problem interface. The group will
  // be initialized to contain the default initial guess for the
  // specified problem.
  NOX::Example::Group grp(broyden);

  // Set up the status tests
  //NOX::StatusTest::NormF statusTestA(grp, 1.0e-12);
  NOX::StatusTest::NormF statusTestA(1.0e-12, 
				     NOX::StatusTest::NormF::Unscaled);
  NOX::StatusTest::MaxIters statusTestB(150);
  NOX::StatusTest::Combo statusTestsCombo(NOX::StatusTest::Combo::OR, statusTestA, statusTestB);

  // Set up the solver parameters
  NOX::Parameter::List solverParameters;

  // Set the level of output (this is the default)
  solverParameters.setParameter("Output Information", 
		     NOX::Utils::Warning + 
		     NOX::Utils::OuterIteration + 
		     NOX::Utils::InnerIteration + 
		     NOX::Utils::Parameters);

  // Choose the solver itself
  solverParameters.setParameter("Nonlinear Solver", "Tensor Based");
  
  // Set up the line search parameters
  NOX::Parameter::List& lineSearchParameters = 
    solverParameters.sublist("Line Search");
  lineSearchParameters.setParameter("Method","Curvilinear");
  lineSearchParameters.setParameter("Max Iters",20);

  // Set up the direction parameters
  NOX::Parameter::List& directionParameters = 
    solverParameters.sublist("Direction");
  directionParameters.setParameter("Method","Tensor");

  // Set up the local solver parameters
  NOX::Parameter::List& localSolverParameters = 
    directionParameters.sublist("Linear Solver");
  localSolverParameters.setParameter("Tolerance",1e-4);
  localSolverParameters.setParameter("Reorthogonalize","As Needed");

  // Create the solver
  NOX::Solver::TensorBased solver(grp, statusTestsCombo, solverParameters);

  // Print the starting point
  cout << "\n" << "-- Starting Point --" << "\n";
  cout << "|| F(x0) || = " << NOX::Utils::sci(grp.getNormF()) << endl;
  // grp.print();

  // Solve the nonlinear system
  NOX::StatusTest::StatusType status = solver.solve();

  // Print the parameters
  cout << "\n" << "-- Parameter List Used in Solver --" << "\n";
  solver.getParameterList().print(cout);

  // Get the answer
  grp = solver.getSolutionGroup();

  // Print the answer
  cout << "\n" << "-- Final Solution From Solver --" << "\n";
  cout << "|| F(x*) || = " << NOX::Utils::sci(grp.getNormF()) << endl;
  // grp.print();

}
#endif
