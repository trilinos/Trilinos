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

   The standard starting point is x(i) = -1, but setting x(i) = 0 tests
   the selected global strategy.

   \author Brett Bader, CU Boulder, 2002
*/


#include "NOX.H"
#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "NOX_LAPACK_Group.H"

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
  
  bool computeJacobian(NOX::LAPACK::Matrix& J, const NOX::LAPACK::Vector & x)
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

//! Main subroutine of Broyden.C
int main()
{
  // Set up the problem interface
  Broyden broyden(100,0.99);
  
  // Create a group which uses that problem interface. The group will
  // be initialized to contain the default initial guess for the
  // specified problem.
  NOX::LAPACK::Group grp(broyden);

  // Create the top level parameter list
  NOX::Parameter::List solverParameters;

  // Set the nonlinear solver method
  //solverParameters.setParameter("Nonlinear Solver", "Tensor-Krylov Based");
  solverParameters.setParameter("Nonlinear Solver", "Tensor Based");
  
  // Sublist for printing parameters
  NOX::Parameter::List& printParams = solverParameters.sublist("Printing");
  //printParams.setParameter("MyPID", 0); 
  printParams.setParameter("Output Precision", 3);
  printParams.setParameter("Output Processor", 0);
  printParams.setParameter("Output Information", 
			NOX::Utils::OuterIteration + 
			NOX::Utils::OuterIterationStatusTest + 
			NOX::Utils::InnerIteration +
			NOX::Utils::Parameters + 
			NOX::Utils::Details + 
			NOX::Utils::Warning);
  NOX::Utils utils(printParams);

  // Sublist for direction parameters
  NOX::Parameter::List& directionParameters = 
    solverParameters.sublist("Direction");
  directionParameters.setParameter("Method","Tensor");
  directionParameters.setParameter("Compute Step","Tensor");

  // Sublist for local solver parameters
  NOX::Parameter::List& localSolverParameters = 
    directionParameters.sublist("Tensor").sublist("Linear Solver");
  //localSolverParameters.setParameter("Tolerance", 1e-4);
  //localSolverParameters.setParameter("Reorthogonalize","Always");
  //localSolverParameters.setParameter("Output Frequency",1);
  //localSolverParameters.setParameter("Max Restarts", 2);
  //localSolverParameters.setParameter("Size of Krylov Subspace", 15);
  //localSolverParameters.setParameter("Preconditioning","Tridiagonal");
  //localSolverParameters.setParameter("Preconditioning Side","None");
  //localSolverParameters.setParameter("Use Shortcut Method",false);

  // Sublist for line search parameters
  NOX::Parameter::List& lineSearchParameters = 
    solverParameters.sublist("Line Search");
  lineSearchParameters.setParameter("Method","Tensor");
  lineSearchParameters.sublist("Tensor").setParameter("Submethod","Curvilinear");
  lineSearchParameters.sublist("Tensor").setParameter("Lambda Selection","Halving");
  lineSearchParameters.sublist("Tensor").setParameter("Max Iters",20);


  // Create the convergence tests
  //NOX::StatusTest::NormF statusTestA(grp, 1.0e-12);
  NOX::StatusTest::NormF statusTestA(1.0e-12, NOX::StatusTest::NormF::Unscaled);
  NOX::StatusTest::MaxIters statusTestB(50);
  NOX::StatusTest::Combo statusTestsCombo(NOX::StatusTest::Combo::OR,
					  statusTestA, statusTestB);

  // Create the solver
  NOX::Solver::Manager solver(grp, statusTestsCombo, solverParameters);

  // Print the starting point
  cout << "\n" << "-- Starting Point --" << "\n";
  cout << "|| F(x0) || = " << utils.sciformat(grp.getNormF()) << endl;
  // grp.print();

  // Solve the nonlinear system
  NOX::StatusTest::StatusType status = solver.solve();

  // Get the answer
  grp = solver.getSolutionGroup();

  // Output the parameter list
  if (utils.isPrintProcessAndType(NOX::Utils::Parameters)) {
    cout << "\n" << "-- Parameter List Used in Solver --" << endl;
    solver.getParameterList().print(cout);
    cout << endl;
  }

  // Print the answer
  if (utils.isPrintProcessAndType(NOX::Utils::Parameters)) {
    cout << "\n" << "-- Final Solution From Solver --" << "\n";
    cout << "|| F(x*) || = " << utils.sciformat(grp.getNormF()) << endl;
    // grp.print();
  }
  
}
