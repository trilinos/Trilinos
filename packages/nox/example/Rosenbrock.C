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

/* 
   This is an example of using NOX with the NOX::Example group and
   vector classes. These are very basic classes intended only to
   illustrate and test NOX. They are based on LAPACK.

   This example is the "Rosenbrock function" from Jorge J. More',
   Burton S. Garbow, and Kenneth E. Hillstrom, Testing Unconstrained
   Optimization Software, ACM TOMS, Vol. 7, No. 1, March 1981,
   pp. 14-41.

   It comes originally from H. H. Rosenbrock, An Automatic Method for
   Finding the Greatest or Least Value of a Function, J. Comput.
   3(1960):175-184.
*/


#include "NOX.H"
#include "NOX_Example_Group.H"

class Rosenbrock : public NOX::Example::Interface {

public:
 
  Rosenbrock() : 
    initialGuess(2),
    solution(2)
  {
    initialGuess(0) = -1.2;
    initialGuess(1) = -1;
    solution(0) = 1;
    solution(1) = 1;
  };

  ~Rosenbrock() {};

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
    f(0) = 10 * (x(1) - x(0) * x(0));
    f(1) = 1 - x(0);
    return true;
  };
  
  bool computeJacobian(NOX::Example::Matrix& J, NOX::Example::Vector & x)
  {
    J(0,0) = -20 * x(0);
    J(0,1) = 10;
    J(1,0) = -1;
    J(1,1) = 0;
    return true;
  };

private:

  NOX::Example::Vector initialGuess;
  NOX::Example::Vector solution;

};

int main()
{
  // Set up the problem interface
  Rosenbrock rosenbrock;;
  
  // Create a group which uses that problem interface. The group will
  // be initialized to contain the default initial guess for the
  // specified problem.
  NOX::Example::Group grp(rosenbrock);

  // Set up the status tests
  NOX::StatusTest::NormF statusTestA(grp, 1.0e-4);
  NOX::StatusTest::MaxIters statusTestB(10);
  NOX::StatusTest::Combo statusTestsCombo(NOX::StatusTest::Combo::OR, statusTestA, statusTestB);

  // Set up the solver parameters
  NOX::Parameter::List solverParameters;

  // Set the level of output (this is the default)
  solverParameters.setParameter("Output Information", 
		     NOX::Utils::Warning + 
		     NOX::Utils::OuterIteration + 
		     NOX::Utils::InnerIteration + 
		     NOX::Utils::Parameters);

  // Choose the solver itself (this is the default)
  solverParameters.setParameter("Solver", "Line Search");

  // Set up the line search parameters
  NOX::Parameter::List& lineSearchParameters = solverParameters.sublist("Line Search");
  lineSearchParameters.setParameter("Method","More'-Thuente");

  // Create the solver
  NOX::Solver::Manager solver(grp, statusTestsCombo, solverParameters);

  // Solve the nonlinesar system
  NOX::StatusTest::StatusType status = solver.solve();

  // Get the answer
  grp = solver.getSolutionGroup();
  
  // Print the answer
  cout << "\n" << "-- Final Solution From Solver --" << "\n";
  grp.print();

  // Print the expected answer
  grp.setX(rosenbrock.getSolution());
  grp.computeF();
  cout << "\n" << "-- Expected Solution --" << "\n";
  grp.print();
}
