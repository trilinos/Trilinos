// $Id$ 
// $Source$ 

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
  \file Linear.C
  \brief A simple linear example based on NOX::LAPACK.
  
   This is an example of using %NOX with the NOX::LAPACK::Group and
   NOX::LAPACK::Vector classes. These are very basic classes intended
   only to illustrate and test NOX. They are based on a combination of
   C++ STL and LAPACK.

   This example is a simple linear problem.
   
   \f[
   F(x) = Ax -b 
   \f]

   If the program is called with the first argument being a 0 (zero),
   then the following problem is used:

   \f[
   A = 
   \left[
   \begin{array}{cccc}
   7 & 1 & 2 & 1\\
   1 & 8 & 3 & 1\\
   1 & 3 & 9 & 1\\
   1 & 1 & 1 & 10
   \end{array}
   \right]

   b = \left[
   \begin{array}{c}
   1 \\
   1 \\
   1 \\
   1
   \end{array}
   \right]

   \f]

   The initial guess is given by
   \f[
   x^{(0)} = \left[
   \begin{array}{c}
   4\\
   4\\
   4\\
   4
   \end{array}
   \right]
   \f]

   The solution is
   \f[
   x^{(*)} \approx \left[
   \begin{array}{c}
   0.1022\\
   0.0783\\
   0.0653\\
   0.0754
   \end{array}
   \right]
   \f]

   Otherwise, if the program is called with an n > 0, a matrix of size
   n is constructed and each element is filledwith a random elements
   in the interval [-1,1).

*/


#include "NOX.H"
#include "NOX_LAPACK_Vector.H"
#include "NOX_LAPACK_Group.H"
#include "NOX_LAPACK_Matrix.H"
#include "NOX_LAPACK_Wrappers.H"

//! An interface to the linear example described in Linear.C
class Linear : public NOX::LAPACK::Interface {

public:
 
  //! Constructor
  Linear(int n, bool useDefaultProblem) : 
    initialGuess(n),
    solution(n),
    A(n, n),
    b(n)
  {
    if (useDefaultProblem)
    {
      if (n != 4)
	throw "Linear::Linear - Problem size mismatch";

      initialGuess.init(4);
      
      b.init(1);
      
      A(0,0) = 7; 
      A(0,1) = 1;
      A(0,2) = 2;
      A(0,3) = 1;
      A(1,0) = 1; 
      A(1,1) = 8;
      A(1,2) = 3;
      A(1,3) = 1;
      A(2,0) = 1; 
      A(2,1) = 3;
      A(2,2) = 9;
      A(2,3) = 1;
      A(3,0) = 1; 
      A(3,1) = 1;
      A(3,2) = 1;
      A(3,3) = 10;

    }
    else
    {
      initialGuess.init(1);
      b.init(1);

      for (int i = 0; i < n; i ++)
	for (int j = 0; j < n; j ++)
	  A(i,j) = drand48();
    }

    int info;
    NOX::LAPACK::Matrix Acopy(A);
    vector<int> ipiv(n,0);
    solution = b;
    DGESV_F77(&n, &NOX::LAPACK::i_one, &Acopy(0,0), &n, &ipiv[0], &solution(0), &n, &info);
  };
  
  //! Destructor
  ~Linear() {};

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
    int n = f.length();
    for (int i = 0; i < n; i ++) 
    {
      f(i) = 0    ;
      for (int j = 0; j < n; j ++) {
	f(i) += A(i,j) * x(j);
      }
      f(i) -= b(i);
    }
    return true;
  };
  
  bool computeJacobian(NOX::LAPACK::Matrix& J, const NOX::LAPACK::Vector& x)
  {
    J = A;
    return true;
  };

private:

  //! Initial guess
  NOX::LAPACK::Vector initialGuess;
  //! Correct solution
  NOX::LAPACK::Vector solution;
  //! Matrix for linear problem
  /*! See Linear.C */
  NOX::LAPACK::Matrix A;
  //! RHS for linear problem
  /*! See Linear.C */
  NOX::LAPACK::Vector b;

};

//! Main subroutine in example Linear.C
int main(int argc, char* argv[])
{
  if (argc < 2)
  {
    cout << "Usage: " << argv[0] << " <linear problem size> [memory size]" << endl;
    exit(1);
  }

  // Set problem size and type
  int n = atoi(argv[1]);
  bool useDefaultProblem = false;
  if (n <= 0)
  {
    n = 4;
    useDefaultProblem = true;
  }


  // Set up the problem interface
  Linear linear(n, useDefaultProblem);
  
  // Create a group which uses that problem interface. The group will
  // be initialized to contain the default initial guess for the
  // specified problem.
  NOX::LAPACK::Group grp(linear);

  // Set up the status tests
  NOX::StatusTest::NormF normf(grp, 1.0e-8);
  NOX::StatusTest::MaxIters maxiters( 5 * n );
  NOX::StatusTest::Combo statusTestsCombo(NOX::StatusTest::Combo::OR, normf, maxiters);

  // Create the list of solver parameters
  NOX::Parameter::List solverParameters;

  // Set the level of output (this is the default)
  solverParameters.sublist("Printing").setParameter("Output Information", 
						    NOX::Utils::Warning + 
						    NOX::Utils::OuterIteration + 
						    NOX::Utils::InnerIteration + 
						    NOX::Utils::Parameters);

  // Set the solver (this is the default)
  solverParameters.setParameter("Nonlinear Solver", "Line Search Based");

  // Create direction sublist
  NOX::Parameter::List& directionParameters = solverParameters.sublist("Direction");
  
  // Set direction
  directionParameters.setParameter("Method", "Quasi-Newton");

  // Set memory size
  int m = 5;
  if (argc >= 3) 
    m = atoi(argv[2]);
  if (m <= 0)
    m = 5;
  directionParameters.sublist("Quasi-Newton").setParameter("Memory", m);

  // Create the line search parameters sublist
  NOX::Parameter::List& lineSearchParameters = solverParameters.sublist("Line Search");

  // Set the line search method
  lineSearchParameters.setParameter("Method", "Polynomial");

  // Create the line search parameters sublist
  NOX::Parameter::List& polynomialParameters = lineSearchParameters.sublist("Polynomial");

  polynomialParameters.setParameter("Force Interpolation", true);
  polynomialParameters.setParameter("Convergence Criteria", "None");
  polynomialParameters.setParameter("Default Step", 1.0);
  polynomialParameters.setParameter("Min Bounds Factor", 0.0);

  // Create the solver
  NOX::Solver::Manager solver(grp, statusTestsCombo, solverParameters);

  // Solve the nonlinesar system
  NOX::StatusTest::StatusType status = solver.solve();

  // Print the answer
  cout << "\n" << "-- Parameter List From Solver --" << "\n";
  solver.getParameterList().print(cout);

  // Get the answer
  grp = solver.getSolutionGroup();
  

  // Print the answer
  cout << "\n" << "-- Final Solution From Solver --" << "\n";
  grp.print();

  // Print the expected answer
  grp.setX(linear.getSolution());
  grp.computeF();
  cout << "\n" << "-- Expected Solution --" << "\n";
  grp.print();
}

#endif
