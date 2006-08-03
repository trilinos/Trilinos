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
  \file Linear.C

  \brief A simple linear example based on NOX::LAPACK. Creates the
  executable named \c linear.
  
  <b>Usage</b>

  <tt>linear \<problem size\> [line search type] [direction type] [memory size] </tt>

  <ul>

  <li> <tt>problem size</tt> - Specify the size of the linear problem

  <li> <tt>line search type</tt> - Optional argument specifying the type
       of line search to be used. Choices are the following:

       <ul>
       <li> None - Always take a step of length 1. Uses NOX::LineSearch::FullStep.

       <li> Exact - Always perform an exact quadratic line search,
            even if the first step decreases the function value.  Uses
            NOX::LineSearch::Polynomial the following settings:

	    <ul>
	    <li>"Force Interpolation" = true
	    <li>"Sufficient Decrease Condition" = "None"
	    <li>"Default Step"= 1.0
	    <li>"Min Bounds Factor" = 0.0
	    </ul>

       <li> Poly - Uses NOX::LineSearch::Polynomial with the default
            settings.

       </ul>

  <li> <tt>direction type</tt> - Optional argument specifying the type
       of direction to be used. Choices are the following:

       <ul>
       <li>Newton - NOX::Direction::Newton [default]
       <li>Quasi-Newton - NOX::Direction::QuasiNewton
       <li>%Broyden - NOX::Direction::Broyden
       </ul>

  <li> <tt>memory size</tt> - Optional argument to specify the memory
       size for the Quasi-Newton and %Broyden methods. (Defaults to
       the default for the given direction.)

  </ul>

  <b>Description</b>

  This is an example of using %NOX with the NOX::LAPACK::Group and
  NOX::LAPACK::Vector classes. These are very basic classes intended
  only to illustrate and test %NOX. They are based on a combination of
  C++ STL and LAPACK.

  The function is linear:
  \f[
  F(x) = Ax -b.
  \f]
  The matrix \f$A\f$ is filled with random values in the interval [-1,1) using NOX::Random.
  The right hand side is specified by
   \f[
   b = \left[
   \begin{array}{c}
   1 \\
   1 \\
   1 \\
   1
   \end{array}
   \right],
   \f]
   and the initial guess is given by
   \f[
   x^{(0)} = \left[
   \begin{array}{c}
   0\\
   0\\
   0\\
   0
   \end{array}
   \right].
   \f]


   <b>Special Case</b>

  If the program is called with problem size set to zero,
  then the following 4 x 4 problem is used instead:

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

*/


#include "NOX.H"
#include "NOX_Random.H"
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
	throw "Linear::Linear (from NOX) - Problem size mismatch";

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
      if (n < 0)
	throw "Linear::Linear (from NOX) - Invalid problem size";

      initialGuess.init(1);
      b.init(1);
      NOX::Random random(1);

      for (int i = 0; i < n; i ++)
	for (int j = 0; j < n; j ++)
	  A(i,j) = random.number();
    }

    int info;
    NOX::LAPACK::Matrix<double> Acopy(A);
    vector<int> ipiv(n,0);
    solution = b;
    DGESV_F77(&n, &NOX::LAPACK::i_one, &Acopy(0,0), &n, &ipiv[0], &solution(0), &n, &info);
  };
  
  //! Destructor
  ~Linear() {};

  // Derived
  const NOX::LAPACK::Vector& getInitialGuess()
  {
    return initialGuess;
  };

  //! Return true solution vector
  const NOX::LAPACK::Vector& getSolution()
  {
    return solution;
  };

  // Derived
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
  
  // Derived
  bool computeJacobian(NOX::LAPACK::Matrix<double>& J, 
		       const NOX::LAPACK::Vector& x)
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
  NOX::LAPACK::Matrix<double> A;
  //! RHS for linear problem
  NOX::LAPACK::Vector b;

};

//! Main subroutine in example Linear.C
int main(int argc, char* argv[])
{
  if (argc < 2)
  {
    cout << "Usage: " << argv[0] << " <problem size> [line search type] [direction type] [memory size]" << endl;
    exit(1);
  }

  // ** Problem Set Up **

  // Set problem size and type
  int n = atoi(argv[1]);
  bool useDefaultProblem = false;

  if (n == 0)
  {
    n = 4;
    useDefaultProblem = true;
  }


  // Set up the problem interface
  Linear linear(n, useDefaultProblem);
  
  // Create a group which uses that problem interface. The group will
  // be initialized to contain the default initial guess for the
  // specified problem.
  Teuchos::RefCountPtr<NOX::LAPACK::Group> grp = 
    Teuchos::rcp(new NOX::LAPACK::Group(linear));
  
  // ** Status Tests **
  Teuchos::RefCountPtr<NOX::StatusTest::NormF> normf =
    Teuchos::rcp(new NOX::StatusTest::NormF(*grp, 1.0e-8));
  Teuchos::RefCountPtr<NOX::StatusTest::MaxIters> maxiters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters( 5 * n ));
  Teuchos::RefCountPtr<NOX::StatusTest::Combo> statusTestsCombo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR, 
					    normf, maxiters));

  // ** Paramter List **
  Teuchos::RefCountPtr<Teuchos::ParameterList> solverParamsPtr =
    Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& solverParams = *solverParamsPtr.get();

  // -- Output Level --
  solverParams.sublist("Printing").set("Output Information", 
						    NOX::Utils::Warning + 
						    NOX::Utils::OuterIteration + 
						    NOX::Utils::InnerIteration + 
						    NOX::Utils::Parameters);

  // -- Solver --
  solverParams.set("Nonlinear Solver", "Line Search Based");

  // -- Line Search --
  Teuchos::ParameterList& lineSearchParams = solverParams.sublist("Line Search");


  string lineSearchType = (argc >= 3) ? argv[2] : "Exact";

  if (lineSearchType == "None")
  {
    lineSearchParams.set("Method", "Full Step");
  }
  else if (lineSearchType == "Poly")
  {
    lineSearchParams.set("Method", "Polynomial");
  }
  else if (lineSearchType == "Backtrack")
  {
    lineSearchParams.set("Method", "Backtrack");
  }
  else // "Exact"
  {
    lineSearchParams.set("Method", "Polynomial");
    Teuchos::ParameterList& polyParams = lineSearchParams.sublist("Polynomial");
    polyParams.set("Force Interpolation", true);
    polyParams.set("Sufficient Decrease Condition", "None");
    polyParams.set("Default Step", 1.0);
    polyParams.set("Min Bounds Factor", 0.0);
  }

  // -- Direction --
  Teuchos::ParameterList& directionParams = solverParams.sublist("Direction");
  
  // Set direction type
  string directionMethod = (argc >= 4) ? argv[3] : "Newton";
  directionParams.set("Method", directionMethod);

  // Set memory size
  if (argc >= 5) 
    directionParams.sublist(directionMethod).set("Memory", atoi(argv[4]));

  // One last detail
  if ((directionMethod == "Broyden") && (lineSearchType == "Poly"))
    directionParams.sublist(directionMethod).set("Compute Jacobian", true);

  // ** Solve **
  
  // Create the solver
  NOX::Solver::Manager solver(grp, statusTestsCombo, solverParamsPtr);

  // Solve the nonlinesar system
  NOX::StatusTest::StatusType status = solver.solve();

  // Warn user if solve failed
  if (status != NOX::StatusTest::Converged)
    cout << "Error: Solve failed to converge!" << endl;

  // ** Output solution **

  // Print the final parameter list from the solver
  cout << "\n" << "-- Parameter List From Solver --" << "\n";
  solver.getList().print(cout);

  // Get the answer from the solver
  NOX::LAPACK::Group solnGrp = 
    dynamic_cast<const NOX::LAPACK::Group&>(solver.getSolutionGroup());
  
  // Print the answer from the solver
  cout << "\n" << "-- Final Solution From Solver --" << "\n";
  solnGrp.print();

  // Print the true answer
  solnGrp.setX(linear.getSolution());
  solnGrp.computeF();
  cout << "\n" << "-- Expected Solution --" << "\n";
  solnGrp.print();
}

#endif
