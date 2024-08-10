// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*!
   \file Rosenbrock.C
   \brief A simple 2D example using %Rosenbrock's function based on NOX::Example

   This is an example of using %NOX with the NOX::LAPACK::Group and
   NOX::LAPACK::Vector classes. These are very basic classes intended
   only to illustrate and test NOX. They are based on a combination of
   C++ STL and LAPACK.

   This example is the "%Rosenbrock function" from Jorge J. More',
   Burton S. Garbow, and Kenneth E. Hillstrom, Testing Unconstrained
   Optimization Software, ACM TOMS, Vol. 7, No. 1, March 1981,
   pp. 14-41.

   It comes originally from H. H. %Rosenbrock, An Automatic Method for
   Finding the Greatest or Least Value of a Function, J. Comput.
   3(1960):175-184.

   The function is defined as
   \f[
   F(x) = \left[
   \begin{array}{c}
   10 (x[2] - x[1]^2) \\
   1 - x[1]
   \end{array}
   \right]
   \f]

   The initial guess is given by
   \f[
   x_0 = \left[
   \begin{array}{c}
   -1.2\\
   1
   \end{array}
   \right]
   \f]

   The solution is
   \f[
   x_* = \left[
   \begin{array}{c}
   1\\
   1
   \end{array}
   \right]
   \f]

*/


#include "NOX.H"
#include "NOX_LAPACK_Group.H"

//! An interface to the example described in Rosenbrock.C
class Rosenbrock : public NOX::LAPACK::Interface {

public:

  //! Constructor
  Rosenbrock() :
    initialGuess(2),
    solution(2)
  {
    initialGuess(0) = -1.2;
    initialGuess(1) = 1;
    solution(0) = 1;
    solution(1) = 1;
  };

  //! Destructor
  ~Rosenbrock() {};

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
    f(0) = 10 * (x(1) - x(0) * x(0));
    f(1) = 1 - x(0);
    return true;
  };

  bool computeJacobian(NOX::LAPACK::Matrix<double>& J,
               const NOX::LAPACK::Vector & x)
  {
    J(0,0) = -20 * x(0);
    J(0,1) = 10;
    J(1,0) = -1;
    J(1,1) = 0;
    return true;
  };

private:

  //! Initial guess
  NOX::LAPACK::Vector initialGuess;
  //! Correct solution
  NOX::LAPACK::Vector solution;

};

//! Main subroutine of Rosenbrock.C
int main()
{
  // Set up the problem interface
  Rosenbrock rosenbrock;

  // Create a group which uses that problem interface. The group will
  // be initialized to contain the default initial guess for the
  // specified problem.
  Teuchos::RCP<NOX::LAPACK::Group> grp =
    Teuchos::rcp(new NOX::LAPACK::Group(rosenbrock));

  // Set up the status tests
  Teuchos::RCP<NOX::StatusTest::NormF> statusTestA =
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-4));
  Teuchos::RCP<NOX::StatusTest::MaxIters> statusTestB =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
  Teuchos::RCP<NOX::StatusTest::Combo> statusTestsCombo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,
                        statusTestA, statusTestB));

  // Create the list of solver parameters
  Teuchos::RCP<Teuchos::ParameterList> solverParametersPtr =
    Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& solverParameters = *solverParametersPtr;

  // Set the level of output (this is the default)
  solverParameters.sublist("Printing").set("Output Information",
                            NOX::Utils::Warning +
                            NOX::Utils::OuterIteration +
                            NOX::Utils::InnerIteration +
                            NOX::Utils::Parameters);

  // Set the solver (this is the default)
  solverParameters.set("Nonlinear Solver", "Line Search Based");

  // Create the line search parameters sublist
  Teuchos::ParameterList& lineSearchParameters = solverParameters.sublist("Line Search");

  // Set the line search method
  lineSearchParameters.set("Method","More'-Thuente");

  // Create the solver
  NOX::Solver::Manager solver(grp, statusTestsCombo, solverParametersPtr);

  // Solve the nonlinesar system
  NOX::StatusTest::StatusType status = solver.solve();

  // Warn user if solve failed
  if (status == NOX::StatusTest::Converged)
    std::cout << "Example Passed!" << std::endl;
  else
    std::cout << "Error: Solve failed to converge!" << std::endl;

  // Print the parameter list
  std::cout << "\n" << "-- Parameter List From Solver --" << "\n";
  solver.getList().print(cout);

  // Get the answer
  NOX::LAPACK::Group solnGrp =
    dynamic_cast<const NOX::LAPACK::Group&>(solver.getSolutionGroup());

  // Print the answer
  std::cout << "\n" << "-- Final Solution From Solver --" << "\n";
  solnGrp.print();

  // Print the expected answer
  solnGrp.setX(rosenbrock.getSolution());
  solnGrp.computeF();
  std::cout << "\n" << "-- Expected Solution --" << "\n";
  solnGrp.print();
}
