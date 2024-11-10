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
#include "NOX_LAPACK_Group.H"
*/

#include "NOX_Common.H"
#include "NOX.H"  // NOX headers
#include "NOX_LAPACK.H" // NOX LAPACK Interface headers
#include "NOX_TestUtils.H" // NOX test parameter input headers
#include "Teuchos_StandardCatchMacros.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#else
#endif


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

int main(int argc, char *argv[]) {

  bool success = true;
  bool verbose = true;
  try {
    // Set up the printing utilities
    Teuchos::RCP<Teuchos::ParameterList> noxParamsPtr =
      Teuchos::rcp(new Teuchos::ParameterList);
    Teuchos::ParameterList& noxParams = *noxParamsPtr.get();
    Teuchos::ParameterList& printParams = noxParams.sublist("Printing");
    printParams.set("Output Precision", 5);

    std::string paramFilename;
    bool   usingParamInputFile = false;

    if (argc > 1) {
      if (argv[1][0]=='-' && argv[1][1]=='v')
        printParams.set("Output Information",
             NOX::Utils::OuterIteration +
             NOX::Utils::OuterIterationStatusTest +
             NOX::Utils::InnerIteration +
             NOX::Utils::Parameters +
             NOX::Utils::Details +
             NOX::Utils::Warning +
             NOX::Utils::TestDetails);
      else if (argv[1][0]=='-' && argv[1][1]=='p') {
        if (argc < 3) {
          std::cout << "Error: An input parameter file was expected but not found. \n" << std::endl;
          printParams.set("Output Information", NOX::Utils::Error);
          NOX::Utils printing(printParams);
          return EXIT_FAILURE;
          }
        paramFilename = argv[2];
        std::cout << "Reading parameter information from file \"" << paramFilename << "\""<< std::endl;
        usingParamInputFile = true;
        }
      else
        printParams.set("Output Information", NOX::Utils::Error);
    }
    NOX::Utils printing(printParams);

    // Identify the test
    if (printing.isPrintType(NOX::Utils::TestDetails)) {
      std::cout << "Starting lapack/NOX_NewTest/NOX_NewTest.exe" << std::endl;
    }

    // Final return value (0 = succefull, non-zero = failure)
    //int status = 0;

    // *** Insert your testing here! ***

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

    // Read parameters from file paramFilename - command line arg#1
    if (usingParamInputFile && !NOX::parseTextInputFile(paramFilename, noxParams))
       std::cout << "Using unchanged parameters " << std::endl;

    // Create the solver
    Teuchos::RCP<NOX::Solver::Generic> solver =
      NOX::Solver::buildSolver(grp, statusTestsCombo, noxParamsPtr);

    // Solve the nonlinesar system
    NOX::StatusTest::StatusType status = solver->solve();

    // Print the answer
    std::cout << "\n" << "-- Parameter List From Solver --" << "\n";
    solver->getList().print(std::cout);

    // Get the answer
    NOX::LAPACK::Group solnGrp =
      dynamic_cast<const NOX::LAPACK::Group&>(solver->getSolutionGroup());

    // Final return value (0 = succefull, non-zero = failure)
    if (status == NOX::StatusTest::Converged) {
      std::cout << "Test passed!" << std::endl;
      success = true;
    }
    else {
      std::cout << "Test failed!" << std::endl;
      success = false;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}

/*
  end of file main.cc
*/
