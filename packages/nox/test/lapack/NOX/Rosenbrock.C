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
  
  bool computeJacobian(NOX::LAPACK::Matrix& J, const NOX::LAPACK::Vector & x)
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

  // Set up the printing utilities
  NOX::Parameter::List noxParams;
  NOX::Parameter::List& printParams = noxParams.sublist("Printing");
  printParams.setParameter("Output Precision", 5);
       
  string paramFilename;     
  bool   usingParamInputFile = false;

  if (argc > 1) { 
    if (argv[1][0]=='-' && argv[1][1]=='v')
       printParams.setParameter("Output Information", 
			NOX::Utils::OuterIteration + 
			NOX::Utils::OuterIterationStatusTest + 
			NOX::Utils::InnerIteration +
			NOX::Utils::Parameters + 
			NOX::Utils::Details + 
			NOX::Utils::Warning +
			NOX::Utils::TestDetails);
    else if (argv[1][0]=='-' && argv[1][1]=='p')
      {

       if (argc < 3)
         {
            cout << "Error: An input parameter file was expected but not found. \n" << endl;
            printParams.setParameter("Output Information", NOX::Utils::Error);
	    NOX::Utils printing(printParams);
            return 0;
         }

       paramFilename = argv[2];
       cout << "Reading parameter information from file \"" << paramFilename << "\""<< endl;
       usingParamInputFile = true;

      }
    else
       printParams.setParameter("Output Information", NOX::Utils::Error);
  }
  NOX::Utils printing(printParams);

  // This test is only for SERIAL
  // Exit if trying to run in parallel
#ifdef HAVE_MPI
  return 0;
#endif
  
  // Identify the test
  if (printing.isPrintProcessAndType(NOX::Utils::TestDetails)) {
    cout << "Starting lapack/NOX_NewTest/NOX_NewTest.exe" << endl;
  }

  // Final return value (0 = succefull, non-zero = failure)
  //int status = 0;

  // *** Insert your testing here! ***

  // Set up the problem interface
  Rosenbrock rosenbrock;
  
  // Create a group which uses that problem interface. The group will
  // be initialized to contain the default initial guess for the
  // specified problem.
  NOX::LAPACK::Group grp(rosenbrock);

  // Set up the status tests
  NOX::StatusTest::NormF statusTestA(grp, 1.0e-4);
  NOX::StatusTest::MaxIters statusTestB(20);
  NOX::StatusTest::Combo statusTestsCombo(NOX::StatusTest::Combo::OR, statusTestA, statusTestB);


  // Create the list of solver parameters
  NOX::Parameter::List solverParameters;

  // Read parameters from file paramFilename - command line arg#1
  if (usingParamInputFile && !NOX::parseTextInputFile(paramFilename, solverParameters))
     cout << "Using unchanged parameters " << endl;
  
  // Create the solver
  NOX::Solver::Manager solver(grp, statusTestsCombo, solverParameters);

  // Solve the nonlinesar system
  NOX::StatusTest::StatusType status = solver.solve();

  // Print the answer
  cout << "\n" << "-- Parameter List From Solver --" << "\n";
  solver.getParameterList().print(cout);

  // Get the answer
  grp = solver.getSolutionGroup();
  


  if (printing.isPrintProcessAndType(NOX::Utils::TestDetails)) {
    if (status == NOX::StatusTest::Converged) 
      cout << "Test was successful!" << endl;
    else 
      cout << "Test Failed!" << endl;
  }

  // Final return value (0 = succefull, non-zero = failure)
  //return status;
  if (status ==  NOX::StatusTest::Converged) 
     return 0;
  else
     return 1;
}

/*
  end of file main.cc
*/
