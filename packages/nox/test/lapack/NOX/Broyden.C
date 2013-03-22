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
#include "NOX_TestUtils.H"

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
  };

  //! Destructor
  ~Broyden() { std::cout << "Function evaluations: " << fevals << "\n"; };

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

//! Main subroutine of Broyden.C
int main(int argc, char *argv[])
{
  std::cout << "Started" << std::endl;

  // Set up the problem interface
  Broyden broyden(100,0.99);
  
  // Create a group which uses that problem interface. The group will
  // be initialized to contain the default initial guess for the
  // specified problem.
  Teuchos::RCP<NOX::LAPACK::Group> grp = 
    Teuchos::rcp(new NOX::LAPACK::Group(broyden));

  // Create the top level parameter list
  Teuchos::RCP<Teuchos::ParameterList> solverParametersPtr =
    Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& solverParameters = *solverParametersPtr;

  // Set the nonlinear solver method
  //solverParameters.set("Nonlinear Solver", "Tensor-Krylov Based");
  solverParameters.set("Nonlinear Solver", "Tensor Based");
  
  // Sublist for printing parameters
  Teuchos::ParameterList& printParams = solverParameters.sublist("Printing");
  //printParams.set("MyPID", 0); 
  printParams.set("Output Precision", 3);
  printParams.set("Output Processor", 0);
  printParams.set("Output Information", 
			NOX::Utils::OuterIteration + 
			NOX::Utils::OuterIterationStatusTest + 
			NOX::Utils::InnerIteration +
			NOX::Utils::Parameters + 
			NOX::Utils::Details + 
			NOX::Utils::Warning);
  NOX::Utils utils(printParams);

  // Sublist for direction parameters
  Teuchos::ParameterList& directionParameters = 
    solverParameters.sublist("Direction");
  directionParameters.set("Method","Tensor");

  // Sublist for local solver parameters
  //Teuchos::ParameterList& localSolverParameters = 
  //directionParameters.sublist("Tensor").sublist("Linear Solver");
  //localSolverParameters.set("Tolerance", 1e-4);
  //localSolverParameters.set("Reorthogonalize","Always");
  //localSolverParameters.set("Output Frequency",1);
  //localSolverParameters.set("Max Restarts", 2);
  //localSolverParameters.set("Size of Krylov Subspace", 15);
  //localSolverParameters.set("Preconditioning","Tridiagonal");
  //localSolverParameters.set("Preconditioning Side","None");
  //localSolverParameters.set("Use Shortcut Method",false);

  // Sublist for line search parameters
  Teuchos::ParameterList& globalStrategyParameters = 
    solverParameters.sublist("Line Search");
  globalStrategyParameters.set("Method","Curvilinear");
  
  // Sublist for line search parameters
  Teuchos::ParameterList& lineSearchParameters =
    globalStrategyParameters.sublist(globalStrategyParameters.
				     get("Method","Curvilinear"));

  lineSearchParameters.set("Lambda Selection","Halving");
  lineSearchParameters.set("Max Iters",20);

  // 
  // Update parameters from an input file if the input file was provided in command
  // line. Usage -p paramFilename 
  //
  std::string paramFilename;     
  bool   usingParamInputFile = false;
  if (argc > 1)
    {
      if (argv[1][0]=='-' && argv[1][1]=='p')
        {

          if (argc < 3)
            {
              std::cout << "Error: A parameter input file was expected but not found. \n" << std::endl;
              printParams.set("Output Information", NOX::Utils::Error);
	      NOX::Utils printing(printParams);
              return 0;
            }

          paramFilename = argv[2];
          std::cout << "Reading parameter information from file \"" << paramFilename << "\""<< std::endl;
          usingParamInputFile = true;

        }
    }

  // Read parameters from file paramFilename - command line arg#1
  if (usingParamInputFile && !NOX::parseTextInputFile(paramFilename, solverParameters))
     std::cout << "Using unchanged parameters " << std::endl;

  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::NormF> statusTestA = 
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-12, 
					    NOX::StatusTest::NormF::Unscaled));
  Teuchos::RCP<NOX::StatusTest::MaxIters> statusTestB = 
    Teuchos::rcp(new NOX::StatusTest::MaxIters(50));
  Teuchos::RCP<NOX::StatusTest::Combo> statusTestsCombo = 
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,
					    statusTestA, statusTestB));
  
  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver = 
    NOX::Solver::buildSolver(grp, statusTestsCombo, solverParametersPtr);

  // Print the starting point
  grp->computeF();
  std::cout << "\n" << "-- Starting Point --" << "\n";
  std::cout << "|| F(x0) || = " << utils.sciformat(grp->getNormF()) << std::endl;
  // grp.print();

  // Solve the nonlinear system
  NOX::StatusTest::StatusType status = solver->solve();

  // Get the answer
  NOX::LAPACK::Group solnGrp = 
    dynamic_cast<const NOX::LAPACK::Group&>(solver->getSolutionGroup());

  // Output the parameter list
  if (utils.isPrintType(NOX::Utils::Parameters)) {
    std::cout << "\n" << "-- Parameter List Used in Solver --" << std::endl;
    solver->getList().print(std::cout);
    std::cout << std::endl;
  }

  // Print the answer
  if (utils.isPrintType(NOX::Utils::Parameters)) {
    std::cout << "\n" << "-- Final Solution From Solver --" << "\n";
    std::cout << "|| F(x*) || = " << utils.sciformat(solnGrp.getNormF()) << std::endl;
    // solnGrp.print();
  }
  
  // Print final status
  int returnValue = 1;
  if (status == NOX::StatusTest::Converged) {
    returnValue = 0;
    std::cout << "Test passed!" << std::endl;    
  }
  else 
    std::cout << "Test failed!" << std::endl;

  return returnValue;
}
