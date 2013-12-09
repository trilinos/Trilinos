// @HEADER
// ************************************************************************
// @HEADER


/*! \file  test_12.cpp
    \brief Check derivative checks.
*/

#include "ROL_StdVector.hpp"
#include "ROL_TestObjectives.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  int errorFlag  = 0;

  // Specify interval on which to generate uniform random numbers.
  RealT left = -1.0, right = 1.0;

  // *** Test body.

  try {

    int dim = 512;
    Teuchos::RCP<std::vector<RealT> > x_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    Teuchos::RCP<std::vector<RealT> > y_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    ROL::StdVector<RealT> x(x_rcp);
    ROL::StdVector<RealT> y(y_rcp);

    // set x,y
    for (int i=0; i<dim; i++) {
      (*x_rcp)[i] = 10.0* (1.0 + (RealT)rand() / (RealT)RAND_MAX);
      (*y_rcp)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
    }

    //ROL::Objective_Rosenbrock<RealT> obj;
    ROL::Objective_PoissonInversion<RealT> obj(dim,1.e-6);
    //ROL::Objective_SumOfSquares<RealT> obj;
    //ROL::Objective_LeastSquares<RealT> obj;

    std::vector<std::vector<RealT> > gCheck = obj.checkGradient(x, y);

    for (unsigned i=0; i<gCheck.size(); i++) {
      if (i==0) {
       	*outStream << std::right
                   << std::setw(20) << "Step size"
                   << std::setw(20) << "grad'*dir"
                   << std::setw(20) << "FD approx"
                   << std::setw(20) << "abs error"
       	       	   << "\n";
      }
      *outStream << std::scientific << std::setprecision(8) << std::right
                 << std::setw(20) << gCheck[i][0]
                 << std::setw(20) << gCheck[i][1]
                 << std::setw(20) << gCheck[i][2]
                 << std::setw(20) << gCheck[i][3]
                 << "\n";
    }

    *outStream << "\n";
    std::vector<std::vector<RealT> > hvCheck = obj.checkHessVec(x, y);

    for (unsigned i=0; i<hvCheck.size(); i++) {
      if (i==0) {
        *outStream << std::right
                   << std::setw(20) << "Step size"
                   << std::setw(20) << "norm(Hess*vec)"
                   << std::setw(20) << "norm(FD approx)"
                   << std::setw(20) << "norm(abs error)"
                   << "\n";
      }
      *outStream << std::scientific << std::setprecision(8) << std::right
                 << std::setw(20) << hvCheck[i][0]
                 << std::setw(20) << hvCheck[i][1]
                 << std::setw(20) << hvCheck[i][2]
                 << std::setw(20) << hvCheck[i][3]
                 << "\n";
    }

  
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  return 0;

}

