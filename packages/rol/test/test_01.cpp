// @HEADER
// ************************************************************************
// @HEADER


/*! \file  test_01.cpp
    \brief Test std::vector interface.
*/

#include "ROL_StdVector.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>

typedef float  RealT;
typedef double ElementT;

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

  int errorFlag  = 0;

  // *** Test body.

  try {

    int dim = 100;
    Teuchos::RCP<std::vector<ElementT> > x_rcp = Teuchos::rcp( new std::vector<ElementT> (dim, 0.0) );
    Teuchos::RCP<std::vector<ElementT> > y_rcp = Teuchos::rcp( new std::vector<ElementT> (dim, 0.0) );
    ROL::StdVector<RealT, ElementT> x(x_rcp);
    ROL::StdVector<RealT, ElementT> y(y_rcp);

    // set x,y
    for (int i=0; i<dim; i++) {
      (*x_rcp)[i] = i;
      (*y_rcp)[i] = 2.0;
    }

    // norm of x
    RealT xnorm = x.norm();
    *outStream << "\nNorm of ROL::StdVector x: " << xnorm << "\n";

    // norm of y
    RealT ynorm = y.norm();
    *outStream << "\nNorm of ROL::StdVector y: " << ynorm << "\n";

    // scale x
    x.scale(0.5);
    xnorm = x.norm();
    *outStream << "\nNorm of half of x: " << xnorm << "\n";

    // clone z from x, deep copy x into z, norm of z
    Teuchos::RCP<ROL::Vector<RealT> > z = x.clone();
    z->set(x);
    RealT znorm = z->norm();
    *outStream << "\nNorm of ROL::Vector z (clone of x): " << znorm << "\n";

    // compute norm of x - x - 0
    z->set(x);
    x.scale(-1.0);
    z->plus(x);
    y.zero();
    z->axpy(-1.0, y);
    znorm = z->norm();
    *outStream << "\nNorm of (x - x) - 0: " << znorm << "\n";

    // set x to dim/2-th basis vector
    z = x.basis(dim/2);
    znorm = z->norm();
    *outStream << "\nNorm of ROL::Vector z (basis vector): " << znorm << "\n";

  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;

}

