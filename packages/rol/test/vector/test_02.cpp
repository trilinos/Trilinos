// @HEADER
// ************************************************************************
// @HEADER


/*! \file  test_02.cpp
    \brief Test Epetra_MultiVector interface.
*/

#include "ROL_EpetraMultiVector.hpp"
#include "ROL_Types.hpp"
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>

typedef double RealT;
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

  double errtol = ROL::ROL_THRESHOLD;

  // *** Test body.

  try {

    int dim = 100;
    Epetra_SerialComm Comm;
    Epetra_Map Map(dim, 0, Comm);

    Teuchos::RCP<Epetra_MultiVector> x_rcp = Teuchos::rcp( new Epetra_MultiVector(Map, 1) );
    Teuchos::RCP<Epetra_MultiVector> y_rcp = Teuchos::rcp( new Epetra_MultiVector(Map, 1) );
    ROL::EpetraMultiVector<RealT> x(x_rcp);
    ROL::EpetraMultiVector<RealT> y(y_rcp);

    // set x,y
    for (int i=0; i<dim; i++) {
      ((*x_rcp)[0])[i] = i;
      ((*y_rcp)[0])[i] = 2.0;
    }

    // norm of x
    RealT xnorm = x.norm();
    *outStream << "\nNorm of ROL::EpetraMultiVector x: " << xnorm << "\n";

    // norm of y
    RealT ynorm = y.norm();
    *outStream << "\nNorm of ROL::EpetraMultiVector y: " << ynorm << "\n";

    // scale x
    x.scale(0.5);
    RealT xnorm2 = x.norm();
    *outStream << "\nNorm of half of x: " << xnorm2 << "\n";
    if ( std::abs(xnorm/xnorm2 - 2.0) > errtol ) {
      *outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    };

    // clone z from x, deep copy x into z, norm of z
    Teuchos::RCP<ROL::Vector<RealT> > z = x.clone();
    z->set(x);
    RealT znorm = z->norm();
    *outStream << "\nNorm of ROL::Vector z (clone of x): " << znorm << "\n";
    if ( std::abs(xnorm2 - znorm) > errtol ) {
      *outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    };

    // compute norm of x - x - 0
    z->set(x);
    x.scale(-1.0);
    z->plus(x);
    y.zero();
    z->axpy(-1.0, y);
    znorm = z->norm();
    *outStream << "\nNorm of (x - x) - 0: " << znorm << "\n";
    if ( std::abs(znorm) > errtol ) {
      *outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    };

    // set x to dim/2-th basis vector
    //z = x.basis(dim/2);
    //znorm = z->norm();
    //*outStream << "\nNorm of ROL::Vector z (basis vector): " << znorm << "\n";

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

