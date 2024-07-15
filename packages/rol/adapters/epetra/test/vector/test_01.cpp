// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_02.cpp
    \brief Test Epetra_MultiVector interface.
*/

#include "ROL_EpetraMultiVector.hpp"
#include "ROL_Types.hpp"
#include "Epetra_Map.h"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <iostream>

typedef double RealT;
typedef double ElementT;

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  Teuchos::GlobalMPISession mpiSession(&argc, &argv,0);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  RealT errtol = ROL::ROL_THRESHOLD<RealT>();

  // *** Test body.

  try {

    int dim = 100;
    Epetra_Map Map(dim, 0, Comm);

    ROL::Ptr<Epetra_MultiVector> x_ptr = ROL::makePtr<Epetra_MultiVector>(Map, 1);
    ROL::Ptr<Epetra_MultiVector> y_ptr = ROL::makePtr<Epetra_MultiVector>(Map, 1);
    ROL::EpetraMultiVector<RealT> x(x_ptr);
    ROL::EpetraMultiVector<RealT> y(y_ptr);

    // set x,y
    for (int i=0; i<dim; i++) {
      x_ptr->Random();
      y_ptr->PutScalar(2.0);
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
    ROL::Ptr<ROL::Vector<RealT> > z = x.clone();
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

    // set x to first basis vector
    z = x.basis(0);
    znorm = z->norm();
    *outStream << "\nNorm of ROL::Vector z (first basis vector): " << znorm << "\n";
    if ( std::abs(znorm-1.0) > errtol ) {
      *outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    };
    // set x to middle basis vector
    z = x.basis(dim/2);
    znorm = z->norm();
    *outStream << "\nNorm of ROL::Vector z ('middle' basis vector): " << znorm << "\n";
    if ( std::abs(znorm-1.0) > errtol ) {
      *outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    };
    // set x to last basis vector
    z = x.basis(dim-1);
    znorm = z->norm();
    *outStream << "\nNorm of ROL::Vector z (last basis vector): " << znorm << "\n";
    if ( std::abs(znorm-1.0) > errtol ) {
      *outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    };

  }
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;

}

