// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_05.cpp
    \brief Test ScaledStdVector interface.
*/


#include "ROL_ScaledStdVector.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

typedef double RealT;
typedef double ElementT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv,0);

  int iprint = argc - 1;
  ROL::nullstream bhs; // outputs nothing
  std::ostream& outStream = (iprint > 0) ? std::cout : bhs;

  int errorFlag = 0;

  RealT errtol = ROL::ROL_THRESHOLD<RealT>();

  try {
    // Dimension of the optimization vector

    int dim = 10; 

    // Create Tpetra::MultiVectors (single vectors) 
    ROL::Ptr<std::vector<ElementT> > x_ptr
      = ROL::makePtr<std::vector<ElementT>>(dim); 
    ROL::Ptr<std::vector<ElementT> > y_ptr
      = ROL::makePtr<std::vector<ElementT>>(dim); 
    ROL::Ptr<std::vector<ElementT> > W_ptr
      = ROL::makePtr<std::vector<ElementT>>(dim,static_cast<ElementT>(2)); 

    // Random elements
    for (int i = 0; i < dim; i++) {
      (*x_ptr)[i] = static_cast<ElementT>(rand())/static_cast<ElementT>(RAND_MAX);
      (*y_ptr)[i] = static_cast<ElementT>(rand())/static_cast<ElementT>(RAND_MAX);
    }

    // Create ROL vectors
    ROL::PrimalScaledStdVector<RealT,ElementT> x(x_ptr,W_ptr);
    ROL::DualScaledStdVector<RealT,ElementT>   y(y_ptr,W_ptr);

    RealT xy = x.dot(y.dual());
    RealT yx = y.dot(x.dual());
    RealT axy = x.apply(y);
    RealT ayx = y.apply(x);

    outStream << "\nAbsolute error between x.dot(y.dual()) and y.dot(x.dual()): "
              << std::abs(xy-yx) << "\n";
    outStream << "x.dot(y.dual()): " << xy << "\n";
    outStream << "y.dot(x.dual()): " << yx << "\n";
    if ( std::abs(xy-yx) > errtol ) {
      outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    }

    outStream << "\nAbsolute error between x.apply(y) and y.apply(x): "
              << std::abs(axy-ayx) << "\n";
    outStream << "x.apply(y): " << axy << "\n";
    outStream << "y.apply(x): " << ayx << "\n";
    if ( std::abs(axy-ayx) > errtol ) {
      outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    }

    outStream << "\nAbsolute error between x.apply(y) and x.dot(y.dual()): "
              << std::abs(axy-ayx) << "\n";
    outStream << "x.apply(y):      " << axy << "\n";
    outStream << "x.dot(y.dual()): " << xy  << "\n";
    if ( std::abs(axy-xy) > errtol ) {
      outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    }

    RealT xx = std::sqrt(x.dot(x)), xnorm = x.norm();
    RealT yy = std::sqrt(y.dot(y)), ynorm = y.norm();

    outStream << "\nAbsolute error between sqrt(x.dot(x)) and x.norm(): "
              << std::abs(xx-xnorm) << "\n";
    outStream << "sqrt(x.dot(x)): " << xx << "\n";
    outStream << "x.norm():       " << xnorm << "\n";
    if ( std::abs(xx-xnorm) > errtol ) {
      outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    }

    outStream << "\nAbsolute error between sqrt(y.dot(y)) and y.norm(): "
              << std::abs(yy-ynorm) << "\n";
    outStream << "sqrt(y.dot(y)): " << yy << "\n";
    outStream << "y.norm():       " << ynorm << "\n";
    if ( std::abs(yy-ynorm) > errtol ) {
      outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    }

    // clone z from x, deep copy x into z, norm of z
    ROL::Ptr<ROL::Vector<RealT> > z = x.clone();
    z->set(x);
    RealT znorm = z->norm();
    outStream << "\nNorm of ROL::Vector z (clone of x): " << znorm << "\n";
    if ( std::abs(xnorm - znorm) > errtol ) {
      outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    }
    ROL::Ptr<ROL::Vector<RealT> > w = y.clone();
    w = y.clone();
    w->set(y);
    RealT wnorm = w->norm();
    outStream << "\nNorm of ROL::Vector w (clone of y): " << wnorm << "\n";
    if ( std::abs(ynorm - wnorm) > errtol ) {
      outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    }

    // Standard tests.
    ROL::Ptr<std::vector<ElementT> > x1_ptr
      = ROL::makePtr<std::vector<ElementT>>(dim);
    ROL::Ptr<std::vector<ElementT> > y1_ptr
      = ROL::makePtr<std::vector<ElementT>>(dim);
    ROL::Ptr<std::vector<ElementT> > z1_ptr
      = ROL::makePtr<std::vector<ElementT>>(dim);

    // Random elements
    for (int i = 0; i < dim; i++) {
      (*x1_ptr)[i] = static_cast<ElementT>(rand())/static_cast<ElementT>(RAND_MAX);
      (*y1_ptr)[i] = static_cast<ElementT>(rand())/static_cast<ElementT>(RAND_MAX);
      (*z1_ptr)[i] = static_cast<ElementT>(rand())/static_cast<ElementT>(RAND_MAX);
    }

    // Create ROL vectors
    ROL::PrimalScaledStdVector<RealT,ElementT> x1(x1_ptr,W_ptr);
    ROL::PrimalScaledStdVector<RealT,ElementT> y1(y1_ptr,W_ptr);
    ROL::PrimalScaledStdVector<RealT,ElementT> z1(z1_ptr,W_ptr);

    std::vector<RealT> consistency = x1.checkVector(y1, z1, true, outStream);
    ROL::StdVector<RealT> checkvec( ROL::makePtrFromRef(consistency) );
    if (checkvec.norm() > std::sqrt(errtol)) {
      errorFlag++;
    }

  }

  catch (std::logic_error& err) {
    outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
