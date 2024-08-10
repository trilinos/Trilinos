// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_01.cpp
    \brief Test Tpetra_MultiVector interface.
*/


#include "ROL_StdVector.hpp"
#include "ROL_ScaledThyraVector.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_VectorSpaceBase_decl.hpp"

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

    int dim = 100; 
  
    Teuchos::RCP<Thyra::VectorSpaceBase<RealT> > map = Thyra::defaultSpmdVectorSpace<RealT>(dim);   

    // Create Teuchos::RCPs to Thyra::Vectors
    Teuchos::RCP<Thyra::VectorBase<RealT> > x_ptr = Thyra::createMember<RealT>(map);
    Teuchos::RCP<Thyra::VectorBase<RealT> > y_ptr = Thyra::createMember<RealT>(map);
    Teuchos::RCP<Thyra::VectorBase<RealT> > W_ptr = Thyra::createMember<RealT>(map);

    // Random elements
    //x_ptr->randomize();
    //y_ptr->randomize();
    Thyra::assign(x_ptr.ptr(),1.0);
    Thyra::assign(y_ptr.ptr(),1.0);

    // Set all values to 2
    Thyra::assign(W_ptr.ptr(),2.0);

    // Create ROL vectors
    ROL::PrimalScaledThyraVector<RealT> x(x_ptr,W_ptr);
    ROL::DualScaledThyraVector<RealT>   y(y_ptr,W_ptr);

    RealT xy = x.dot(y.dual());
    RealT yx = y.dot(x.dual());

    outStream << "\nAbsolute error between x.dot(y.dual()) and y.dot(x.dual()): "
              << std::abs(xy-yx) << "\n";
    outStream << "x.dot(y.dual()): " << xy << "\n";
    outStream << "y.dot(x.dual()): " << yx << "\n";
    if ( std::abs(xy-yx) > errtol ) {
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
    Teuchos::RCP<ROL::Vector<RealT> > z = x.clone();
    z->set(x);
    RealT znorm = z->norm();
    outStream << "\nNorm of ROL::Vector z (clone of x): " << znorm << "\n";
    if ( std::abs(xnorm - znorm) > errtol ) {
      outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    }
    Teuchos::RCP<ROL::Vector<RealT> > w = y.clone();
    w = y.clone();
    w->set(y);
    RealT wnorm = w->norm();
    outStream << "\nNorm of ROL::Vector w (clone of y): " << wnorm << "\n";
    if ( std::abs(ynorm - wnorm) > errtol ) {
      outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    }

    // Standard tests.
    // Create Tpetra::MultiVectors (single vectors) 
    Teuchos::RCP<Thyra::VectorBase<RealT> > x1_ptr = Thyra::createMember<RealT>(map);
    Teuchos::RCP<Thyra::VectorBase<RealT> > y1_ptr = Thyra::createMember<RealT>(map);
    Teuchos::RCP<Thyra::VectorBase<RealT> > z1_ptr = Thyra::createMember<RealT>(map);
    ROL::PrimalScaledThyraVector<RealT> x1(x1_ptr,W_ptr);
    ROL::PrimalScaledThyraVector<RealT> y1(y1_ptr,W_ptr);
    ROL::PrimalScaledThyraVector<RealT> z1(z1_ptr,W_ptr);
    Thyra::randomize(-1.0,1.0,x1_ptr.ptr());
    Thyra::randomize(-1.0,1.0,y1_ptr.ptr());
    Thyra::randomize(-1.0,1.0,z1_ptr.ptr());

    std::vector<RealT> consistency = x1.checkVector(y1, z1, true, outStream);
    ROL::StdVector<RealT> checkvec(Teuchos::rcpFromRef(consistency));
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
