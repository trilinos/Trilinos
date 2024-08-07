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


#include "ROL_StdVector.hpp"
#include "ROL_RieszVector.hpp"
#include "ROL_RandomVector.hpp"

#include "ROL_DiagonalOperator.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

   

  typedef std::vector<RealT>    vector;
  typedef ROL::Vector<RealT>    V;
  typedef ROL::StdVector<RealT> SV;
  
  typedef ROL::LinearOperator<RealT>     LinearOperator;
  typedef ROL::DiagonalOperator<RealT>   DiagonalOperator;  

  typedef ROL::RieszPrimalVector<RealT>  PrimalVector;
  typedef ROL::RieszDualVector<RealT>    DualVector;

  
  Teuchos::GlobalMPISession mpiSession(&argc, &argv,0);

  int iprint = argc - 1;
  ROL::nullstream bhs; // outputs nothing
  ROL::Ptr<std::ostream> outStream;
  if  (iprint > 0) 
    outStream = ROL::makePtrFromRef(std::cout);
  else  
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag = 0;

  RealT errtol = ROL::ROL_THRESHOLD<RealT>();

  try {
    // Dimension of the optimization vector

    int dim = 10; 

    // Create Tpetra::MultiVectors (single vectors) 
    ROL::Ptr<vector> x_ptr = ROL::makePtr<vector>(dim);
    ROL::Ptr<vector> y_ptr = ROL::makePtr<vector>(dim);
    ROL::Ptr<vector> w_ptr = ROL::makePtr<vector>(dim,2.0);

    ROL::Ptr<V> xs = ROL::makePtr<SV>( x_ptr );
    ROL::Ptr<V> ys = ROL::makePtr<SV>( y_ptr );
    
    SV ws( w_ptr );

    ROL::RandomizeVector(*xs);
    ROL::RandomizeVector(*ys);

    ROL::Ptr<LinearOperator> W = ROL::makePtr<DiagonalOperator>(ws);
 
    PrimalVector x(xs,W);
    DualVector   y(ys,W);    

    RealT xy = x.dot(y.dual());
    RealT yx = y.dot(x.dual());

    *outStream << "\nAbsolute error between x.dot(y.dual()) and y.dot(x.dual()): "
              << std::abs(xy-yx) << "\n";
    *outStream << "x.dot(y.dual()): " << xy << "\n";
    *outStream << "y.dot(x.dual()): " << yx << "\n";
    if ( std::abs(xy-yx) > errtol ) {
      *outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    }

    RealT xx = std::sqrt(x.dot(x)), xnorm = x.norm();
    RealT yy = std::sqrt(y.dot(y)), ynorm = y.norm();

    *outStream << "\nAbsolute error between sqrt(x.dot(x)) and x.norm(): "
              << std::abs(xx-xnorm) << "\n";
    *outStream << "sqrt(x.dot(x)): " << xx << "\n";
    *outStream << "x.norm():       " << xnorm << "\n";
    if ( std::abs(xx-xnorm) > errtol ) {
      *outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    }

    *outStream << "\nAbsolute error between sqrt(y.dot(y)) and y.norm(): "
              << std::abs(yy-ynorm) << "\n";
    *outStream << "sqrt(y.dot(y)): " << yy << "\n";
    *outStream << "y.norm():       " << ynorm << "\n";
    if ( std::abs(yy-ynorm) > errtol ) {
      *outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    }

    // clone z from x, deep copy x into z, norm of z
    ROL::Ptr<ROL::Vector<RealT> > z = x.clone();
    z->set(x);
    RealT znorm = z->norm();
    *outStream << "\nNorm of ROL::Vector z (clone of x): " << znorm << "\n";
    if ( std::abs(xnorm - znorm) > errtol ) {
      *outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    }
    ROL::Ptr<ROL::Vector<RealT> > w = y.clone();
    w = y.clone();
    w->set(y);
    RealT wnorm = w->norm();
    *outStream << "\nNorm of ROL::Vector w (clone of y): " << wnorm << "\n";
    if ( std::abs(ynorm - wnorm) > errtol ) {
      *outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    }

    // Standard tests.
    ROL::Ptr<vector> x1_ptr = ROL::makePtr<vector>(dim);
    ROL::Ptr<vector> y1_ptr = ROL::makePtr<vector>(dim);
    ROL::Ptr<vector> z1_ptr = ROL::makePtr<vector>(dim);

    ROL::Ptr<V> x1s = ROL::makePtr<SV>( x1_ptr );
    ROL::Ptr<V> y1s = ROL::makePtr<SV>( y1_ptr );
    ROL::Ptr<V> z1s = ROL::makePtr<SV>( z1_ptr );

    ROL::RandomizeVector(*x1s);
    ROL::RandomizeVector(*y1s);
    ROL::RandomizeVector(*z1s);

    // Create ROL vectors
    PrimalVector x1(x1s,W);
    PrimalVector y1(y1s,W);
    PrimalVector z1(z1s,W);

    std::vector<RealT> consistency = x1.checkVector(y1, z1, true, *outStream);
    ROL::StdVector<RealT> checkvec( ROL::makePtrFromRef(consistency) );
    if (checkvec.norm() > std::sqrt(errtol)) {
      errorFlag++;
    }

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
