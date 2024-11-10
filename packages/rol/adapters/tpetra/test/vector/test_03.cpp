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


#include "ROL_ScaledTpetraMultiVector.hpp"
#include "ROL_StdVector.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Tpetra_Core.hpp"

typedef double RealT;
typedef double ElementT;

typedef Tpetra::Map<>::local_ordinal_type LO;
typedef Tpetra::Map<>::global_ordinal_type GO;
typedef Tpetra::Map<>::node_type Node;
typedef Tpetra::Map<LO, GO, Node> Map;
typedef Tpetra::MultiVector<RealT, LO, GO, Node> MV;
typedef Tpetra::Vector<RealT, LO, GO, Node> V;
typedef ROL::Ptr<MV> MVP;
typedef ROL::Ptr<V> VP;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv,0);
  ROL::Ptr<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int iprint = argc - 1;
  ROL::nullstream bhs; // outputs nothing
  std::ostream& outStream = (iprint > 0) ? std::cout : bhs;

  int errorFlag = 0;

  RealT errtol = ROL::ROL_THRESHOLD<RealT>();

  try {
    // Dimension of the optimization vector

    int dim = 10; 
  
    ROL::Ptr<Map> map = ROL::makePtr<Map>(dim,0,comm);

    // Create Tpetra::MultiVectors (single vectors) 
    MVP x_ptr = ROL::makePtr<MV>(map,1,true); 
    MVP y_ptr = ROL::makePtr<MV>(map,1,true); 
    VP W_ptr = ROL::makePtr<V>(map,true);

    // Random elements
    //x_ptr->randomize();
    //y_ptr->randomize();
    x_ptr->putScalar(1.0);
    y_ptr->putScalar(1.0);

    // Set all values to 2
    W_ptr->putScalar(2.0);

    // Create ROL vectors
    ROL::PrimalScaledTpetraMultiVector<RealT,LO,GO,Node> x(x_ptr,W_ptr);
    ROL::DualScaledTpetraMultiVector<RealT,LO,GO,Node>   y(y_ptr,W_ptr);

//    const ROL::Vector<RealT> &g = x.dual();
//    const ROL::Vector<RealT> &h = x.dual();
//    RealT hnorm = h.norm();
//    RealT gnorm = g.norm();

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
    // Create Tpetra::MultiVectors (single vectors) 
    MVP x1_ptr = ROL::makePtr<MV>(map,1,true); 
    MVP y1_ptr = ROL::makePtr<MV>(map,1,true); 
    MVP z1_ptr = ROL::makePtr<MV>(map,1,true); 
    ROL::PrimalScaledTpetraMultiVector<RealT,LO,GO,Node> x1(x1_ptr,W_ptr);
    ROL::PrimalScaledTpetraMultiVector<RealT,LO,GO,Node> y1(y1_ptr,W_ptr);
    ROL::PrimalScaledTpetraMultiVector<RealT,LO,GO,Node> z1(z1_ptr,W_ptr);
    x1_ptr->randomize();
    y1_ptr->randomize();
    z1_ptr->randomize();

    std::vector<RealT> consistency = x1.checkVector(y1, z1, true, outStream);
    ROL::StdVector<RealT> checkvec(ROL::makePtrFromRef(consistency));
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
