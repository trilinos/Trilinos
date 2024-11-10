// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "ROL_StdVector.hpp"
#include "ROL_ThyraVector.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_Types.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_VectorSpaceBase_decl.hpp"

typedef double RealT;
typedef double ElementT;

int main(int argc, char *argv[]) {

  using namespace Teuchos;

  GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcpFromRef(std::cout);
  else
    outStream = Teuchos::rcpFromRef(bhs);

  int errorFlag  = 0;

  RealT errtol = ROL::ROL_THRESHOLD<RealT>();

  // *** Test body.

  try {

    int dim = 100; 

    Teuchos::RCP<Thyra::VectorSpaceBase<RealT> > euclidean = Thyra::defaultSpmdVectorSpace<RealT>(dim);   

    // Create Teuchos::RCPs to Thyra::Vectors
    Teuchos::RCP<Thyra::VectorBase<RealT> > x_ptr = Thyra::createMember<RealT>(euclidean);
    Teuchos::RCP<Thyra::VectorBase<RealT> > y_ptr = Thyra::createMember<RealT>(euclidean);
    Teuchos::RCP<Thyra::VectorBase<RealT> > z_ptr = Thyra::createMember<RealT>(euclidean);

    // Create ROL::ThyraVectors
  
    ROL::ThyraVector<RealT> x(x_ptr);
    ROL::ThyraVector<RealT> y(y_ptr);
    ROL::ThyraVector<RealT> z(z_ptr);

    RealT left = -1e0, right = 1e0;
     
    // Randomize x,y,z vectors
    ROL::RandomizeVector(x,left,right);
    ROL::RandomizeVector(y,left,right);
    ROL::RandomizeVector(z,left,right);

    // Standard tests.
    std::vector<RealT> consistency = x.checkVector(y, z, true, *outStream);
    ROL::StdVector<RealT, ElementT> checkvec(Teuchos::rcpFromRef(consistency));
    if (checkvec.norm() > std::sqrt(ROL::ROL_EPSILON<RealT>())) {
      errorFlag++;
    }      

    // Basis tests.
    // set x to first basis vector
    Teuchos::RCP<ROL::Vector<RealT> > zp = x.clone();
    zp = x.basis(0);
    RealT znorm = zp->norm();
    *outStream << "Norm of ROL::Vector z (first basis vector): " << znorm << "\n";
    if ( std::abs(znorm-1.0) > errtol ) {
      *outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    };
    // set x to middle basis vector
    zp = x.basis(dim/2);
    znorm = zp->norm();
    *outStream << "\nNorm of ROL::Vector z ('middle' basis vector): " << znorm << "\n";
    if ( std::abs(znorm-1.0) > errtol ) {
      *outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    };
    // set x to last basis vector
    zp = x.basis(dim-1);
    znorm = zp->norm();
    *outStream << "\nNorm of ROL::Vector z (last basis vector): " << znorm << "\n";
    if ( std::abs(znorm-1.0) > errtol ) {
      *outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    };

    // Repeat the checkVector tests with a zero vector.
    x.scale(0.0);
    consistency = x.checkVector(x, x, true, *outStream);
    if (checkvec.norm() > 0.0) {
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



