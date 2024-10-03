// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_01.cpp
    \brief Test ROL::ArrayFireVector interface.
*/

#include "ROL_ArrayFireVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Types.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>

typedef float RealT;
typedef float ElementT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  // *** Test body.

  try {

    /***** Display ArrayFire info. *****/
    af::dtype afType = f32;
    if (sizeof(RealT) == sizeof(double)) {
      afType = f64;
    }
    int device = 0;
    af::setDevice(device);
    std::string afinfo;
    afinfo = af::infoString(true);
    *outStream << std::endl << afinfo << std::endl;

    /***** Define ROL::Ptrs to random arrays. *****/
    ROL::Ptr<af::array> A_ptr = ROL::makePtr<af::array>(5,1,afType);
    ROL::Ptr<af::array> B_ptr = ROL::makePtr<af::array>(5,1,afType);
    ROL::Ptr<af::array> C_ptr = ROL::makePtr<af::array>(5,1,afType);
    *A_ptr = af::randu(A_ptr->dims(),afType);
    *B_ptr = af::randu(B_ptr->dims(),afType);
    *C_ptr = af::randu(C_ptr->dims(),afType);
    //af_print(*A_ptr);
    //af_print(*B_ptr);
    //af_print(*C_ptr);

    /***** Wrap ArrayFire ROL::Ptrs as ROL::ArrayFireVectors. *****/
    ROL::ArrayFireVector<RealT, ElementT> A(A_ptr);
    ROL::ArrayFireVector<RealT, ElementT> B(B_ptr);
    ROL::ArrayFireVector<RealT, ElementT> C(C_ptr);

    /***** Run checkVector tests. *****/
    // with random vector
    std::vector<RealT> consistency = A.checkVector(B, C, true, *outStream);
    ROL::StdVector<RealT, ElementT> checkvec(ROL::makePtrFromRef(consistency));
    if (checkvec.norm() > std::sqrt(ROL::ROL_EPSILON<RealT>())) {
      errorFlag++;
    }

    // with zero vectors
    RealT zero(0);
    A.scale(zero); B.scale(zero); C.scale(zero);
    consistency = A.checkVector(B, C, true, *outStream);
    if (checkvec.norm() > zero) {
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

