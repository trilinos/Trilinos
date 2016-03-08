// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER


/*! \file  test_01.cpp
    \brief Test ROL::ArrayFireVector interface.
*/

#include "ROL_ArrayFireVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Types.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>

typedef float RealT;
typedef float ElementT;

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

    /***** Define RCPs to random arrays. *****/
    Teuchos::RCP<af::array> A_rcp = Teuchos::rcp( new af::array(5,1,afType));
    Teuchos::RCP<af::array> B_rcp = Teuchos::rcp( new af::array(5,1,afType));
    Teuchos::RCP<af::array> C_rcp = Teuchos::rcp( new af::array(5,1,afType));
    *A_rcp = af::randu(A_rcp->dims(),afType);
    *B_rcp = af::randu(B_rcp->dims(),afType);
    *C_rcp = af::randu(C_rcp->dims(),afType);
    //af_print(*A_rcp);
    //af_print(*B_rcp);
    //af_print(*C_rcp);

    /***** Wrap ArrayFire RCPs as ROL::ArrayFireVectors. *****/
    ROL::ArrayFireVector<RealT, ElementT> A(A_rcp);
    ROL::ArrayFireVector<RealT, ElementT> B(B_rcp);
    ROL::ArrayFireVector<RealT, ElementT> C(C_rcp);

    /***** Run checkVector tests. *****/
    // with random vector
    std::vector<RealT> consistency = A.checkVector(B, C, true, *outStream);
    ROL::StdVector<RealT, ElementT> checkvec(Teuchos::rcpFromRef(consistency));
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

