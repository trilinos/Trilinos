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


/*! \file  stdvector.cpp
    \brief Test std::vector interface.
*/

#include "ROL2.hpp"

int main(int argc, char *argv[]) {

  using RealT = double;

  auto os_ptr = ROL2::makeStreamPtr(std::cout, argc);
  auto& os = *os_ptr;

  int errorFlag  = 0;

  auto errtol = ROL2::default_tolerance<RealT>();

  try {

    // Number of elements in StdVector objects
    int dim = 100; 

    ROL2::StdVector<RealT> x(dim), y(dim), z(dim);

    // Lower and upper bounds on randomized vector element values
    RealT left = -1e0, right = 1e0;

    // set x,y,z
    x.randomize(left,right);
    y.randomize(left,right);
    z.randomize(left,right);

    // Standard tests.
    int consistencyError = x.checkVector(y, z, true, os);

    errorFlag += consistencyError;

    // Basis tests.
    // set x to first basis vector
    auto zp = x.basis(0);
    RealT znorm = zp->norm();
    os << "\nNorm of ROL::Vector z (first basis vector): " << znorm << "\n";
    if ( std::abs(znorm-1.0) > errtol ) {
      os << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    };

    // set x to middle basis vector
    zp = x.basis(dim/2);
    znorm = zp->norm();
    os << "\nNorm of ROL::Vector z ('middle' basis vector): " << znorm << "\n";

    if ( std::abs(znorm-1.0) > errtol ) {
      os << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    };

    // set x to last basis vector
    zp = x.basis(dim-1);
    znorm = zp->norm();
    os << "\nNorm of ROL::Vector z (last basis vector): " << znorm << "\n";
    if ( std::abs(znorm-1.0) > errtol ) {
      os << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    };

    os << "\n";

    // Repeat the checkVector tests with a zero vector.
    x.scale(0.0);
    consistencyError = x.checkVector(x, x, true, os);
    errorFlag += consistencyError;
  }
  catch (std::logic_error& err) {
    os << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0) std::cout << "End Result: TEST FAILED\n";
  else                std::cout << "End Result: TEST PASSED\n";

  return 0;

}

