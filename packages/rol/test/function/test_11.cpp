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

/*! \file  test_11.cpp
    \brief Checks that Coleman-Li BoundConstraint functions agree across implementations.

*/

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_Bounds.hpp"
#include "ROL_StdBoundConstraint.hpp"
#include "ROL_StdVector.hpp"

typedef double RealT;

RealT calcError(ROL::Vector<RealT> &a, const ROL::Vector<RealT> &b) {
  RealT one(1);
  a.axpy(-one, b);
  a.applyUnary(ROL::Elementwise::AbsoluteValue<RealT>());
  return a.reduce(ROL::Elementwise::ReductionMax<RealT>());
}

int testRandomInputs(int numPoints, RealT tol, ROL::Ptr<std::ostream> outStream) {

  RealT invInftyNorm, jacInftyNorm;

  // Generate standard vectors that hold data.
  ROL::Ptr<std::vector<RealT>> vp = ROL::makePtr<std::vector<RealT>>(numPoints, 0.0);
  ROL::Ptr<std::vector<RealT>> xp = ROL::makePtr<std::vector<RealT>>(numPoints, 0.0);
  ROL::Ptr<std::vector<RealT>> gp = ROL::makePtr<std::vector<RealT>>(numPoints, 0.0);
  ROL::Ptr<std::vector<RealT>> lp = ROL::makePtr<std::vector<RealT>>(numPoints, 0.0);
  ROL::Ptr<std::vector<RealT>> up = ROL::makePtr<std::vector<RealT>>(numPoints, 0.0);
  // Include space for storing results.
  ROL::Ptr<std::vector<RealT>> result1p
    = ROL::makePtr<std::vector<RealT>>(numPoints, 0.0);
  ROL::Ptr<std::vector<RealT>> result2p
    = ROL::makePtr<std::vector<RealT>>(numPoints, 0.0);

  // Use the standard vectors above to define ROL::StdVectors (or, in the
  // case of l and u, pointers to ROL::Vectors).
  ROL::StdVector<RealT> v(vp);
  ROL::StdVector<RealT> x(xp);
  ROL::StdVector<RealT> g(gp);
  ROL::Ptr<ROL::Vector<RealT>> l = ROL::makePtr<ROL::StdVector<RealT>>(lp);
  ROL::Ptr<ROL::Vector<RealT>> u = ROL::makePtr<ROL::StdVector<RealT>>(up);
  ROL::StdVector<RealT> result1(result1p);
  ROL::StdVector<RealT> result2(result2p);

  // Initialize.
  v.setScalar(  2.0);
  g.randomize(-20.0, 20.0);
  x.randomize(  0.0, 10.0);
  l->setScalar( 0.0);
  u->randomize( 0.0, 10.0);
  u->plus(x);

  ROL::StdBoundConstraint<RealT> standardVecBC(*lp,*up);
  ROL::Bounds<RealT>             elementwiseBC( l,  u );

  standardVecBC.applyInverseScalingFunction(result1, v, x, g);
  elementwiseBC.applyInverseScalingFunction(result2, v, x, g);
  invInftyNorm = calcError(result1, result2);

  standardVecBC.applyScalingFunctionJacobian(result1, v, x, g);
  elementwiseBC.applyScalingFunctionJacobian(result2, v, x, g);
  jacInftyNorm = calcError(result1, result2);

  *outStream << std::endl;
  *outStream << "|StdBoundConstraint - Bounds| at " << numPoints
             << " Randomly Sampled Points (Infinity Norm): " << std::endl
             << "  Inverse          = " << invInftyNorm << std::endl
             << "  Jacobian         = " << jacInftyNorm << std::endl;
  *outStream << std::endl;

  return (invInftyNorm + jacInftyNorm) > tol;
}

int testCases(RealT tol, ROL::Ptr<std::ostream> outStream) {

  int numCases = 3;

  std::vector<RealT> ewErrors, svErrors;

  // Generate standard vectors that hold data.
  ROL::Ptr<std::vector<RealT>> vp = ROL::makePtr<std::vector<RealT>>(numCases, 0.0);
  ROL::Ptr<std::vector<RealT>> xp = ROL::makePtr<std::vector<RealT>>(numCases, 0.0);
  ROL::Ptr<std::vector<RealT>> gp = ROL::makePtr<std::vector<RealT>>(numCases, 0.0);
  ROL::Ptr<std::vector<RealT>> lp = ROL::makePtr<std::vector<RealT>>(numCases, 0.0);
  ROL::Ptr<std::vector<RealT>> up = ROL::makePtr<std::vector<RealT>>(numCases, 0.0);
  // Include space for storing results.
  ROL::Ptr<std::vector<RealT>> resultp
    = ROL::makePtr<std::vector<RealT>>(numCases, 0.0);
  ROL::Ptr<std::vector<RealT>> targetp
    = ROL::makePtr<std::vector<RealT>>(numCases, 0.0);

  // Use the standard vectors above to define ROL::StdVectors (or, in the
  // case of l and u, pointers to ROL::Vectors).
  ROL::StdVector<RealT> v(vp);
  ROL::StdVector<RealT> x(xp);
  ROL::StdVector<RealT> g(gp);
  ROL::Ptr<ROL::Vector<RealT>> l = ROL::makePtr<ROL::StdVector<RealT>>(lp);
  ROL::Ptr<ROL::Vector<RealT>> u = ROL::makePtr<ROL::StdVector<RealT>>(up);
  ROL::StdVector<RealT> result(resultp);
  ROL::StdVector<RealT> target(targetp);

  // Problem 1
  (*vp)[0] =  4.0;
  (*xp)[0] =  1.9;
  (*gp)[0] =  0.5;
  (*lp)[0] =  0.0;
  (*up)[0] =  2.0;

  // Problem 2
  (*vp)[1] = -1.0;
  (*xp)[1] = 10.0;
  (*gp)[1] =  0.0002;
  (*lp)[1] = ROL::ROL_NINF<RealT>();
  (*up)[1] = ROL::ROL_INF<RealT>();

  // Problem 3
  (*vp)[2] = -0.0002;
  (*xp)[2] =  1.0;
  (*gp)[2] =  0.5;
  (*lp)[2] =  0.0;
  (*up)[2] = ROL::ROL_INF<RealT>();

  ROL::StdBoundConstraint<RealT> standardVecBC(*lp,*up);
  ROL::Bounds<RealT>             elementwiseBC( l,  u );

  // Expected results when applying the scaling function to v.
  (*targetp)[0] = (*vp)[0]*(*gp)[0];
  (*targetp)[1] = (*vp)[1]*1.0;
  (*targetp)[2] = (*vp)[2]*1.0;

  for (unsigned long i = 0; i < targetp->size(); i++) {
    (*targetp)[i] = (*vp)[i]*(*vp)[i]/(*targetp)[i];
  }
  standardVecBC.applyInverseScalingFunction(result, v, x, g);
  svErrors.push_back(calcError(result, target));
  elementwiseBC.applyInverseScalingFunction(result, v, x, g);
  ewErrors.push_back(calcError(result, target));

  (*targetp)[0] = (*vp)[0]*(*gp)[0];
  (*targetp)[1] = 0.0;
  (*targetp)[2] = 0.0;
  standardVecBC.applyScalingFunctionJacobian(result, v, x, g);
  svErrors.push_back(calcError(result, target));
  elementwiseBC.applyScalingFunctionJacobian(result, v, x, g);
  ewErrors.push_back(calcError(result, target));

  *outStream << "StdBoundConstraint Test Case Errors (Infinity Norm):" << std::endl
    << "  Inverse          = " << svErrors[1] << std::endl
    << "  Jacobian         = " << svErrors[2] << std::endl;
  *outStream << "Bounds             Test Case Errors (Infinity Norm):" << std::endl
    << "  Inverse          = " << ewErrors[1] << std::endl
    << "  Jacobian         = " << ewErrors[2] << std::endl;
  *outStream << std::endl;

  RealT maxError = std::max(*std::max_element(svErrors.begin(), svErrors.end()),
                            *std::max_element(ewErrors.begin(), ewErrors.end()));
  return maxError > tol;
}

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a
  // (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag = 0;

  // *** Test body.

  try {
    RealT tol = 1e-8;  // tolerance
    int   n   = 1e+3;  // number of random test points
    errorFlag += testRandomInputs(n, tol, outStream);
    errorFlag += testCases(tol, outStream);
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
