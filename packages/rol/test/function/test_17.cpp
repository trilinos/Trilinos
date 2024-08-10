// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_17.cpp
    \brief Checks that Coleman-Li BoundConstraint functions agree across implementations.

*/

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_Bounds.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_StdBoundConstraint.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_UnaryFunctions.hpp"

typedef double RealT;

RealT calcError(ROL::Vector<RealT> &a, const ROL::Vector<RealT> &b) {
  RealT one(1);
  a.axpy(-one, b);
  a.applyUnary(ROL::Elementwise::AbsoluteValue<RealT>());
  return a.reduce(ROL::Elementwise::ReductionMax<RealT>());
}

int testRandomInputs(int numPoints, RealT tol,
                     ROL::Ptr<std::ostream> outStream) {

  int   errorFlag = 0;
  RealT errorInftyNorm;

  // Generate standard vectors that hold data.
  ROL::Ptr<std::vector<RealT>> vsv
    = ROL::makePtr<std::vector<RealT>>(numPoints, 0.0);
  ROL::Ptr<std::vector<RealT>> xsv
    = ROL::makePtr<std::vector<RealT>>(numPoints, 0.0);
  ROL::Ptr<std::vector<RealT>> gsv
    = ROL::makePtr<std::vector<RealT>>(numPoints, 0.0);
  ROL::Ptr<std::vector<RealT>> lsv
    = ROL::makePtr<std::vector<RealT>>(numPoints, 0.0);
  ROL::Ptr<std::vector<RealT>> usv
    = ROL::makePtr<std::vector<RealT>>(numPoints, 0.0);
  // Include space for storing results.
  ROL::Ptr<std::vector<RealT>> out1sv
    = ROL::makePtr<std::vector<RealT>>(numPoints, 0.0);
  ROL::Ptr<std::vector<RealT>> out2sv
    = ROL::makePtr<std::vector<RealT>>(numPoints, 0.0);

  // Use these standard vectors to define ROL::StdVectors (or, in the case of lp
  // and up, pointers to ROL::Vectors).
  ROL::StdVector<RealT> v(vsv);
  ROL::StdVector<RealT> x(xsv);
  ROL::StdVector<RealT> g(gsv);
  ROL::Ptr<ROL::Vector<RealT>> lp = ROL::makePtr<ROL::StdVector<RealT>>(lsv);
  ROL::Ptr<ROL::Vector<RealT>> up = ROL::makePtr<ROL::StdVector<RealT>>(usv);
  ROL::StdVector<RealT> out1(out1sv);
  ROL::StdVector<RealT> out2(out2sv);

  // Initialize.
  lp->randomize(-10.0, 10.0);
  up->randomize(  0.0, 10.0);
  up->plus(*lp);
  x.randomize(-20.0, 20.0);
  v.randomize(- 3.0,  3.0);
  g.randomize(-20.0, 20.0);

  ROL::Bounds<RealT>             boundsBC( lp , up );
  ROL::StdBoundConstraint<RealT> stdvecBC(*lsv,*usv);
  boundsBC.projectInterior(x);

  // Test 1 - Check that the Elementwise applyInverseScalingFunction does
  //          indeed apply a scale factor to the vector v.

  boundsBC.applyInverseScalingFunction(out1, v, x, g);
  out1.applyUnary(ROL::Elementwise::Reciprocal<RealT>());
  boundsBC.applyInverseScalingFunction(out2, out1, x, g);
  out2.applyUnary(ROL::Elementwise::Reciprocal<RealT>());
  errorInftyNorm = calcError(out2, v);
  errorFlag += errorInftyNorm > tol;

  *outStream << std::endl;
  *outStream << "Scaling Check at " << numPoints
             << " Randomly Sampled Points -- " << std::endl
             << " Infinity Norm of | v - 1/f(1/f(v)) | = "
             << errorInftyNorm << std::endl;
  *outStream << std::endl;

  // Test 2 - Use finite differences to check that the Elementwise
  //          applyScalingFunctionJacobian and applyInverseScalingFunction are
  //          consistent with each other.
  //          This test is meant to be visually inspected; it cannot cause the
  //          cpp file to fail.

  class TestWrapper : public ROL::Constraint<RealT> {
   private:
    ROL::Bounds<RealT> boundsBC_;
    ROL::StdVector<RealT> g_;
    ROL::StdVector<RealT> v_;

   public:
    TestWrapper(ROL::Bounds<RealT>& boundsBC, ROL::StdVector<RealT>& g)
      : boundsBC_(boundsBC), g_(g), v_(g) {
        RealT one(1);
        v_.setScalar(one);
      }

    void value(ROL::Vector<RealT>& c, const ROL::Vector<RealT>& x,
               RealT& tol) override {
      boundsBC_.applyInverseScalingFunction(c, v_, x, g_);
      c.applyUnary( ROL::Elementwise::Reciprocal<RealT>());
      c.applyBinary(ROL::Elementwise::Multiply<RealT>(), g_);
    }

    void applyJacobian(ROL::Vector<RealT>& jv, const ROL::Vector<RealT>& v,
                       const ROL::Vector<RealT>& x, RealT& tol) override {
      boundsBC_.applyScalingFunctionJacobian(jv, v, x, g_);
    }
  } testWrapper(boundsBC, g);

  // Use out1 and and out2 as working arrays to build a point comfortably
  // within our bounds (so that the finite difference increments don't all step
  // out). Larger values of gamma => larger separation between this point
  // and our bounds.
  RealT gamma = 1e-8;
  out1.randomize(-1.0, 1.0);
  out2.set(*up);
  out2.axpy(-1,*lp);
  out2.scale(1 - gamma);
  out1.applyBinary(ROL::Elementwise::Multiply<RealT>(), out2);
  out1.plus(*lp);
  out1.plus(*up);
  out1.scale(0.5);  // the point at which we check the Jacobian

  *outStream << "Elementwise Jacobian Check:" << std::endl;
  testWrapper.checkApplyJacobian(out1, v, out2, true, *outStream, 15);
  *outStream << std::endl;

  // Test 3 - Check that applyInverseScalingFunction and
  //          applyScalingFunctionJacobian agree between the Elementwise and
  //          StdVector implementations.

  boundsBC.applyInverseScalingFunction(out1, v, x, g);
  stdvecBC.applyInverseScalingFunction(out2, v, x, g);
  errorInftyNorm = 100;
  errorInftyNorm = calcError(out1, out2);
  errorFlag += errorInftyNorm > tol;

  *outStream << "Consistency Check at " << numPoints
             << " Randomly Sampled Points -- " << std::endl
             << " Infinity Norm of | StdBoundConstraint - Elementwise |:"
             << std::endl
             << "  Inverse  = " << errorInftyNorm << std::endl;

  boundsBC.applyScalingFunctionJacobian(out1, v, x, g);
  stdvecBC.applyScalingFunctionJacobian(out2, v, x, g);
  errorInftyNorm = calcError(out1, out2);
  errorFlag += errorInftyNorm > tol;

  *outStream << "  Jacobian = " << errorInftyNorm << std::endl;
  *outStream << std::endl;

  return errorFlag;
}

int testCases(RealT tol, ROL::Ptr<std::ostream> outStream) {

  // Test 4 - Check the Elementwise and StdVector implementations of
  //          applyInverseScalingFunction and applyScalingFunctionJacobian on
  //          specific test cases.

  int numCases = 3;

  std::vector<RealT> ewErrors, svErrors;

  // Generate standard vectors that hold data.
  ROL::Ptr<std::vector<RealT>> vsv
    = ROL::makePtr<std::vector<RealT>>(numCases, 0.0);
  ROL::Ptr<std::vector<RealT>> xsv
    = ROL::makePtr<std::vector<RealT>>(numCases, 0.0);
  ROL::Ptr<std::vector<RealT>> gsv
    = ROL::makePtr<std::vector<RealT>>(numCases, 0.0);
  ROL::Ptr<std::vector<RealT>> lsv
    = ROL::makePtr<std::vector<RealT>>(numCases, 0.0);
  ROL::Ptr<std::vector<RealT>> usv
    = ROL::makePtr<std::vector<RealT>>(numCases, 0.0);
  // Include space for storing results.
  ROL::Ptr<std::vector<RealT>> resultp
    = ROL::makePtr<std::vector<RealT>>(numCases, 0.0);
  ROL::Ptr<std::vector<RealT>> targetp
    = ROL::makePtr<std::vector<RealT>>(numCases, 0.0);

  // Use the standard vectors above to define ROL::StdVectors (or, in the
  // case of lp and up, pointers to ROL::Vectors).
  ROL::StdVector<RealT> v(vsv);
  ROL::StdVector<RealT> x(xsv);
  ROL::StdVector<RealT> g(gsv);
  ROL::Ptr<ROL::Vector<RealT>> lp = ROL::makePtr<ROL::StdVector<RealT>>(lsv);
  ROL::Ptr<ROL::Vector<RealT>> up = ROL::makePtr<ROL::StdVector<RealT>>(usv);
  ROL::StdVector<RealT> result(resultp);
  ROL::StdVector<RealT> target(targetp);

  // Problem 1
  (*vsv)[0] =  4.0;
  (*xsv)[0] =  1.9;
  (*gsv)[0] =  0.5;
  (*lsv)[0] =  0.0;
  (*usv)[0] =  2.0;

  // Problem 2
  (*vsv)[1] = -1.0;
  (*xsv)[1] = 10.0;
  (*gsv)[1] =  0.0002;
  (*lsv)[1] = ROL::ROL_NINF<RealT>();
  (*usv)[1] = ROL::ROL_INF<RealT>();

  // Problem 3
  (*vsv)[2] = -0.0002;
  (*xsv)[2] =  1.0;
  (*gsv)[2] =  0.5;
  (*lsv)[2] =  0.0;
  (*usv)[2] = ROL::ROL_INF<RealT>();

  ROL::Bounds<RealT>             boundsBC( lp,  up );
  ROL::StdBoundConstraint<RealT> stdvecBC(*lsv,*usv);

  // Expected results when applying the scaling function to v.
  (*targetp)[0] = (*vsv)[0]*(*gsv)[0];
  (*targetp)[1] = (*vsv)[1]*1.0;
  (*targetp)[2] = (*vsv)[2]*1.0;

  target.applyBinary(ROL::Elementwise::DivideAndInvert<RealT>(), v);
  target.applyBinary(ROL::Elementwise::Multiply<RealT>(), v);
  boundsBC.applyInverseScalingFunction(result, v, x, g);
  ewErrors.push_back(calcError(result, target));
  stdvecBC.applyInverseScalingFunction(result, v, x, g);
  svErrors.push_back(calcError(result, target));

  // Expected results when applying the Jacobian to v.
  (*targetp)[0] = (*vsv)[0]*(*gsv)[0];
  (*targetp)[1] = 0.0;
  (*targetp)[2] = 0.0;
  boundsBC.applyScalingFunctionJacobian(result, v, x, g);
  ewErrors.push_back(calcError(result, target));
  stdvecBC.applyScalingFunctionJacobian(result, v, x, g);
  svErrors.push_back(calcError(result, target));

  *outStream << "Elementwise Test Case Errors (Infinity Norm):" << std::endl
    << "  Inverse  = " << ewErrors[1] << std::endl
    << "  Jacobian = " << ewErrors[2] << std::endl;
  *outStream << "StdBoundConstraint Test Case Errors (Infinity Norm):" << std::endl
    << "  Inverse  = " << svErrors[1] << std::endl
    << "  Jacobian = " << svErrors[2] << std::endl;
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
