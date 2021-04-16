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

/*! \file  test_15.cpp
    \brief Validate polyhedral projection infrastructure.

*/

#include "ROL_NullSpaceOperator.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_ScaledStdVector.hpp"
#include "ROL_StdConstraint.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

template<typename Real>
class con2d : public ROL::StdConstraint<Real> {
public:
  void value(std::vector<Real> &c, const std::vector<Real> &x, Real &tol) {
    c[0] = x[0]+x[1];
  }
  void applyJacobian(std::vector<Real> &jv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol) {
    jv[0] = v[0]+v[1];
  }
  void applyAdjointJacobian(std::vector<Real> &ajv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol) {
    ajv[0] = v[0];
    ajv[1] = v[0];
  }
};

typedef double RealT;

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

  try {
    RealT tol = std::sqrt(ROL::ROL_EPSILON<RealT>());
    RealT err(0);
    ROL::Ptr<con2d<RealT>> con = ROL::makePtr<con2d<RealT>>();

    ROL::Ptr<std::vector<RealT>> yptr = ROL::makePtr<std::vector<RealT>>(2);
    (*yptr)[0] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
    (*yptr)[1] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
    ROL::StdVector<RealT> y(yptr);

    ROL::StdVector<RealT> r(1);

    ROL::Ptr<std::vector<RealT>> xptr = ROL::makePtr<std::vector<RealT>>(2);
    (*xptr)[0] = (*yptr)[0];
    (*xptr)[1] = (*yptr)[1];
    ROL::StdVector<RealT> x(xptr);

    ROL::Ptr<std::vector<RealT>> Pxptr = ROL::makePtr<std::vector<RealT>>(2,0.0);
    ROL::StdVector<RealT> Px(Pxptr);

    ROL::NullSpaceOperator<RealT> ns0(con,x,r);
    ns0.apply(Px,x,tol);
 
    ROL::Ptr<std::vector<RealT>> x0ptr = ROL::makePtr<std::vector<RealT>>(2);
    (*x0ptr)[0] = ((*yptr)[0]-(*yptr)[1])/static_cast<RealT>(2);
    (*x0ptr)[1] = -(*x0ptr)[0];
    ROL::StdVector<RealT> x0(x0ptr);

    ROL::StdVector<RealT> e0(2);

    *outStream << std::setprecision(6) << std::scientific << std::endl;
    *outStream << "   x[0] = " <<  (*xptr)[0] << "   x[1] = " <<  (*xptr)[1] << std::endl;
    *outStream << "  Px[0] = " << (*Pxptr)[0] << "  Px[1] = " << (*Pxptr)[1] << std::endl;
    *outStream << "  x*[0] = " << (*x0ptr)[0] << "  x*[1] = " << (*x0ptr)[1] << std::endl;

    e0.set(x0); e0.axpy(static_cast<RealT>(-1),Px);
    err = e0.norm();
    *outStream << "  Error in Euclidean Projection: " << err << std::endl;

    e0.set(x); e0.axpy(static_cast<RealT>(-1),x0);
    *outStream << "  ||x*-x||^2 = " << e0.norm() << std::endl;

    e0.set(x); e0.axpy(static_cast<RealT>(-1),Px);
    *outStream << "  ||Px-x||^2 = " << e0.norm() << std::endl << std::endl;

    errorFlag += (err > tol);

    ROL::Ptr<std::vector<RealT>> dptr = ROL::makePtr<std::vector<RealT>>(2);
    (*dptr)[0] = static_cast<RealT>(1)+static_cast<RealT>(2)*static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
    (*dptr)[1] = static_cast<RealT>(1)+static_cast<RealT>(5)*static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
 
    ROL::Ptr<std::vector<RealT>> x1ptr = ROL::makePtr<std::vector<RealT>>(2);
    (*x1ptr)[0] = ((*dptr)[0]*(*yptr)[0]-(*dptr)[1]*(*yptr)[1])/((*dptr)[0]+(*dptr)[1]);
    (*x1ptr)[1] = -(*x1ptr)[0];
    ROL::PrimalScaledStdVector<RealT> x1(x1ptr,dptr);

    ROL::Ptr<std::vector<RealT>> zptr = ROL::makePtr<std::vector<RealT>>(2);
    (*zptr)[0] = (*yptr)[0];
    (*zptr)[1] = (*yptr)[1];
    ROL::PrimalScaledStdVector<RealT> z(zptr,dptr);

    ROL::Ptr<std::vector<RealT>> Pzptr = ROL::makePtr<std::vector<RealT>>(2,0.0);
    ROL::PrimalScaledStdVector<RealT> Pz(Pzptr,dptr);

    ROL::NullSpaceOperator<RealT> ns1(con,z,r);
    ns1.apply(Pz,z,tol);

    ROL::Ptr<std::vector<RealT>> e1ptr = ROL::makePtr<std::vector<RealT>>(2);
    ROL::PrimalScaledStdVector<RealT> e1(e1ptr,dptr);

    *outStream << std::endl;
    *outStream << "   x[0] = " <<  (*zptr)[0] << "   x[1] = " <<  (*zptr)[1] << std::endl;
    *outStream << "  Px[0] = " << (*Pzptr)[0] << "  Px[1] = " << (*Pzptr)[1] << std::endl;
    *outStream << "  x*[0] = " << (*x1ptr)[0] << "  x*[1] = " << (*x1ptr)[1] << std::endl;

    e1.set(x1); e1.axpy(static_cast<RealT>(-1),Pz);
    err = e1.norm();
    *outStream << "  Error in Euclidean Projection: " << err << std::endl;

    e1.set(z); e1.axpy(static_cast<RealT>(-1),x1);
    *outStream << "  ||x*-x||^2 = " << e1.norm() << std::endl;

    e1.set(z); e1.axpy(static_cast<RealT>(-1),Pz);
    *outStream << "  ||Px-x||^2 = " << e1.norm() << std::endl << std::endl;

    errorFlag += (err > tol);
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

