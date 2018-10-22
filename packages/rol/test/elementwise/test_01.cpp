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
    \brief Test std::vector interface for elementwise functions
*/

#include "ROL_StdVector.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"


// Unary function
template<class Real>
Real Reciprocal( const Real &val ) {
  return 1.0/val;
}     

// Binary function
template<class Real>
Real Product( const Real &x, const Real &y ) {
  return x*y;
}

// Unary Functor
template<class Real>
class Threshold {
public:
  Threshold(const Real &value) : value_(value) {}
  Real operator()(const Real &x) const {
    return x > value_ ? value_ : x;
  }
private:
  Real value_;
};

// Unary Functor using inheritance
template<class Real>
class Square : public ROL::Elementwise::UnaryFunction<Real> {
public:
  Real apply( const Real &x ) const {
    return x*x;
  }
};





typedef double RealT;

int main(int argc, char *argv[]) {

  using namespace Teuchos;
 
  typedef std::vector<RealT>     vec;
  typedef ROL::StdVector<RealT>  V;

  GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  RealT errtol = ROL::ROL_THRESHOLD<RealT>();

  // *** Test body.

  try {

    int k = 5;
    int dim = k*k;
    RealT threshValue = 4.0;

    ROL::Ptr<vec> w_ptr        = ROL::makePtr<vec>(dim, 0.0);
    ROL::Ptr<vec> w2_ptr       = ROL::makePtr<vec>(dim, 0.0);
    ROL::Ptr<vec> x_ptr        = ROL::makePtr<vec>(dim, 0.0);
    ROL::Ptr<vec> x_recip_ptr  = ROL::makePtr<vec>(dim, 0.0);
    ROL::Ptr<vec> y_ptr        = ROL::makePtr<vec>(dim, 0.0);
    ROL::Ptr<vec> z_ptr        = ROL::makePtr<vec>(dim, 0.0);
    ROL::Ptr<vec> z_thresh_ptr = ROL::makePtr<vec>(dim, 0.0);

    V w(w_ptr);
    V w2(w2_ptr);
    V x(x_ptr);
    V x_recip(x_recip_ptr);
    V y(y_ptr);
    V z(z_ptr);
    V z_thresh(z_thresh_ptr);

    for(int i=0; i<dim; ++i) {
      (*w_ptr)[i] = 1.0+i;

      (*w2_ptr)[i] = std::pow(1.0+i,2);   
      (*x_recip_ptr)[i]  = 1.0/(1.0+i);
      (*z_thresh_ptr)[i] = std::min(1.0+i,threshValue);
    }

    x.set(w);
    y.set(w);
    z.set(w);

    // Test Unary Function with inheritance
    Square<RealT> square;
    w.applyUnary(square);

    w.axpy(-1.0,w2);

    errorFlag += ( w.norm() > errtol ) ? 1 : 0;

    // Test Unary Function with wrapper
    ROL::Elementwise::applyUnaryInPlace(x,Reciprocal<RealT>);
 
    x.axpy(-1.0,x_recip);

    errorFlag += ( x.norm() > errtol ) ? 1 : 0;

    // Test Binary Function with wrapper
    ROL::Elementwise::applyBinaryInPlace(y,x_recip,Product<RealT>);    
    
    errorFlag += ( std::abs(y.norm()-static_cast<RealT>(k)) > errtol ) ? 1 : 0; 

    // Test Unary Functor with wrapper
    Threshold<RealT> threshold(threshValue);
    ROL::Elementwise::applyUnaryInPlace(z,threshold);

    z.axpy(-1.0,z_thresh);

    errorFlag += ( z.norm() > errtol ) ? 1 : 0;
   
    // Test reduce 
    ROL::Elementwise::ReductionMax<RealT> maximum;

    errorFlag += std::abs(w2.reduce(maximum)-pow(dim,2)) > errtol ? 1 : 0;

  }
  catch (std::logic_error err) {
    *outStream << err.what() << std::endl;
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;

}

