// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_03.cpp
    \brief Test ChainRuleConstraint class
*/

#include <cmath>

#include "ROL_RandomVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_ChainRuleConstraint.hpp"
#include "ROL_StdConstraint.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using RealT   = double;
using VectorT = std::vector<RealT>;


/* constraint function with (range) = 1 + dim(domain)
 */
class InnerConstraint : public ROL::StdConstraint<RealT> {
public:
  InnerConstraint( int n ) : y_(n+1), sin_(n+1), cos_(n+1), n_(n), m_(n+1) {}

  void update( const VectorT&        x, 
                     ROL::UpdateType type,
                     int             iter ) override {
    for( int k=0; k<m_; ++k ) {
      y_[k] = 0;
      if( k<n_ ) y_[k] += x[k]*x[k];
      if( k>0 )  y_[k] -= 2*x[k-1];
      sin_[k] = std::sin(y_[k]);
      cos_[k] = std::cos(y_[k]);
    }
  }    

  void value(       VectorT& c, 
              const VectorT& x, 
                    RealT& tol ) override { 
    for( int k=0; k<m_; ++k ) c[k] = std::cos(y_[k]);
  }

  void applyJacobian(       VectorT& jv, 
                      const VectorT& v, 
                      const VectorT& x, 
                            RealT&   tol ) override { 
    for( int k=0; k<m_; ++k ) {
      jv[k] = 0;
      if( k<n_ ) jv[k] -= 2*v[k]*x[k]*sin_[k];
      if( k>0 )  jv[k] += 2*v[k-1]*sin_[k];
    }
  }

  void applyAdjointJacobian(       VectorT& ajl, 
                             const VectorT& l, 
                             const VectorT& x, 
                                   RealT&   tol ) override { 
    for( int i=0; i<n_; ++i ) {
      ajl[i] = 2*(l[i+1]*sin_[i+1] - l[i]*x[i]*sin_[i]);
    }
  }

  void applyAdjointHessian(        VectorT& ahlv, 
                             const VectorT& l, 
                             const VectorT& v, 
                             const VectorT& x, 
                                   RealT&   tol ) override { 
    ahlv[0] = -4*cos_[0]*l[0]*v[0]*x[0]*x[0] - 4*cos_[1]*l[1]*v[0] + 4*cos_[1]*l[1]*v[1]*x[1] - 2*l[0]*sin_[0]*v[0];
    for(int i=0; i<n_-1; ++i) {
      ahlv[i] = 4*cos_[i]*l[i]*v[i-1]*x[i] - 4*cos_[i]*l[i]*v[i]*x[i]*x[i] - 4*cos_[i+1]*l[i+1]*v[i] + 4*cos_[i+1]*l[i+1]*v[i+1]*x[i+1] - 2*l[i]*sin_[i]*v[i];
    }
    ahlv[n_-1] = 4*cos_[n_-1]*l[n_-1]*v[n_-2]*x[n_-1] - 4*cos_[n_-1]*l[n_-1]*v[n_-1]*x[n_-1]*x[n_-1] - 4*cos_[n_]*l[n_]*v[n_-1] - 2*l[n_-1]*sin_[n_-1]*v[n_-1];
  }  
private:
  VectorT y_, sin_, cos_;
  int n_, m_;
  
};

/* constraint function with dim(range) = dim(domain) - 2
 */
class OuterConstraint : public ROL::StdConstraint<RealT> {
public:
  OuterConstraint( int n ): d_(n-2), dv_(n-2), n_(n), m_(n-2) {}

  void update( const VectorT&        x, 
                     ROL::UpdateType type,
                     int             iter ) override {
    for(int i=0; i<m_; ++i) 
      d_[i] = -x[i] + 2*x[i+1] - x[i+2];
  }

  void value(       VectorT& c, 
              const VectorT& x, 
                    RealT& tol ) override { 
    for(int i=0; i<m_; ++i) 
      c[i] = d_[i]*d_[i]; 
  }

  void applyJacobian(       VectorT& jv, 
                      const VectorT& v, 
                      const VectorT& x, 
                            RealT&   tol ) override { 
    for(int i=0; i<m_; ++i) 
      jv[i] = 2*d_[i]*(-v[i]+2*v[i+1]-v[i+2]);
  }

  void applyAdjointJacobian(       VectorT& ajl, 
                             const VectorT& l, 
                             const VectorT& x, 
                                   RealT&   tol ) override { 
    for(int i=0; i<n_; ++i) {
      ajl[i] = 0;
      for( int j=0; j<m_; ++j ) {
        if( j == i-1 ) ajl[i] += 4*d_[j]*l[j];
        else if( std::abs(i-j-1) == 1)  ajl[i] -= 2*d_[j]*l[j];
      }
    }
  }

  void applyAdjointHessian(        VectorT& ahlv, 
                             const VectorT& l, 
                             const VectorT& v, 
                             const VectorT& x, 
                                   RealT&   tol ) override { 
    for(int i=0; i<m_; ++i) 
      dv_[i] = -v[i] + 2*v[i+1] - v[i+2];
    for(int i=0; i<n_; ++i) {
      ahlv[i] = 0;
      for( int j=0; j<m_; ++j ) {
        if( j == i-1 ) ahlv[i] += 4*dv_[j]*l[j];
        else if( std::abs(i-j-1) == 1)  ahlv[i] -= 2*dv_[j]*l[j];
      }
    }
  }  

private:
  VectorT d_, dv_;
  int n_, m_;
  
};




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

  // Save the format state of the original std::cout.
  ROL::nullstream oldFormatState;
  oldFormatState.copyfmt(std::cout);

//  RealT errtol = std::sqrt(ROL::ROL_THRESHOLD<RealT>());

  int errorFlag  = 0;

  // *** Test body.

  try {

    int x_dim = 7;          // Inner constraint domain space dimension
    int ci_dim = x_dim + 1; // Inner constraint range space dimension (Outer constraint domain space dimension)
    int co_dim = x_dim - 1; // Outer constraint range space dimension

    auto inner_con = ROL::makePtr<InnerConstraint>(x_dim);
    auto outer_con = ROL::makePtr<OuterConstraint>(ci_dim);

    auto x  = ROL::makePtr<ROL::StdVector<RealT>>(x_dim);
    auto v  = ROL::makePtr<ROL::StdVector<RealT>>(x_dim);
    auto y  = ROL::makePtr<ROL::StdVector<RealT>>(ci_dim);
    auto w  = ROL::makePtr<ROL::StdVector<RealT>>(ci_dim);
    auto ci = ROL::makePtr<ROL::StdVector<RealT>>(ci_dim);
    auto co = ROL::makePtr<ROL::StdVector<RealT>>(co_dim);

    auto cr_con = ROL::makePtr<ROL::ChainRuleConstraint<RealT>>(outer_con,inner_con,*x,*ci);

    ROL::RandomizeVector(*x);
    ROL::RandomizeVector(*y);
    ROL::RandomizeVector(*v);
    ROL::RandomizeVector(*w);
    ROL::RandomizeVector(*ci);
    ROL::RandomizeVector(*co);

    *outStream << "\n\nInner Constraint Check:\n\n";

    inner_con->checkApplyJacobian(*x,*v,*ci,true,*outStream,7,4);
    inner_con->checkAdjointConsistencyJacobian(*ci,*v,*x,true,*outStream);
    inner_con->checkApplyAdjointHessian(*x,*ci,*v,*x,true,*outStream,7,4);

    *outStream << "\n\nOuter Constraint Check:\n\n";
    outer_con->checkApplyJacobian(*y,*w,*co,true,*outStream,7,4);
    outer_con->checkAdjointConsistencyJacobian(*co,*w,*y,true,*outStream);
    outer_con->checkApplyAdjointHessian(*y,*co,*w,*y,true,*outStream,7,4);

    *outStream << "\n\nChain Rule Constraint Check:\n\n";
    cr_con->checkApplyJacobian(*x,*v,*co,true,*outStream,7,4);
    cr_con->checkAdjointConsistencyJacobian(*co,*v,*x,true,*outStream);
    cr_con->checkApplyAdjointHessian(*x,*co,*v,*x,true,*outStream,7,4);

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

