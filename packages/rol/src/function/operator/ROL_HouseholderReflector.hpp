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

#ifndef ROL_HOUSEHOLDERREFLECTOR_H
#define ROL_HOUSEHOLDERREFLECTOR_H

#include "ROL_LinearOperator.hpp"

/** @ingroup func_group
    \class ROL::HouseholderReflector
    \brief Provides the interface to create a Householder reflector 
           operator, that when applied to a vector x, produces a 
           vector parallel to y

    Forms the linear operator H such that 

    \f[
    Hv = (I-2u\otimes u) x = \frac{||x||}{||y||} y
    \f]

    where 
    \f[
    \hat u = x + \sign(\langle x,u\rangle)\frac{||x||}{||y||}y,\quad u = \frac{\hat u}{||\hat u||}
    \f]

    If y is not specified, it is taken to be the first canonical vector

    ---
*/


namespace ROL {

template <class Real>
class HouseholderReflector : public LinearOperator<Real> {

  typedef Vector<Real>  V;
  
private:

  const ROL::Ptr<const V> x_;
  const ROL::Ptr<const V> y_;

  ROL::Ptr<V> u_;

public:
  
  HouseholderReflector( const ROL::Ptr<const Vector<Real> > &x, 
                        const ROL::Ptr<const Vector<Real> > &y) : x_(x), y_(y), u_(x->clone()) {}
  

  HouseholderReflector( const ROL::Ptr<const Vector<Real> > &x,
                        const ROL::Ptr<const Vector<Real> > &y,
                        ROL::Ptr<Vector<Real> > &scratch ) : x_(x), y_(y), u_(scratch) {}

  HouseholderReflector( const ROL::Ptr<const Vector<Real> > &x ) : x_(x), y_(x->basis(0)), u_(x->clone()) {}
  
  HouseholderReflector( const ROL::Ptr<const Vector<Real> > &x,
                        ROL::Ptr<Vector<Real> > &scratch ) : x_(x), y_(x->basis(0)), u_(scratch) {}
  
                      

  void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {

    Real xdoty = x_->dot(*y_);
    Real xnorm = x_->norm();
    Real ynorm = y_->norm();
    Real sgn = xdoty/std::abs(xdoty);

    Real alpha = sgn*xnorm/ynorm;

    u_->set(*x_);
    u_->axpy(alpha,*y_);

    Real beta  = -2.0*u_->dot(v)/u_->dot(*u_);
   
    Hv.set(v);
    Hv.axpy(beta,*u_);
  }

  void applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    apply(Hv,v,tol); 
  }

}; // class HouseholderReflector

} // namespace ROL

#endif // ROL_HOUSEHOLDERREFLECTOR_H
