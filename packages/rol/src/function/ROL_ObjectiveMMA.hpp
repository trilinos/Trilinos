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

#ifndef ROL_OBJECTIVEMMA_H
#define ROL_OBJECTIVEMMA_H

#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"

/** @ingroup func_group
    \class ROL::ObjectiveMMA
    \brief Provides the interface to to Method of Moving Asymptotes 
           Objective function

    ---
*/


namespace ROL {

template<class Real>
class ObjectiveMMA : public Objective<Real> {

  template <typename T> using ROL::Ptr = ROL::Ptr<T>;

  typedef Objective<Real>       OBJ;
  typedef BoundConstraint<Real> BND;

private:

  const ROL::Ptr<OBJ> obj_;
  const ROL::Ptr<BND> bnd_;
  

  ROL::Ptr<V> l_; // Lower bound
  ROL::Ptr<V> u_; // Upper bound
  
  ROL::Ptr<V> p_; // First MMA numerator
  ROL::Ptr<V> q_; // Second MMA numerator

  ROL::Ptr<V> d_; // Scratch vector

  Real fval_; // Original objective value

  Real tol_;

public:

  ObjectiveMMA( const ROL::Ptr<Objective<Real> > &obj,
                const ROL::Ptr<BoundConstraint<Real> > &bnd,
                const Vector<Real> &x,
                Real tol=std::sqrt(ROL_EPSILON<Real>()) ) : 
    obj_(obj), bnd_(bnd), tol_(tol) {

    l_   = bnd_->getLowerBound();
    u_   = bnd_->getUpperBound();

    p_   = x.clone();
    q_   = x.clone();
    d_   = x.clone();
 
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
  
    Elementwise::ThresholdUpper<Real> positive(0.0);
    Elementwise::Power<Real>          square(2.0);
    Elementwise::Multiply<Real>       mult;

    obj_->update(x,flag,iter);
   
    fval_ = obj_->value(x,tol);
    obj_->gradient(*p_,x,tol);
    q_->set(*p_);

    p_->applyUnary(positive);
    q_->applyUnary(negative);

    d_->set(x);
    d_->axpy(-1.0,*l_);
    d_->applyUnary(square);
    p_->applyBinary(mult,*d_);

    d_->set(*u_);
    d_->axpy(-1.0,x);
    d_->applyUnary(square);
    q_->applyBinary(mult,*d_);

  }

  /*
     \f[ F(x) \approx F(x^0) + \sum\limit_{i=1}^n \left( \frac{p_i}{U_i-x_i} + \frac{q_i}{x_i-L_i}\right) \f] 
  */
  Real value( const Vector<Real> &x, Real &tol ) {
  
    Elementwise::ReductionSum<Real>    sum;
    Elementwise::DivideAndInvert<Real> divinv;
    Real fval = fval_;

    d_->set(*u_);
    d_->axpy(-1.0,x);
    d_->applyBinary(divinv,*p_);  

    fval += d_->reduce(sum);

    d_->set(x);
    d_->axpy(-1.0,*l_);
    d_->applyBinary(divinv,*q_);

    fval += d_->reduce(sum);   

    return fval;

  }

  /*
     \f[ \frac{F(x)}{\partial x_j} =  \frac{p_j}{(U_j-x_j)^2} - \frac{q_j}({x_j-L_j)^2}\ \f] 
  */
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {

    Elementwise::DivideAndInvert<Real> divinv;
    Elementwise::Power<Real>           square(2.0);

    d_->set(*u_);
    d_->axpy(-1.0,x);
    d_->applyUnary(square);
    d_->applyBinary(divinv,*p_);            

    g.set(*d_);

    d_->set(x);
    d_->axpy(-1.0,*l_);
    d_->applyUnary(square);
    d_->applyBinary(divinv,*q_);

    g.plus(*d_);

  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  
    Elementwise::DivideAndInvert<Real> divinv;
    Elementwise::Multiply<Real>        mult;
    Elementwise::Power<Real>           cube(3.0);

    d_->set(*u_);
    d_->axpy(-1.0,x);
    d_->applyUnary(cube);          
    d_->applyBinary(divinv,*p_);           
    d_->scale(-2.0);

    hv.set(*d_); 

    d_->set(x);
    d_->axpy(-1.0,*l_);
    d_->applyUnary(cube);
    d_->applyBinary(divinv,*q_);
    d_->scale(2.0);
    
    hv.plus(*d_);
    hv.applyBinary(mult,v);

  }

  void invHessVec( Vector<Real> &h, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {

    Elementwise::DivideAndInvert<Real> divinv;
    Elementwise::Multiply<Real>        mult;
    Elementwise::Power<Real>           cube(3.0);

    d_->set(*u_);
    d_->axpy(-1.0,x);
    d_->applyUnary(cube);          
    d_->applyBinary(divinv,*p_);           
    d_->scale(-2.0);

    hv.set(*d_); 

    d_->set(x);
    d_->axpy(-1.0,*l_);
    d_->applyUnary(cube);
    d_->applyBinary(divinv,*q_);
    d_->scale(2.0);
    
    hv.plus(*d_);
    hv.applyBinary(divinv,v);
      
  }

}; // class ObjectiveMMA 

} // namespace ROL





#endif // ROL_OBJECTIVEMMA_H

