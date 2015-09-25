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

#ifndef ROL_BOUND_CONSTRAINT_H
#define ROL_BOUND_CONSTRAINT_H

#include "ROL_Vector.hpp"
#include "ROL_Types.hpp"
#include <iostream>

/** @ingroup func_group
    \class ROL::BoundConstraint
    \brief Provides the interface to apply upper and lower bound constraints.

    ROL's bound constraint class is to designed to handle point wise bound 
    constraints on optimization variables.  That is, let \f$\mathcal{X}\f$ 
    be a Banach space of functions from \f$\Xi\f$ into \f$\mathbb{R}\f$ 
    (for example, \f$\Xi\subset\mathbb{R}^d\f$ for some positive integer \f$d\f$
    and \f$\mathcal{X}=L^2(\Xi)\f$ or \f$\Xi = \{1,\ldots,n\}\f$ and 
    \f$\mathcal{X}=\mathbb{R}^n\f$).  For any \f$x\in\mathcal{X}\f$, we consider 
    bounds of the form 
    \f[
        a(\xi) \le x(\xi) \le b(\xi) \quad \text{for almost every }\xi\in\Xi.
    \f] 
    Here, \f$a(\xi)\le b(\xi)\f$ for almost every \f$\xi\in\Xi\f$ and \f$a,b\in \mathcal{X}\f$.
*/


namespace ROL {

template <class Real>
class BoundConstraint {
private:
  int dim_;
  Teuchos::RCP<Vector<Real> > x_lo_;
  Teuchos::RCP<Vector<Real> > x_up_;
  Real scale_;  
  bool activated_; ///< Flag that determines whether or not the constraints are being used.
  Real min_diff_;

  Elementwise::ReductionMin<Real> minimum_;

  // Add a constant to every element in a vector
  class Shift : public Elementwise::UnaryFunction<Real> {
    public:
      Shift(Real offset) : offset_(offset) {}
      Real apply( const Real &x ) const { return offset_ + x; }
    private:
      Real offset_;   
  };
  
  class LogicalOr : public Elementwise::BinaryFunction<Real> {
    public:
      Real apply( const Real &x, const Real &y ) const {
        return (x==0.0 && y==0.0) ? 0.0 : 1.0;
      }
  };

  class Product : public Elementwise::BinaryFunction<Real> {
    public:
      Real apply( const Real &x, const Real &y ) const { return x*y; }
  } product;


public:

  virtual ~BoundConstraint() {}

  BoundConstraint(void) : activated_(true) {}  

  /** \brief Default constructor.

      The default constructor automatically turns the constraints on.
  */
  BoundConstraint(const Teuchos::RCP<Vector<Real> > &x_lo,
                  const Teuchos::RCP<Vector<Real> > &x_up,
                  Real scale = 1.0) : 
                    x_lo_(x_lo), x_up_(x_up),scale_(scale),activated_(true) {

    // Compute difference between upper and lower bounds
    Teuchos::RCP<Vector<Real> > diff = x_up_->clone();
    diff->set(*x_up_);
    diff->axpy(-1.0,*x_lo_);  

    // Compute minimum difference 
    min_diff_ = diff->reduce(minimum_);
    min_diff_ *= 0.5; 
    
  }

  /** \brief Update bounds. 

      The update function allows the user to update the bounds at each new iterations. 
          @param[in]      x      is the optimization variable.
          @param[in]      flag   is set to true if control is changed.
          @param[in]      iter   is the outer algorithm iterations count.
  */
  virtual void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {}

  /** \brief Project optimization variables onto the bounds.

      This function implements the projection of \f$x\f$ onto the bounds, i.e., 
      \f[
         (P_{[a,b]}(x))(\xi) = \min\{b(\xi),\max\{a(\xi),x(\xi)\}\} \quad \text{for almost every }\xi\in\Xi. 
      \f]
       @param[in,out]      x is the optimization variable.
  */
  virtual void project( Vector<Real> &x ) {

    struct Lesser : public Elementwise::BinaryFunction<Real> {
      Real apply(const Real &x, const Real &y) const { return x<y ? x : y; }   
    } lesser;
   
    struct Greater : public Elementwise::BinaryFunction<Real> {
      Real apply(const Real &x, const Real &y) const { return x>y ? x : y; }   
    } greater;

    x.applyBinary(lesser, *x_up_); // Set x to the elementwise minimum of x and x_up_
    x.applyBinary(greater,*x_lo_); // Set x to the elementwise minimum of x and x_lo_

  }

  /** \brief Set variables to zero if they correspond to the upper \f$\epsilon\f$-active set.
  
      This function sets \f$v(\xi)=0\f$ if \f$\xi\in\mathcal{A}^+_\epsilon(x)\f$.  Here, 
      the upper \f$\epsilon\f$-active set is defined as 
      \f[
         \mathcal{A}^+_\epsilon(x) = \{\,\xi\in\Xi\,:\,x(\xi) = b(\xi)-\epsilon\,\}.
      \f]
      @param[out]      v   is the variable to be pruned.
      @param[in]       x   is the current optimization variable.
      @param[in]       eps is the active-set tolerance \f$\epsilon\f$.
  */
  virtual void pruneUpperActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0.0 ) {

    Real epsn = std::min(scale_*eps,this->min_diff_);          

    // Find the indices where x - u + epsn >= 0

    Shift shift(epsn);    

    Teuchos::RCP<Vector<Real> > mask = x.clone();
    mask->set(x);
    mask->axpy(-1.0,*x_up_);
    mask->applyUnary(shift); 
    
    struct Condition : public Elementwise::UnaryFunction<Real> {
      Real apply(const Real &x) const {
        return x>=0 ? 0.0 : 1.0;
      }
    } condition;   

    mask->applyUnary(condition); // mask = 0 anywhere x - u + epsn >= 0

    v.applyBinary(product, *mask); // zeros elements of v where mask=0 
  }

  /** \brief Set variables to zero if they correspond to the upper \f$\epsilon\f$-binding set.
  
      This function sets \f$v(\xi)=0\f$ if \f$\xi\in\mathcal{B}^+_\epsilon(x)\f$.  Here, 
      the upper \f$\epsilon\f$-binding set is defined as 
      \f[
         \mathcal{B}^+_\epsilon(x) = \{\,\xi\in\Xi\,:\,x(\xi) = b(\xi)-\epsilon,\; 
                g(\xi) < 0 \,\}.
      \f]
      @param[out]      v   is the variable to be pruned.
      @param[in]       x   is the current optimization variable.
      @param[in]       g   is the negative search direction.
      @param[in]       eps is the active-set tolerance \f$\epsilon\f$.
  */
  virtual void pruneUpperActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0.0 ) {

    Real epsn = std::min(scale_*eps,this->min_diff_);          

    // Find the indices where x - u + epsn >= 0

    Shift shift(epsn);    

    Teuchos::RCP<Vector<Real> > mask1 = x.clone();
    mask1->set(x);
    mask1->axpy(-1.0,*x_up_);
    mask1->applyUnary(shift); 
    
    struct Condition1 : public Elementwise::UnaryFunction<Real> {
      Real apply(const Real &x) const {
        return x>=0 ? 0.0 : 1.0;
      }
    } condition1;   

    mask1->applyUnary(condition1);

    struct Condition2 : public Elementwise::UnaryFunction<Real> {  
      Real apply(const Real &x) const {
        return x < 0.0 ? 0.0 : 1.0;
      }
    } condition2;

    Teuchos::RCP<Vector<Real> > mask2 = g.clone();
    mask2->set(g);    

    mask2->applyUnary(condition2); 

    LogicalOr logicalOr;

    mask1->applyBinary(logicalOr,*mask2); 

    v.applyBinary(product,*mask1);
 
  }
 
  /** \brief Set variables to zero if they correspond to the lower \f$\epsilon\f$-active set.
  
      This function sets \f$v(\xi)=0\f$ if \f$\xi\in\mathcal{A}^-_\epsilon(x)\f$.  Here, 
      the lower \f$\epsilon\f$-active set is defined as 
      \f[
         \mathcal{A}^-_\epsilon(x) = \{\,\xi\in\Xi\,:\,x(\xi) = a(\xi)+\epsilon\,\}.
      \f]
      @param[out]      v   is the variable to be pruned.
      @param[in]       x   is the current optimization variable.
      @param[in]       eps is the active-set tolerance \f$\epsilon\f$.
  */
  virtual void pruneLowerActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0.0 ) {

    Real epsn = std::min(scale_*eps,this->min_diff_);          

    // Find the indices where -x + l + epsn >= 0

    Shift shift(epsn);    

    Teuchos::RCP<Vector<Real> > mask = x_lo_->clone();
    mask->set(*x_lo_);
    mask->axpy(-1.0,x);
    mask->applyUnary(shift); 
    
    struct Condition: public Elementwise::UnaryFunction<Real> {
      Real apply(const Real &x) const {
        return x>=0 ? 0.0 : 1.0;
      }
    } condition;   

    mask->applyUnary(condition); // mask = 0 anywhere -x + l + epsn >= 0

    v.applyBinary(product, *mask); // zeros elements of v where mask=0 

  }

  /** \brief Set variables to zero if they correspond to the lower \f$\epsilon\f$-binding set.
  
      This function sets \f$v(\xi)=0\f$ if \f$\xi\in\mathcal{B}^-_\epsilon(x)\f$.  Here, 
      the lower \f$\epsilon\f$-binding set is defined as 
      \f[
         \mathcal{B}^-_\epsilon(x) = \{\,\xi\in\Xi\,:\,x(\xi) = a(\xi)+\epsilon,\; 
                g(\xi) > 0 \,\}.
      \f]
      @param[out]      v   is the variable to be pruned.
      @param[in]       x   is the current optimization variable.
      @param[in]       g   is the negative search direction.
      @param[in]       eps is the active-set tolerance \f$\epsilon\f$.
  */
  virtual void pruneLowerActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0.0 ) {

    Real epsn = std::min(scale_*eps,this->min_diff_);          

    // Find the indices where -x + l + epsn >= 0

    Shift shift(epsn);    

    Teuchos::RCP<Vector<Real> > mask1 = x_lo_->clone();
    mask1->set(*x_lo_);
    mask1->axpy(-1.0,x);
    mask1->applyUnary(shift); 
    
    struct Condition1: public Elementwise::UnaryFunction<Real> {
      Real apply(const Real &x) const {
        return x>=0 ? 0.0 : 1.0;
      }
    } condition1;   

    mask1->applyUnary(condition1); // mask = 0 anywhere -x + l + epsn >= 0

    struct Condition2 : public Elementwise::UnaryFunction<Real> {  
      Real apply(const Real &x) const {
        return x > 0.0 ? 0.0 : 1.0;
      }
    } condition2;
    
    Teuchos::RCP<Vector<Real> > mask2 = g.clone();
    mask2->set(g); 

    mask2->applyUnary(condition2); 

    LogicalOr logicalOr;

    mask1->applyBinary(logicalOr,*mask2); 

    v.applyBinary(product,*mask1);

  }
 
  /** \brief Set the input vector to the upper bound.

      This function sets the input vector \f$u\f$ to the upper bound \f$b\f$.
      @param[out]    u   is the vector to be set to the upper bound.
  */ 
  virtual void setVectorToUpperBound( Vector<Real> &u ) {
    u.set(*x_up_);
  }

  /** \brief Set the input vector to the lower bound.

      This function sets the input vector \f$l\f$ to the lower bound \f$a\f$.
      @param[out]    l   is the vector to be set to the lower bound.
  */ 
  virtual void setVectorToLowerBound( Vector<Real> &l ) {
    l.set(*x_lo_);
  }

  /** \brief Set variables to zero if they correspond to the \f$\epsilon\f$-active set.
  
      This function sets \f$v(\xi)=0\f$ if \f$\xi\in\mathcal{A}_\epsilon(x)\f$.  Here, 
      the \f$\epsilon\f$-active set is defined as 
      \f[
         \mathcal{A}_\epsilon(x) = \mathcal{A}^+_\epsilon(x)\cap\mathcal{A}^-_\epsilon(x).
      \f]
      @param[out]      v   is the variable to be pruned.
      @param[in]       x   is the current optimization variable.
      @param[in]       eps is the active-set tolerance \f$\epsilon\f$.
  */
  virtual void pruneActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0.0 ) {
    this->pruneUpperActive(v,x,eps);
    this->pruneLowerActive(v,x,eps);
  }

  /** \brief Set variables to zero if they correspond to the \f$\epsilon\f$-binding set.
  
      This function sets \f$v(\xi)=0\f$ if \f$\xi\in\mathcal{B}_\epsilon(x)\f$.  Here, 
      the \f$\epsilon\f$-binding set is defined as 
      \f[
         \mathcal{B}^+_\epsilon(x) = \mathcal{B}^+_\epsilon(x)\cap\mathcal{B}^-_\epsilon(x).
      \f]
      @param[out]      v   is the variable to be pruned.
      @param[in]       x   is the current optimization variable.
      @param[in]       g   is the negative search direction.
      @param[in]       eps is the active-set tolerance \f$\epsilon\f$.
  */
  virtual void pruneActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0.0 ) {
    this->pruneUpperActive(v,g,x,eps);
    this->pruneLowerActive(v,g,x,eps);
  }

  /** \brief Check if the vector, v, is feasible.

      This function returns true if \f$v = P_{[a,b]}(v)\f$.
      @param[in]    v   is the vector to be checked.
  */
  virtual bool isFeasible( const Vector<Real> &v ) { 
    if ( this->activated_ ) {
      Teuchos::RCP<Vector<Real> > s = v.clone();
      s->set(*x_up_);

      s->axpy(-1.0,v); // u-v

      Real uminusv = s->reduce(minimum_);

      s->set(v);
      s->axpy(-1.0,*x_lo_);

      Real vminusl = s->reduce(minimum_);

      if( (uminusv<0) || (vminusl<0) ) {
        return false;
      } 
      else {
        return true;
      }

    }
    else {
      return true; 
    }
  }

  /** \brief Turn on bounds.
   
      This function turns the bounds on. 
  */
  void activate(void)    { this->activated_ = true;  }

  /** \brief Turn off bounds.

      This function turns the bounds off.
  */
  void deactivate(void)  { this->activated_ = false; }

  /** \brief Check if bounds are on.

      This function returns true if the bounds are turned on.
  */
  bool isActivated(void) { return this->activated_;  }

  /** \brief Set variables to zero if they correspond to the \f$\epsilon\f$-inactive set.
  
      This function sets \f$v(\xi)=0\f$ if \f$\xi\in\Xi\setminus\mathcal{A}_\epsilon(x)\f$.  Here, 
      @param[out]      v   is the variable to be pruned.
      @param[in]       x   is the current optimization variable.
      @param[in]       eps is the active-set tolerance \f$\epsilon\f$.
  */
  void pruneInactive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0.0 ) { 
    Teuchos::RCP<Vector<Real> > tmp = v.clone(); 
    tmp->set(v);
    this->pruneActive(*tmp,x,eps);
    v.axpy(-1.0,*tmp);
  }

  /** \brief Set variables to zero if they correspond to the \f$\epsilon\f$-nonbinding set.
  
      This function sets \f$v(\xi)=0\f$ if \f$\xi\in\Xi\setminus\mathcal{B}_\epsilon(x)\f$.  
      @param[out]      v   is the variable to be pruned.
      @param[in]       x   is the current optimization variable.
      @param[in]       g   is the negative search direction.
      @param[in]       eps is the active-set tolerance \f$\epsilon\f$.
  */
  void pruneInactive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0.0 ) { 
    Teuchos::RCP<Vector<Real> > tmp = v.clone(); 
    tmp->set(v);
    this->pruneActive(*tmp,g,x,eps);
    v.axpy(-1.0,*tmp);
  }
 
  /** \brief Compute projected gradient.

      This function projects the gradient \f$g\f$ onto the tangent cone.
           @param[in,out]    g  is the gradient of the objective function at x.
           @param[in]        x  is the optimization variable
  */
  void computeProjectedGradient( Vector<Real> &g, const Vector<Real> &x ) {
    Teuchos::RCP<Vector<Real> > tmp = g.clone();
    tmp->set(g);
    this->pruneActive(g,*tmp,x);
  }
 
  /** \brief Compute projected step.

      This function computes the projected step \f$P_{[a,b]}(x+v) - x\f$.
      @param[in,out]         v  is the step variable.
      @param[in]             x is the optimization variable.
  */
  void computeProjectedStep( Vector<Real> &v, const Vector<Real> &x ) { 
    v.plus(x);
    this->project(v);
    v.axpy(-1.0,x);
  }

}; // class BoundConstraint

} // namespace ROL

#endif
