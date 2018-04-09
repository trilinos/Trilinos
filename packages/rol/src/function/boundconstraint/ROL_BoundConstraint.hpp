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
  bool Lactivated_; ///< Flag that determines whether or not the lower bounds are being used.
  bool Uactivated_; ///< Flag that determines whether or not the upper bounds are being used.
  Ptr<Vector<Real> > lower_;
  Ptr<Vector<Real> > upper_;
  bool hasSetScalar_;

public:

  virtual ~BoundConstraint() {}

  BoundConstraint(void) : Lactivated_(true), Uactivated_(true), hasSetScalar_(false) {}

  BoundConstraint(const Vector<Real> &x) : Lactivated_(false), Uactivated_(false), hasSetScalar_(false) {
    try {
      lower_ = x.clone(); lower_->setScalar(ROL_NINF<Real>());
      upper_ = x.clone(); upper_->setScalar(ROL_INF<Real>());
      hasSetScalar_ = true;
    }
    catch(std::exception &e) {
      hasSetScalar_ = false;
      // Do nothing.  If someone calls getLowerBound or getUpperBound,
      // an exception will be thrown.
    }
  }


  // REQUIRED FUNCTIONS (VIRTUAL)

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
    if (isActivated()) {
      throw Exception::NotImplemented(">>> ROL::BoundConstraint::project: Not Implemented!");
    }
  }

  /** \brief Project optimization variables into the interior of the feasible set.

      This function implements the projection of \f$x\f$ into the interior of the
      feasible set, i.e.,
      \f[
         (\bar{P}_{[a,b]}(x))(\xi) \in (a(\xi),b(\xi))
           \quad \text{for almost every }\xi\in\Xi.
      \f]
       @param[in,out]      x is the optimization variable.
  */
  virtual void projectInterior( Vector<Real> &x ) {
    if (isActivated()) {
      throw Exception::NotImplemented(">>> ROL::BoundConstraint::projectInterior: Not Implemented!");
    }
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
  virtual void pruneUpperActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0 ) {
    if (isUpperActivated()) {
      throw Exception::NotImplemented(">>> ROL::BoundConstraint::pruneUpperActive: Not Implemented!");
    }
  }

  /** \brief Set variables to zero if they correspond to the upper \f$\epsilon\f$-binding set.

      This function sets \f$v(\xi)=0\f$ if \f$\xi\in\mathcal{B}^+_\epsilon(x)\f$.  Here,
      the upper \f$\epsilon\f$-binding set is defined as
      \f[
         \mathcal{B}^+_\epsilon(x) = \{\,\xi\in\Xi\,:\,x(\xi) = b(\xi)-\epsilon,\;
                g(\xi) < 0 \,\}.
      \f]
      @param[out]      v   is the variable to be pruned.
      @param[in]       g   is the negative search direction.
      @param[in]       x   is the current optimization variable.
      @param[in]       eps is the active-set tolerance \f$\epsilon\f$.
  */
  virtual void pruneUpperActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0 ) {
    if (isUpperActivated()) {
      throw Exception::NotImplemented(">>> ROL::BoundConstraint::pruneUpperActive: Not Implemented!");
    }
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
  virtual void pruneLowerActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0 ) {
    if (isLowerActivated()) {
      throw Exception::NotImplemented(">>> ROL::BoundConstraint::pruneLowerActive: Not Implemented!");
    }
  }

  /** \brief Set variables to zero if they correspond to the \f$\epsilon\f$-binding set.

      This function sets \f$v(\xi)=0\f$ if \f$\xi\in\mathcal{B}^-_\epsilon(x)\f$.  Here,
      the lower \f$\epsilon\f$-binding set is defined as
      \f[
         \mathcal{B}^-_\epsilon(x) = \{\,\xi\in\Xi\,:\,x(\xi) = a(\xi)+\epsilon,\;
                g(\xi) > 0 \,\}.
      \f]
      @param[out]      v   is the variable to be pruned.
      @param[in]       g   is the negative search direction.
      @param[in]       x   is the current optimization variable.
      @param[in]       eps is the active-set tolerance \f$\epsilon\f$.
  */
  virtual void pruneLowerActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0 ) {
    if (isLowerActivated()) {
      throw Exception::NotImplemented(">>> ROL::BoundConstraint::pruneLowerActive: Not Implemented!");
    }
  }


  // QUERY FUNCTIONS (VIRTUAL AND NONVIRTUAL)

  /** \brief Return the ref count pointer to the lower bound vector */
  virtual const ROL::Ptr<const Vector<Real> > getLowerBound( void ) const {
    if (hasSetScalar_) {
      return lower_;
    }
    throw Exception::NotImplemented(">>> ROL::BoundConstraint::getLowerBound: Not implemented!");
  }

  /** \brief Return the ref count pointer to the upper bound vector */
  virtual const ROL::Ptr<const Vector<Real> > getUpperBound( void ) const {
    if (hasSetScalar_) {
      return upper_;
    }
    throw Exception::NotImplemented(">>> ROL::BoundConstraint::getUpperBound: Not implemented!");
  }

  /** \brief Check if the vector, v, is feasible.

      This function returns true if \f$v = P_{[a,b]}(v)\f$.
      @param[in]    v   is the vector to be checked.
  */
  virtual bool isFeasible( const Vector<Real> &v ) { 
    if (isActivated()) {
      throw Exception::NotImplemented(">>> ROL::BoundConstraint::isFeasible: Not implemented!");
    }
    return true;
  }

  /** \brief Turn on lower bound.

      This function turns on lower bounds.
  */
  void activateLower(void) {
    Lactivated_ = true;
  }

  /** \brief Turn on upper bound.

      This function turns on upper bounds.
  */
  void activateUpper(void) {
    Uactivated_ = true;
  }

  /** \brief Turn on bounds.
   
      This function turns the bounds on. 
  */
  void activate(void) {
    activateLower();
    activateUpper();
  }

  /** \brief Turn off lower bound.

      This function turns the lower bounds off.
  */
  void deactivateLower(void) {
    Lactivated_ = false;
  }

  /** \brief Turn off upper bound.

      This function turns the upper bounds off.
  */
  void deactivateUpper(void) {
    Uactivated_ = false;
  }

  /** \brief Turn off bounds.

      This function turns the bounds off.
  */
  void deactivate(void) {
    deactivateLower();
    deactivateUpper();
  }

  /** \brief Check if lower bound are on.

      This function returns true if the lower bounds are turned on.
  */
  bool isLowerActivated(void) const {
    return Lactivated_;
  }

  /** \brief Check if upper bound are on.

      This function returns true if the upper bounds are turned on.
  */
  bool isUpperActivated(void) const {
    return Uactivated_;
  }

  /** \brief Check if bounds are on.

      This function returns true if the bounds are turned on.
  */
  bool isActivated(void) const {
    return (isLowerActivated() || isUpperActivated());
  }


  // HELPER FUNCTIONS (NONVIRTUAL)

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
  void pruneActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0 ) {
    if (isActivated()) {
      pruneUpperActive(v,x,eps);
      pruneLowerActive(v,x,eps);
    }
  }

  /** \brief Set variables to zero if they correspond to the \f$\epsilon\f$-binding set.
  
      This function sets \f$v(\xi)=0\f$ if \f$\xi\in\mathcal{B}_\epsilon(x)\f$.  Here, 
      the \f$\epsilon\f$-binding set is defined as 
      \f[
         \mathcal{B}^+_\epsilon(x) = \mathcal{B}^+_\epsilon(x)\cap\mathcal{B}^-_\epsilon(x).
      \f]
      @param[out]      v   is the variable to be pruned.
      @param[in]       g   is the negative search direction.
      @param[in]       x   is the current optimization variable.
      @param[in]       eps is the active-set tolerance \f$\epsilon\f$.
  */
  void pruneActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0 ) {
    if (isActivated()) {
      pruneUpperActive(v,g,x,eps);
      pruneLowerActive(v,g,x,eps);
    }
  }

  /** \brief Set variables to zero if they correspond to the \f$\epsilon\f$-inactive set.
  
      This function sets \f$v(\xi)=0\f$ if \f$\xi\in\Xi\setminus\mathcal{A}_\epsilon(x)\f$.  Here, 
      @param[out]      v   is the variable to be pruned.
      @param[in]       x   is the current optimization variable.
      @param[in]       eps is the active-set tolerance \f$\epsilon\f$.
  */
  void pruneLowerInactive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0 ) {
    if (isLowerActivated()) {
      const Real one(1);
      ROL::Ptr<Vector<Real> > tmp = v.clone(); 
      tmp->set(v);
      pruneLowerActive(*tmp,x,eps);
      v.axpy(-one,*tmp);
    }
  }

  /** \brief Set variables to zero if they correspond to the \f$\epsilon\f$-inactive set.
  
      This function sets \f$v(\xi)=0\f$ if \f$\xi\in\Xi\setminus\mathcal{A}_\epsilon(x)\f$.  Here, 
      @param[out]      v   is the variable to be pruned.
      @param[in]       x   is the current optimization variable.
      @param[in]       eps is the active-set tolerance \f$\epsilon\f$.
  */
  void pruneUpperInactive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0 ) { 
    if (isUpperActivated()) {
      const Real one(1);
      ROL::Ptr<Vector<Real> > tmp = v.clone(); 
      tmp->set(v);
      pruneUpperActive(*tmp,x,eps);
      v.axpy(-one,*tmp);
    }
  }

  /** \brief Set variables to zero if they correspond to the \f$\epsilon\f$-nonbinding set.
  
      This function sets \f$v(\xi)=0\f$ if \f$\xi\in\Xi\setminus\mathcal{B}_\epsilon(x)\f$.  
      @param[out]      v   is the variable to be pruned.
      @param[in]       x   is the current optimization variable.
      @param[in]       g   is the negative search direction.
      @param[in]       eps is the active-set tolerance \f$\epsilon\f$.
  */
  void pruneLowerInactive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0 ) { 
    if (isLowerActivated()) {
      const Real one(1);
      ROL::Ptr<Vector<Real> > tmp = v.clone(); 
      tmp->set(v);
      pruneLowerActive(*tmp,g,x,eps);
      v.axpy(-one,*tmp);
    }
  }

  /** \brief Set variables to zero if they correspond to the \f$\epsilon\f$-nonbinding set.
  
      This function sets \f$v(\xi)=0\f$ if \f$\xi\in\Xi\setminus\mathcal{B}_\epsilon(x)\f$.  
      @param[out]      v   is the variable to be pruned.
      @param[in]       x   is the current optimization variable.
      @param[in]       g   is the negative search direction.
      @param[in]       eps is the active-set tolerance \f$\epsilon\f$.
  */
  void pruneUpperInactive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0 ) { 
    if (isUpperActivated()) {
      const Real one(1);
      ROL::Ptr<Vector<Real> > tmp = v.clone(); 
      tmp->set(v);
      pruneUpperActive(*tmp,g,x,eps);
      v.axpy(-one,*tmp);
    }
  }

  /** \brief Set variables to zero if they correspond to the \f$\epsilon\f$-inactive set.
  
      This function sets \f$v(\xi)=0\f$ if \f$\xi\in\Xi\setminus\mathcal{A}_\epsilon(x)\f$.  Here, 
      @param[out]      v   is the variable to be pruned.
      @param[in]       x   is the current optimization variable.
      @param[in]       eps is the active-set tolerance \f$\epsilon\f$.
  */
  void pruneInactive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0 ) { 
    if (isActivated()) {
      const Real one(1);
      ROL::Ptr<Vector<Real> > tmp = v.clone(); 
      tmp->set(v);
      pruneActive(*tmp,x,eps);
      v.axpy(-one,*tmp);
    }
  }

  /** \brief Set variables to zero if they correspond to the \f$\epsilon\f$-nonbinding set.
  
      This function sets \f$v(\xi)=0\f$ if \f$\xi\in\Xi\setminus\mathcal{B}_\epsilon(x)\f$.  
      @param[out]      v   is the variable to be pruned.
      @param[in]       x   is the current optimization variable.
      @param[in]       g   is the negative search direction.
      @param[in]       eps is the active-set tolerance \f$\epsilon\f$.
  */
  void pruneInactive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0 ) { 
    if (isActivated()) {
      const Real one(1);
      ROL::Ptr<Vector<Real> > tmp = v.clone(); 
      tmp->set(v);
      pruneActive(*tmp,g,x,eps);
      v.axpy(-one,*tmp);
    }
  }
 
  /** \brief Compute projected gradient.

      This function projects the gradient \f$g\f$ onto the tangent cone.
           @param[in,out]    g  is the gradient of the objective function at x.
           @param[in]        x  is the optimization variable
  */
  void computeProjectedGradient( Vector<Real> &g, const Vector<Real> &x ) {
    if (isActivated()) {
      ROL::Ptr<Vector<Real> > tmp = g.clone();
      tmp->set(g);
      pruneActive(g,*tmp,x);
    }
  }
 
  /** \brief Compute projected step.

      This function computes the projected step \f$P_{[a,b]}(x+v) - x\f$.
      @param[in,out]         v  is the step variable.
      @param[in]             x is the optimization variable.
  */
  void computeProjectedStep( Vector<Real> &v, const Vector<Real> &x ) { 
    if (isActivated()) {
      const Real one(1);
      v.plus(x);
      project(v);
      v.axpy(-one,x);
    }
  }

}; // class BoundConstraint

} // namespace ROL

#endif
