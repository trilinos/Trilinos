// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_THYRA_BOUND_CONSTRAINT_H
#define ROL_THYRA_BOUND_CONSTRAINT_H

#include "TOpEleWiseBoundsHelpers.hpp"
#include "ROL_BoundConstraint.hpp"

/** @ingroup func_group
    \class ROL::Thyra_BoundConstraint
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
class Thyra_BoundConstraint: public BoundConstraint<Real> {
private:
 Teuchos::RCP<Thyra::VectorBase<Real> > thyra_x_lo_, thyra_x_up_;
  Real min_diff_;

public:

  virtual ~Thyra_BoundConstraint() {}

  Thyra_BoundConstraint(Teuchos::RCP<Thyra::VectorBase<Real> > p_min, Teuchos::RCP<Thyra::VectorBase<Real> > p_max, Real min_diff) : BoundConstraint<Real>(), thyra_x_lo_(p_min), thyra_x_up_(p_max), min_diff_(min_diff) {
  // Safety check, in case the user passes empty pointers
  if (thyra_x_lo_!=Teuchos::null && thyra_x_up_!=Teuchos::null)
    this->activate();
};

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
	  ThyraVector<Real>  & thyra_x = dynamic_cast<ThyraVector<Real>&>(x);
    Thyra::ele_wise_bound(*thyra_x_lo_, *thyra_x_up_, thyra_x.getVector().ptr() );
  }




  /** \brief Set variables to zero if they correspond to the upper \f$\epsilon\f$-active set.

      This function sets \f$v(\xi)=0\f$ if \f$\xi\in\mathcal{A}^+_\epsilon(x)\f$.  Here,
      the upper \f$\epsilon\f$-active set is defined as
      \f[
         \mathcal{A}^+_\epsilon(x) = \{\,\xi\in\Xi\,:\,x(\xi) \ge b(\xi)-\epsilon\,\}.
      \f]
      @param[out]      v   is the variable to be pruned.
      @param[in]       x   is the current optimization variable.
      @param[in]       eps is the active-set tolerance \f$\epsilon\f$.
  */
  virtual void pruneUpperActive( Vector<Real> &v, const Vector<Real> &x, Real eps = Real(0) ) {
	  const ThyraVector<Real>  & thyra_x = dynamic_cast<const ThyraVector<Real>&>(x);
	  ThyraVector<Real>  & thyra_v = dynamic_cast<ThyraVector<Real>&>(v);

	  Real epsn = std::min(eps,min_diff_);
    Thyra::ele_wise_prune_upper(*thyra_x.getVector(), *thyra_x_up_, thyra_v.getVector().ptr(), epsn );

}

  /** \brief Set variables to zero if they correspond to the upper \f$\epsilon\f$-binding set.

      This function sets \f$v(\xi)=0\f$ if \f$\xi\in\mathcal{B}^+_\epsilon(x)\f$.  Here,
      the upper \f$\epsilon\f$-binding set is defined as
      \f[
         \mathcal{B}^+_\epsilon(x) = \{\,\xi\in\Xi\,:\,x(\xi) \ge b(\xi)-\epsilon_x,\;
                g(\xi) < -\epsilon_g \,\}.
      \f]
      @param[out]      v    is the variable to be pruned.
      @param[in]       x    is the current optimization variable.
      @param[in]       g    is the negative search direction.
      @param[in]       xeps is the active-set tolerance \f$\epsilon_x\f$.
      @param[in]       geps is the binding-set tolerance \f$\epsilon_g\f$.
  */
  virtual void pruneUpperActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0) ) {
	const ThyraVector<Real>  & thyra_x = dynamic_cast<const ThyraVector<Real>&>(x);
	ThyraVector<Real>  & thyra_v = dynamic_cast<ThyraVector<Real>&>(v);
	const ThyraVector<Real>  & thyra_g = dynamic_cast<const ThyraVector<Real>&>(g);

	Real epsn = std::min(xeps,min_diff_);

  Thyra::ele_wise_prune_upper(*thyra_x.getVector(), *thyra_x_up_, *thyra_g.getVector(), thyra_v.getVector().ptr(), epsn );

  }

  /** \brief Set variables to zero if they correspond to the lower \f$\epsilon\f$-active set.

      This function sets \f$v(\xi)=0\f$ if \f$\xi\in\mathcal{A}^-_\epsilon(x)\f$.  Here,
      the lower \f$\epsilon\f$-active set is defined as
      \f[
         \mathcal{A}^-_\epsilon(x) = \{\,\xi\in\Xi\,:\,x(\xi) \le a(\xi)+\epsilon\,\}.
      \f]
      @param[out]      v   is the variable to be pruned.
      @param[in]       x   is the current optimization variable.
      @param[in]       eps is the active-set tolerance \f$\epsilon\f$.
  */
  virtual void pruneLowerActive( Vector<Real> &v, const Vector<Real> &x, Real eps = Real(0) ) {
	const ThyraVector<Real>  & thyra_x = dynamic_cast<const ThyraVector<Real>&>(x);
  ThyraVector<Real>  & thyra_v = dynamic_cast<ThyraVector<Real>&>(v);

	Real epsn = std::min(eps,min_diff_);

  Thyra::ele_wise_prune_lower(*thyra_x.getVector(), *thyra_x_lo_, thyra_v.getVector().ptr(), epsn );

  }

  /** \brief Set variables to zero if they correspond to the lower \f$\epsilon\f$-binding set.

      This function sets \f$v(\xi)=0\f$ if \f$\xi\in\mathcal{B}^-_\epsilon(x)\f$.  Here,
      the lower \f$\epsilon\f$-binding set is defined as
      \f[
         \mathcal{B}^-_\epsilon(x) = \{\,\xi\in\Xi\,:\,x(\xi) \le a(\xi)+\epsilon_x,\;
                g(\xi) > \epsilon_g \,\}.
      \f]
      @param[out]      v    is the variable to be pruned.
      @param[in]       x    is the current optimization variable.
      @param[in]       g    is the negative search direction.
      @param[in]       xeps is the active-set tolerance \f$\epsilon_x\f$.
      @param[in]       geps is the binding-set tolerance \f$\epsilon_g\f$.
  */
  virtual void pruneLowerActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0) ) {
	const ThyraVector<Real>  & thyra_x = dynamic_cast<const ThyraVector<Real>&>(x);
  ThyraVector<Real>  & thyra_v = dynamic_cast<ThyraVector<Real>&>(v);
	const ThyraVector<Real>  & thyra_g = dynamic_cast<const ThyraVector<Real>&>(g);

	Real epsn = std::min(xeps,min_diff_);

  Thyra::ele_wise_prune_lower(*thyra_x.getVector(), *thyra_x_lo_, *thyra_g.getVector(), thyra_v.getVector().ptr(), epsn );
  }

  /** \brief Set the input vector to the upper bound.

      This function sets the input vector \f$u\f$ to the upper bound \f$b\f$.
      @param[out]    u   is the vector to be set to the upper bound.
  */
  virtual void setVectorToUpperBound( Vector<Real> &u ) {
	  ThyraVector<Real>  & thyra_up = dynamic_cast<ThyraVector<Real>&>(u);
	  thyra_up.set(ThyraVector<Real>(thyra_x_up_));
  }

  /** \brief Set the input vector to the lower bound.

      This function sets the input vector \f$l\f$ to the lower bound \f$a\f$.
      @param[out]    l   is the vector to be set to the lower bound.
  */
  virtual void setVectorToLowerBound( Vector<Real> &l ) {
	  ThyraVector<Real>  & thyra_lo = dynamic_cast<ThyraVector<Real>&>(l);
	  thyra_lo.set(ThyraVector<Real>(thyra_x_lo_));
  }

}; // class Thyra_BoundConstraint

} // namespace ROL

#endif
