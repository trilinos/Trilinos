// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BOUND_CONSTRAINT_SIMOPT_H
#define ROL_BOUND_CONSTRAINT_SIMOPT_H

#include "ROL_BoundConstraint.hpp"
#include "ROL_Vector_SimOpt.hpp"
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
class BoundConstraint_SimOpt : public BoundConstraint<Real> {
private:
  const Ptr<BoundConstraint<Real>> bnd1_, bnd2_;

public:
  ~BoundConstraint_SimOpt() {}

  /** \brief Default constructor.

      The default constructor automatically turns the constraints on.
  */
  BoundConstraint_SimOpt(const Ptr<BoundConstraint<Real>> &bnd1,
                         const Ptr<BoundConstraint<Real>> &bnd2)
    : bnd1_(bnd1), bnd2_(bnd2) {
    if ( bnd1_->isActivated() || bnd2_->isActivated() )
      BoundConstraint<Real>::activate();
    else
      BoundConstraint<Real>::deactivate();
  }

  /** \brief Constructor for single bound constraint.

      Construct a SimOpt bound constraint with either a bound on
      the Sim or Opt variables.
      @param[in]     bnd is the bound constraint.
      @param[in]     optBnd is true if the bound is on the Opt variable.
  */
  BoundConstraint_SimOpt(const Ptr<BoundConstraint<Real>> &bnd,
                         bool optBnd = true)
    : bnd1_(optBnd ? makePtr<BoundConstraint<Real>>() : bnd),
      bnd2_(optBnd ? bnd : makePtr<BoundConstraint<Real>>()) {
    if ( bnd1_->isActivated() || bnd2_->isActivated() )
      BoundConstraint<Real>::activate();
    else
      BoundConstraint<Real>::deactivate();
  }

  /** \brief Constructor for single bound constraint.

      Construct a SimOpt bound constraint with either a bound on
      the Sim or Opt variables.
      @param[in]     bnd is the bound constraint.
      @param[in]     x is a Sim variable if optBnd, an Opt variable otherwise.
      @param[in]     optBnd is true if the bound is on the Opt variable.
  */
  BoundConstraint_SimOpt(const Ptr<BoundConstraint<Real>> &bnd,
                         const Vector<Real> &x,
                         bool optBnd = true)
    : bnd1_(optBnd ? makePtr<BoundConstraint<Real>>(x) : bnd),
      bnd2_(optBnd ? bnd : makePtr<BoundConstraint<Real>>(x)) {
    if ( bnd1_->isActivated() || bnd2_->isActivated() )
      BoundConstraint<Real>::activate();
    else
      BoundConstraint<Real>::deactivate();
  }

  /** \brief Project optimization variables onto the bounds.

      This function implements the projection of \f$x\f$ onto the bounds, i.e.,
      \f[
         (P_{[a,b]}(x))(\xi) = \min\{b(\xi),\max\{a(\xi),x(\xi)\}\} \quad \text{for almost every }\xi\in\Xi.
      \f]
       @param[in,out]      x is the optimization variable.
  */
  void project( Vector<Real> &x ) {
    Vector_SimOpt<Real> &xs = dynamic_cast<Vector_SimOpt<Real>&>(x);
    if (bnd1_->isActivated()) bnd1_->project(*(xs.get_1()));
    if (bnd2_->isActivated()) bnd2_->project(*(xs.get_2()));
  }

  /** \brief Project optimization variables into the interior of the feasible set.

      This function implements the projection of \f$x\f$ into the interior of the
      feasible set, i.e.,
      \f[
         (P_{[a,b]}(x))(\xi) \in (a(\xi),b(\xi))
             \quad \text{for almost every }\xi\in\Xi.
      \f]
       @param[in,out]      x is the optimization variable.
  */
  void projectInterior( Vector<Real> &x ) {
    Vector_SimOpt<Real> &xs = dynamic_cast<Vector_SimOpt<Real>&>(x);
    if (bnd1_->isActivated()) bnd1_->projectInterior(*(xs.get_1()));
    if (bnd2_->isActivated()) bnd2_->projectInterior(*(xs.get_2()));
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
  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &x, Real eps = Real(0) ) {
    Vector_SimOpt<Real> &vs = dynamic_cast<Vector_SimOpt<Real>&>(v);
    const Vector_SimOpt<Real> &xs = dynamic_cast<const Vector_SimOpt<Real>&>(x);
    if (bnd1_->isActivated()) bnd1_->pruneUpperActive(*(vs.get_1()),*(xs.get_1()),eps);
    if (bnd2_->isActivated()) bnd2_->pruneUpperActive(*(vs.get_2()),*(xs.get_2()),eps);
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
  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0) ) {
    Vector_SimOpt<Real> &vs = dynamic_cast<Vector_SimOpt<Real>&>(v);
    const Vector_SimOpt<Real> &gs = dynamic_cast<const Vector_SimOpt<Real>&>(g);
    const Vector_SimOpt<Real> &xs = dynamic_cast<const Vector_SimOpt<Real>&>(x);
    if (bnd1_->isActivated()) bnd1_->pruneUpperActive(*(vs.get_1()),*(gs.get_1()),*(xs.get_1()),xeps,geps);
    if (bnd2_->isActivated()) bnd2_->pruneUpperActive(*(vs.get_2()),*(gs.get_2()),*(xs.get_2()),xeps,geps);
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
  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &x, Real eps = Real(0) ) {
    Vector_SimOpt<Real> &vs = dynamic_cast<Vector_SimOpt<Real>&>(v);
    const Vector_SimOpt<Real> &xs = dynamic_cast<const Vector_SimOpt<Real>&>(x);
    if (bnd1_->isActivated()) bnd1_->pruneLowerActive(*(vs.get_1()),*(xs.get_1()),eps);
    if (bnd2_->isActivated()) bnd2_->pruneLowerActive(*(vs.get_2()),*(xs.get_2()),eps);
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
  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0) ) {
    Vector_SimOpt<Real> &vs = dynamic_cast<Vector_SimOpt<Real>&>(v);
    const Vector_SimOpt<Real> &gs = dynamic_cast<const Vector_SimOpt<Real>&>(g);
    const Vector_SimOpt<Real> &xs = dynamic_cast<const Vector_SimOpt<Real>&>(x);
    if (bnd1_->isActivated()) bnd1_->pruneLowerActive(*(vs.get_1()),*(gs.get_1()),*(xs.get_1()),xeps,geps);
    if (bnd2_->isActivated()) bnd2_->pruneLowerActive(*(vs.get_2()),*(gs.get_2()),*(xs.get_2()),xeps,geps);
  }

  const Ptr<const Vector<Real>> getLowerBound( void ) const {
    const Ptr<const Vector<Real>> l1 = bnd1_->getLowerBound();
    const Ptr<const Vector<Real>> l2 = bnd2_->getLowerBound();
    return makePtr<Vector_SimOpt<Real>>( constPtrCast<Vector<Real>>(l1),
                                         constPtrCast<Vector<Real>>(l2) );
  }

  const Ptr<const Vector<Real>> getUpperBound(void) const {
    const Ptr<const Vector<Real>> u1 = bnd1_->getUpperBound();
    const Ptr<const Vector<Real>> u2 = bnd2_->getUpperBound();
    return makePtr<Vector_SimOpt<Real>>( constPtrCast<Vector<Real>>(u1),
                                         constPtrCast<Vector<Real>>(u2) );
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
  void pruneActive( Vector<Real> &v, const Vector<Real> &x, Real eps = Real(0) ) {
    Vector_SimOpt<Real> &vs = dynamic_cast<Vector_SimOpt<Real>&>(v);
    const Vector_SimOpt<Real> &xs = dynamic_cast<const Vector_SimOpt<Real>&>(x);
    if (bnd1_->isActivated()) bnd1_->pruneActive(*(vs.get_1()),*(xs.get_1()),eps);
    if (bnd2_->isActivated()) bnd2_->pruneActive(*(vs.get_2()),*(xs.get_2()),eps);
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
  void pruneActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0) ) {
    Vector_SimOpt<Real> &vs = dynamic_cast<Vector_SimOpt<Real>&>(v);
    const Vector_SimOpt<Real> &gs = dynamic_cast<const Vector_SimOpt<Real>&>(g);
    const Vector_SimOpt<Real> &xs = dynamic_cast<const Vector_SimOpt<Real>&>(x);
    if (bnd1_->isActivated()) bnd1_->pruneActive(*(vs.get_1()),*(gs.get_1()),*(xs.get_1()),xeps,geps);
    if (bnd2_->isActivated()) bnd2_->pruneActive(*(vs.get_2()),*(gs.get_2()),*(xs.get_2()),xeps,geps);
  }

  /** \brief Check if the vector, v, is feasible.

      This function returns true if \f$v = P_{[a,b]}(v)\f$.
      @param[in]    v   is the vector to be checked.
  */
  bool isFeasible( const Vector<Real> &v ) {
    const Vector_SimOpt<Real> &vs = dynamic_cast<const Vector_SimOpt<Real>&>(v);
    return (bnd1_->isFeasible(*(vs.get_1()))) && (bnd2_->isFeasible(*(vs.get_2())));
  }

  /** \brief Apply inverse scaling function.

      This function applies the inverse scaling function \f$d(x,g)\f$ to
      a vector \f$v\f$, i.e., the output is \f$\mathrm{diag}(d(x,g)^{-1})v\f$.
      The scaling function must satisfy:
      (i) \f$d(x,g)_i = 0\f$ if \f$x_i = a_i\f$ and \f$g_i \ge 0\f$;
      (ii) \f$d(x,g)_i = 0\f$ if \f$x_i = b_i\f$ and \f$g_i \le 0\f$; and
      (iii) \f$d(x,g)_i > 0\f$ otherwise.
      @param[out] dv   is the inverse scaling function applied to v.
      @param[in]   v   is the vector being scaled.
      @param[in]   x   is the primal vector at which the scaling function is evaluated.
      @param[in]   g   is the dual vector at which the scaling function is evaluated.
  */
  void applyInverseScalingFunction(Vector<Real> &dv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &g) const{
    Vector_SimOpt<Real> &dvs = dynamic_cast<Vector_SimOpt<Real>&>(dv);
    const Vector_SimOpt<Real> &vs = dynamic_cast<const Vector_SimOpt<Real>&>(v);
    const Vector_SimOpt<Real> &xs = dynamic_cast<const Vector_SimOpt<Real>&>(x);
    const Vector_SimOpt<Real> &gs = dynamic_cast<const Vector_SimOpt<Real>&>(g);
    if (bnd1_->isActivated()) bnd1_->applyInverseScalingFunction(*(dvs.get_1()),*(vs.get_1()),*(xs.get_1()),*(gs.get_1()));
    if (bnd2_->isActivated()) bnd2_->applyInverseScalingFunction(*(dvs.get_2()),*(vs.get_2()),*(xs.get_2()),*(gs.get_2()));
  }

  /** \brief Apply scaling function Jacobian.

      This function applies the Jacobian of the scaling function \f$d(x,g)\f$ to
      a vector \f$v\f$.  The output is \f$\mathrm{diag}(d_x(x,g)g)v\f$.  The
      scaling function must satisfy:
      (i) \f$d(x,g)_i = 0\f$ if \f$x_i = a_i\f$ and \f$g_i \ge 0\f$;
      (ii) \f$d(x,g)_i = 0\f$ if \f$x_i = b_i\f$ and \f$g_i \le 0\f$; and
      (iii) \f$d(x,g)_i > 0\f$ otherwise.
      @param[out] dv   is the scaling function Jacobian applied to v.
      @param[in]   v   is the vector being scaled.
      @param[in]   x   is the primal vector at which the scaling function is evaluated.
      @param[in]   g   is the dual vector at which the scaling function is evaluated.
  */
  void applyScalingFunctionJacobian(Vector<Real> &dv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &g) const {
    Vector_SimOpt<Real> &dvs = dynamic_cast<Vector_SimOpt<Real>&>(dv);
    const Vector_SimOpt<Real> &vs = dynamic_cast<const Vector_SimOpt<Real>&>(v);
    const Vector_SimOpt<Real> &xs = dynamic_cast<const Vector_SimOpt<Real>&>(x);
    const Vector_SimOpt<Real> &gs = dynamic_cast<const Vector_SimOpt<Real>&>(g);
    if (bnd1_->isActivated()) bnd1_->applyScalingFunctionJacobian(*(dvs.get_1()),*(vs.get_1()),*(xs.get_1()),*(gs.get_1()));
    if (bnd2_->isActivated()) bnd2_->applyScalingFunctionJacobian(*(dvs.get_2()),*(vs.get_2()),*(xs.get_2()),*(gs.get_2()));
  }

}; // class BoundConstraint_SimOpt

} // namespace ROL

#endif
