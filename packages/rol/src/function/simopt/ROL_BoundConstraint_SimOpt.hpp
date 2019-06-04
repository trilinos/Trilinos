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
  ROL::Ptr<BoundConstraint<Real> > bnd1_;
  ROL::Ptr<BoundConstraint<Real> > bnd2_;

public:
  ~BoundConstraint_SimOpt() {}

  /** \brief Default constructor.

      The default constructor automatically turns the constraints on.
  */
  BoundConstraint_SimOpt(const ROL::Ptr<BoundConstraint<Real> > &bnd1,
                         const ROL::Ptr<BoundConstraint<Real> > &bnd2)
    : bnd1_(bnd1), bnd2_(bnd2) {
    if ( bnd1_->isActivated() || bnd2_->isActivated() ) {
      BoundConstraint<Real>::activate();
    }
    else {
      BoundConstraint<Real>::deactivate();
    }
  }

  /** \brief Update bounds. 

      The update function allows the user to update the bounds at each new iterations. 
          @param[in]      x      is the optimization variable.
          @param[in]      flag   is set to true if control is changed.
          @param[in]      iter   is the outer algorithm iterations count.
  */
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    const ROL::Vector_SimOpt<Real> &xs = dynamic_cast<const ROL::Vector_SimOpt<Real>&>(
      dynamic_cast<const ROL::Vector<Real>&>(x));
    if ( bnd1_->isActivated() ) {
      bnd1_->update(*(xs.get_1()),flag,iter);
    }
    if ( bnd2_->isActivated() ) {
      bnd2_->update(*(xs.get_2()),flag,iter);
    }
  }

  /** \brief Project optimization variables onto the bounds.

      This function implements the projection of \f$x\f$ onto the bounds, i.e., 
      \f[
         (P_{[a,b]}(x))(\xi) = \min\{b(\xi),\max\{a(\xi),x(\xi)\}\} \quad \text{for almost every }\xi\in\Xi. 
      \f]
       @param[in,out]      x is the optimization variable.
  */
  void project( Vector<Real> &x ) {
    ROL::Vector_SimOpt<Real> &xs = dynamic_cast<ROL::Vector_SimOpt<Real>&>(
      dynamic_cast<ROL::Vector<Real>&>(x));
    if ( bnd1_->isActivated() ) {
      ROL::Ptr<Vector<Real> > x1 = xs.get_1()->clone(); x1->set(*(xs.get_1()));
      bnd1_->project(*x1);
      xs.set_1(*x1);
    }
    if ( bnd2_->isActivated() ) {
      ROL::Ptr<Vector<Real> > x2 = xs.get_2()->clone(); x2->set(*(xs.get_2()));
      bnd2_->project(*x2);
      xs.set_2(*x2);
    }
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
    ROL::Vector_SimOpt<Real> &xs = dynamic_cast<ROL::Vector_SimOpt<Real>&>(
      dynamic_cast<ROL::Vector<Real>&>(x));
    if ( bnd1_->isActivated() ) {
      ROL::Ptr<Vector<Real> > x1 = xs.get_1()->clone(); x1->set(*(xs.get_1()));
      bnd1_->projectInterior(*x1);
      xs.set_1(*x1);
    }
    if ( bnd2_->isActivated() ) {
      ROL::Ptr<Vector<Real> > x2 = xs.get_2()->clone(); x2->set(*(xs.get_2()));
      bnd2_->projectInterior(*x2);
      xs.set_2(*x2);
    }
  }

  /** \brief Determine if a vector of Lagrange multipliers is nonnegative components.
  
      This function returns true if components of \f$l\f$ corresponding to the components of \f$x\f$ 
      that are active at the upper bound are nonpositive or the components of \f$l\f$ corresponding
      to the components of \f$x\f$ that are active at the lower bound are nonnegative.
  */
  bool checkMultipliers( const Vector<Real> &l, const Vector<Real> &x ) {
    const ROL::Vector_SimOpt<Real> &ls = dynamic_cast<const ROL::Vector_SimOpt<Real>&>(
      dynamic_cast<const ROL::Vector<Real>&>(l));
    const ROL::Vector_SimOpt<Real> &xs = dynamic_cast<const ROL::Vector_SimOpt<Real>&>(
      dynamic_cast<const ROL::Vector<Real>&>(x));
    bool nn1 = true;
    if ( bnd1_->isActivated() ) {
      nn1 = bnd1_->checkMultipliers(*(ls.get_1()),*(xs.get_1()));
    }
    bool nn2 = true;
    if ( bnd2_->isActivated() ) {
      nn2 = bnd2_->checkMultipliers(*(ls.get_2()),*(xs.get_2()));
    }
    return (nn1 && nn2);
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
  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0.0 ) {
    ROL::Vector_SimOpt<Real> &vs = dynamic_cast<ROL::Vector_SimOpt<Real>&>(
      dynamic_cast<ROL::Vector<Real>&>(v));
    const ROL::Vector_SimOpt<Real> &xs = dynamic_cast<const ROL::Vector_SimOpt<Real>&>(
      dynamic_cast<const ROL::Vector<Real>&>(x));
    if ( bnd1_->isActivated() ) {
      ROL::Ptr<Vector<Real> > v1 = vs.get_1()->clone(); v1->set(*(vs.get_1()));
      bnd1_->pruneUpperActive(*v1,*(xs.get_1()),eps);
      vs.set_1(*v1);
    }
    if ( bnd2_->isActivated() ) {
      ROL::Ptr<Vector<Real> > v2 = vs.get_2()->clone(); v2->set(*(vs.get_2()));
      bnd2_->pruneUpperActive(*v2,*(xs.get_2()),eps);
      vs.set_2(*v2);
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
      @param[in]       x   is the current optimization variable.
      @param[in]       g   is the negative search direction.
      @param[in]       eps is the active-set tolerance \f$\epsilon\f$.
  */
  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0.0 ) {
    ROL::Vector_SimOpt<Real> &vs = dynamic_cast<ROL::Vector_SimOpt<Real>&>(
      dynamic_cast<ROL::Vector<Real>&>(v));
    const ROL::Vector_SimOpt<Real> &gs = dynamic_cast<const ROL::Vector_SimOpt<Real>&>(
      dynamic_cast<const ROL::Vector<Real>&>(g));
    const ROL::Vector_SimOpt<Real> &xs = dynamic_cast<const ROL::Vector_SimOpt<Real>&>(
      dynamic_cast<const ROL::Vector<Real>&>(x));
    if ( bnd1_->isActivated() ) {
      ROL::Ptr<Vector<Real> > v1 = vs.get_1()->clone(); v1->set(*(vs.get_1()));
      bnd1_->pruneUpperActive(*v1,*(gs.get_1()),*(xs.get_1()),eps);
      vs.set_1(*v1);
    }
    if ( bnd2_->isActivated() ) {
      ROL::Ptr<Vector<Real> > v2 = vs.get_2()->clone(); v2->set(*(vs.get_2()));
      bnd2_->pruneUpperActive(*v2,*(gs.get_2()),*(xs.get_2()),eps);
      vs.set_2(*v2);
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
  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0.0 ) {
    ROL::Vector_SimOpt<Real> &vs = dynamic_cast<ROL::Vector_SimOpt<Real>&>(
      dynamic_cast<ROL::Vector<Real>&>(v));
    const ROL::Vector_SimOpt<Real> &xs = dynamic_cast<const ROL::Vector_SimOpt<Real>&>(
      dynamic_cast<const ROL::Vector<Real>&>(x));
    if ( bnd1_->isActivated() ) {
      ROL::Ptr<Vector<Real> > v1 = vs.get_1()->clone(); v1->set(*(vs.get_1()));
      bnd1_->pruneLowerActive(*v1,*(xs.get_1()),eps);
      vs.set_1(*v1);
    }
    if ( bnd2_->isActivated() ) {
      ROL::Ptr<Vector<Real> > v2 = vs.get_2()->clone(); v2->set(*(vs.get_2()));
      bnd2_->pruneLowerActive(*v2,*(xs.get_2()),eps);
      vs.set_2(*v2);
    }
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
  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0.0 ) {
    ROL::Vector_SimOpt<Real> &vs = dynamic_cast<ROL::Vector_SimOpt<Real>&>(
      dynamic_cast<ROL::Vector<Real>&>(v));
    const ROL::Vector_SimOpt<Real> &gs = dynamic_cast<const ROL::Vector_SimOpt<Real>&>(
      dynamic_cast<const ROL::Vector<Real>&>(g));
    const ROL::Vector_SimOpt<Real> &xs = dynamic_cast<const ROL::Vector_SimOpt<Real>&>(
      dynamic_cast<const ROL::Vector<Real>&>(x));
    if ( bnd1_->isActivated() ) {
      ROL::Ptr<Vector<Real> > v1 = vs.get_1()->clone(); v1->set(*(vs.get_1()));
      bnd1_->pruneLowerActive(*v1,*(gs.get_1()),*(xs.get_1()),eps);
      vs.set_1(*v1);
    }
    if ( bnd2_->isActivated() ) {
      ROL::Ptr<Vector<Real> > v2 = vs.get_2()->clone(); v2->set(*(vs.get_2()));
      bnd2_->pruneLowerActive(*v2,*(gs.get_2()),*(xs.get_2()),eps);
      vs.set_2(*v2);
    }
  }
 
  const ROL::Ptr<const Vector<Real>> getLowerBound( void ) const {
    const ROL::Ptr<const Vector<Real>> l1 = bnd1_->getLowerBound();
    const ROL::Ptr<const Vector<Real>> l2 = bnd2_->getLowerBound();
    return ROL::makePtr<Vector_SimOpt<Real>>( ROL::constPtrCast<Vector<Real>>(l1),
                                                 ROL::constPtrCast<Vector<Real>>(l2) );
  }

  const ROL::Ptr<const Vector<Real>> getUpperBound(void) const {
    const ROL::Ptr<const Vector<Real>> u1 = bnd1_->getUpperBound();
    const ROL::Ptr<const Vector<Real>> u2 = bnd2_->getUpperBound();
    return ROL::makePtr<Vector_SimOpt<Real>>( ROL::constPtrCast<Vector<Real>>(u1),
                                                 ROL::constPtrCast<Vector<Real>>(u2) );
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
  void pruneActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0.0 ) {
    ROL::Vector_SimOpt<Real> &vs = dynamic_cast<ROL::Vector_SimOpt<Real>&>(
      dynamic_cast<ROL::Vector<Real>&>(v));
    const ROL::Vector_SimOpt<Real> &xs = dynamic_cast<const ROL::Vector_SimOpt<Real>&>(
      dynamic_cast<const ROL::Vector<Real>&>(x));
    if ( bnd1_->isActivated() ) {
      ROL::Ptr<Vector<Real> > v1 = vs.get_1()->clone(); v1->set(*(vs.get_1()));
      bnd1_->pruneActive(*v1,*(xs.get_1()),eps);
      vs.set_1(*v1);
    }
    if ( bnd2_->isActivated() ) {
      ROL::Ptr<Vector<Real> > v2 = vs.get_2()->clone(); v2->set(*(vs.get_2()));
      bnd2_->pruneActive(*v2,*(xs.get_2()),eps);
      vs.set_2(*v2);
    }
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
  void pruneActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0.0 ) {
    ROL::Vector_SimOpt<Real> &vs = dynamic_cast<ROL::Vector_SimOpt<Real>&>(v);
    const ROL::Vector_SimOpt<Real> &gs = dynamic_cast<const ROL::Vector_SimOpt<Real>&>(g);
    const ROL::Vector_SimOpt<Real> &xs = dynamic_cast<const ROL::Vector_SimOpt<Real>&>(x);
    if ( bnd1_->isActivated() ) {
      ROL::Ptr<Vector<Real> > v1 = vs.get_1()->clone(); v1->set(*(vs.get_1()));
      bnd1_->pruneActive(*v1,*(gs.get_1()),*(xs.get_1()),eps);
      vs.set_1(*v1);
    }
    if ( bnd2_->isActivated() ) {
      ROL::Ptr<Vector<Real> > v2 = vs.get_2()->clone(); v2->set(*(vs.get_2()));
      bnd2_->pruneActive(*v2,*(gs.get_2()),*(xs.get_2()),eps);
      vs.set_2(*v2);
    }
  }

  /** \brief Check if the vector, v, is feasible.

      This function returns true if \f$v = P_{[a,b]}(v)\f$.
      @param[in]    v   is the vector to be checked.
  */
  bool isFeasible( const Vector<Real> &v ) { 
    const ROL::Vector_SimOpt<Real> &vs = dynamic_cast<const ROL::Vector_SimOpt<Real>&>(v);
    return (bnd1_->isFeasible(*(vs.get_1()))) && (bnd2_->isFeasible(*(vs.get_2())));
  }

}; // class BoundConstraint

} // namespace ROL

#endif
