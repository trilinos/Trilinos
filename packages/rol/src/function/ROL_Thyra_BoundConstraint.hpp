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

#ifndef ROL_THYRA_BOUND_CONSTRAINT_H
#define ROL_THYRA_BOUND_CONSTRAINT_H

#include "ROL_BoundConstraint.hpp"
#include "TOpEleWiseBoundsHelpers.hpp"

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

  virtual bool isFeasible( const Vector<Real> &v ) {
	  bool check = BoundConstraint<Real>::isFeasible(v);
	  if(!check)
	    std::cout << "It's not feasible!!" << std::endl;
	  return check;
  }

  Thyra_BoundConstraint(Teuchos::RCP<Thyra::VectorBase<Real> > p_min, Teuchos::RCP<Thyra::VectorBase<Real> > p_max, Real min_diff) : BoundConstraint<Real>(), thyra_x_lo_(p_min), thyra_x_up_(p_max), min_diff_(min_diff) {};

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
	  ThyraVector<Real>  & thyra_x = Teuchos::dyn_cast<ThyraVector<Real> >(x);

	  //Teuchos::RCP<Vector<Real> > x_lo = x.clone();
	  //Teuchos::RCP<Vector<Real> > x_up = x.clone();
	//  ThyraVector<Real> thyra_x_lo = Teuchos::dyn_cast<ThyraVector<Real> >(*x_lo);
	//  ThyraVector<Real> thyra_x_up = Teuchos::dyn_cast<ThyraVector<Real> >(*x_up);
	//  thyra_x_lo.putScalar(p_min_);
	//  thyra_x_up.putScalar(p_max_);
	  Thyra::ele_wise_bound(*thyra_x_lo_, *thyra_x_up_, thyra_x.getNonConstVector().ptr() );

//	  {
//		  Thyra::DetachedVectorView<Real> thyra_x_view(thyra_x.getNonConstVector());
//		  for (int i = 0; i < thyra_x_view.subDim(); ++i)
//			  thyra_x_view[i] = std::max(p_min_,std::min(p_max_,thyra_x_view[i]));
//	  }
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
	  const ThyraVector<Real>  & thyra_x = Teuchos::dyn_cast<const ThyraVector<Real> >(x);
	  ThyraVector<Real>  & thyra_v = Teuchos::dyn_cast<ThyraVector<Real> >(v);

	Real epsn = std::min(eps,min_diff_);

//  Teuchos::RCP<Vector<Real> > x_up = x.clone();
//  ThyraVector<Real> thyra_x_up = Teuchos::dyn_cast<ThyraVector<Real> >(*x_up);
//  thyra_x_up.putScalar(p_max_);
  Thyra::ele_wise_prune_upper(*thyra_x.getVector(), *thyra_x_up_, thyra_v.getNonConstVector().ptr(), epsn );


/*	{
      Thyra::DetachedVectorView<Real> thyra_x_view(thyra_x.getNonConstVector());
	  Thyra::DetachedVectorView<Real> thyra_v_view(thyra_v.getNonConstVector());

	  for (int i = 0; i < thyra_v_view.subDim(); ++i) {
		if (thyra_x_view[i] >= p_max_-epsn)
		  thyra_v_view[i] = 0.0;
	  }
	}
	*/
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
	const ThyraVector<Real>  & thyra_x = Teuchos::dyn_cast<const ThyraVector<Real> >(x);
	ThyraVector<Real>  & thyra_v = Teuchos::dyn_cast<ThyraVector<Real> >(v);
	const ThyraVector<Real>  & thyra_g = Teuchos::dyn_cast<const ThyraVector<Real> >(g);

	Real epsn = std::min(eps,min_diff_);

  //Teuchos::RCP<Vector<Real> > x_up = x.clone();
 // ThyraVector<Real> thyra_x_up = Teuchos::dyn_cast<ThyraVector<Real> >(*x_up);
  //thyra_x_up.putScalar(p_max_);
  Thyra::ele_wise_prune_upper(*thyra_x.getVector(), *thyra_x_up_, *thyra_g.getVector(), thyra_v.getNonConstVector().ptr(), epsn );

/*	{
	  Thyra::DetachedVectorView<Real> thyra_x_view(thyra_x.getNonConstVector());
	  Thyra::DetachedVectorView<Real> thyra_v_view(thyra_v.getNonConstVector());
	  Thyra::DetachedVectorView<Real> thyra_g_view(thyra_g.getNonConstVector());
	  for (int i = 0; i < thyra_v_view.subDim(); ++i) {
		if ( (thyra_x_view[i] >= p_max_-epsn) && (thyra_g_view[i] < 0.0))
		  thyra_v_view[i] = 0.0;
	  }
	}
	*/
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
	const ThyraVector<Real>  & thyra_x = Teuchos::dyn_cast<const ThyraVector<Real> >(x);
  ThyraVector<Real>  & thyra_v = Teuchos::dyn_cast<ThyraVector<Real> >(v);

	Real epsn = std::min(eps,min_diff_);

  //Teuchos::RCP<Vector<Real> > x_lo = x.clone();
 // ThyraVector<Real> thyra_x_lo = Teuchos::dyn_cast<ThyraVector<Real> >(*x_lo);
 // thyra_x_lo.putScalar(p_min_);
  Thyra::ele_wise_prune_lower(*thyra_x.getVector(), *thyra_x_lo_, thyra_v.getNonConstVector().ptr(), epsn );
/*	{
	  Thyra::DetachedVectorView<Real> thyra_x_view(thyra_x.getNonConstVector());
	  Thyra::DetachedVectorView<Real> thyra_v_view(thyra_v.getNonConstVector());
	  for (int i = 0; i < thyra_v_view.subDim(); ++i) {
	    if (thyra_x_view[i] <= p_min_+epsn)
	      thyra_v_view[i] = 0.0;
	  	}
	}*/
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
	const ThyraVector<Real>  & thyra_x = Teuchos::dyn_cast<const ThyraVector<Real> >(x);
  ThyraVector<Real>  & thyra_v = Teuchos::dyn_cast<ThyraVector<Real> >(v);
	const ThyraVector<Real>  & thyra_g = Teuchos::dyn_cast<const ThyraVector<Real> >(g);

	Real epsn = std::min(eps,min_diff_);

 // Teuchos::RCP<Vector<Real> > x_lo = x.clone();
 // ThyraVector<Real> thyra_x_lo = Teuchos::dyn_cast<ThyraVector<Real> >(*x_lo);
 // thyra_x_lo.putScalar(p_min_);
  Thyra::ele_wise_prune_lower(*thyra_x.getVector(), *thyra_x_lo_, *thyra_g.getVector(), thyra_v.getNonConstVector().ptr(), epsn );

	/*{
	  Thyra::DetachedVectorView<Real> thyra_x_view(thyra_x.getNonConstVector());
	  Thyra::DetachedVectorView<Real> thyra_v_view(thyra_v.getNonConstVector());
      Thyra::DetachedVectorView<Real> thyra_g_view(thyra_g.getNonConstVector());
      for (int i = 0; i < thyra_v_view.subDim(); ++i) {
		if ( (thyra_x_view[i] <= p_min_+epsn) && (thyra_g_view[i] > 0.0))
		  thyra_v_view[i] = 0.0;
	  }
	}*/
  }

  /** \brief Set the input vector to the upper bound.

      This function sets the input vector \f$u\f$ to the upper bound \f$b\f$.
      @param[out]    u   is the vector to be set to the upper bound.
  */
  virtual void setVectorToUpperBound( Vector<Real> &u ) {
	  ThyraVector<Real>  & thyra_up = Teuchos::dyn_cast<ThyraVector<Real> >(u);
	  thyra_up.set(ThyraVector<Real>(thyra_x_up_));
  }

  /** \brief Set the input vector to the lower bound.

      This function sets the input vector \f$l\f$ to the lower bound \f$a\f$.
      @param[out]    l   is the vector to be set to the lower bound.
  */
  virtual void setVectorToLowerBound( Vector<Real> &l ) {
	  ThyraVector<Real>  & thyra_lo = Teuchos::dyn_cast<ThyraVector<Real> >(l);
	  thyra_lo.set(ThyraVector<Real>(thyra_x_lo_));
  }

}; // class Thyra_BoundConstraint

} // namespace ROL

#endif
