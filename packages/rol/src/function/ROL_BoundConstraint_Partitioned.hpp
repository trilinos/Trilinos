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

#ifndef ROL_BOUND_CONSTRAINT_PARTITIONED_H
#define ROL_BOUND_CONSTRAINT_PARTITIONED_H

#include "ROL_BoundConstraint.hpp"
#include "ROL_Vector_Partitioned.hpp"
#include "ROL_Types.hpp"
#include <iostream>

/** @ingroup func_group
    \class ROL::BoundConstraint_Partitioned
    \brief A composite composite BoundConstraint formed
           from bound constraints on subvectors of a PartitionedVector

*/
namespace ROL {

template <class Real>
class BoundConstraint_Partitioned : public BoundConstraint<Real> {

  typedef Vector<Real>                          V;
  typedef PartitionedVector<Real>               PV;
  typedef typename std::vector<Real>::size_type uint;

private:
  std::vector<Teuchos::RCP<BoundConstraint<Real> > > bnd_;

  uint dim_;

public:
  ~BoundConstraint_Partitioned() {}

  /** \brief Default constructor.

      The default constructor automatically turns the constraints on.
  */
  BoundConstraint_Partitioned(std::vector<const Teuchos::RCP<BoundConstraint<Real> > > &bnd )
    : bnd_(bnd) {
    bool active = false;

    dim_ = bnd.size();

    for( uint k=0; k<dim_; ++k ) {
      active = active || bnd_[k]->isActivated();
    }    
  
    if( active ) {
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

    const PV &xpv = Teuchos::dyn_cast<const PV>(x);

    for( uint k=0; k<dim_; ++k ) {
      if( bnd_[k]->isActivate() {
        bnd_[k]->update(*(xpv.get(k)),flag,iter);   
      }
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

    PV &xpv = Teuchos::dyn_cast<PV>(x);
    for( uint k=0; k<dim_; ++k ) {
      Teuchos::RCP<V> xk = xpv.get(k)->clone(); 
      xk->set(*(xpv.get(k)));
      bnd_[k]->project(*xk);
      xpv.set(k,*xs);
    }

  /** \brief Determine if a vector of Lagrange multipliers is nonnegative components.
  
      This function returns true if components of \f$l\f$ corresponding to the components of \f$x\f$ 
      that are active at the upper bound are nonpositive or the components of \f$l\f$ corresponding
      to the components of \f$x\f$ that are active at the lower bound are nonnegative.
  */
  bool checkMultipliers( const Vector<Real> &l, const Vector<Real> &x ) {

    const PV &lpv = Teuchos::dyn_cast<const PV>(l);
    const PV &xpv = Teuchos::dyn_cast<const PV>(x);

    bool nonneg = true;

    for( uint k=0; k<dim_; ++k ) {
      if( bnd_[k]->isActivated() ) {
        nonneg = nonneg && bnd_[k]->checkMultipliers(*(lpv.get(k)),*(xpv.get(k)));
      }
    } 
    return nonneg;
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

          PV &vpv = Teuchos::dyn_cast<PV>(v);
    const PV &xpv = Teuchos::dyn_cast<const PV>(x);

    for( uint k=0; k<dim_; ++k ) {
      if( bnd_[k]->isActivated() ) {
        Teuchos::RCP<V> vk = vpv.get(k)->clone(); vk->set(*(vpv.get(k)));
        bnd_[k]->pruneUpperActive(*vk,*(xpv.get(k)),eps); 
        vpv.set(k,*vk);
      }
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

          PV &vpv = Teuchos::dyn_cast<PV>(v);
    const PV &gpv = Teuchos::dyn_cast<const PV>(g);
    const PV &xpv = Teuchos::dyn_cast<const PV>(x);

    for( uint k=0; k<dim_; ++k ) {
      Teuchos::RCP<V> vk = vpv.get(k)->clone(); vk->set(*(vpv.get(k)));
      bnd_[k]->pruneUpperActive(*vk,*(xpv.get(k)),eps);
      vpv.set(k,*vk);
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
     
          PV &vpv = Teuchos::dyn_cast<PV>(v);
    const PV &xpv = Teuchos::dyn_cast<const PV>(x);
 
   for( uint k=0; k<dim_; ++k ) {
      if( bnd_[k]->isActivated() ) {
        Teuchos::RCP<V> vk = vpv.get(k)->clone(); vk->set(*(vpv.get(k)));
        bnd_[k]->pruneLowerActive(*vk,*(xpv.get(k)),eps); 
        vpv.set(k,*vk);
      }
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
          PV &vpv = Teuchos::dyn_cast<PV>(v);
    const PV &gpv = Teuchos::dyn_cast<const PV>(g);
    const PV &xpv = Teuchos::dyn_cast<const PV>(x);

    for( uint k=0; k<dim_; ++k ) {
      Teuchos::RCP<V> vk = vpv.get(k)->clone(); vk->set(*(vpv.get(k)));
      bnd_[k]->pruneLowerActive(*vk,*(xpv.get(k)),eps);
      vpv.set(k,*vk);
    }
  }
 
  const Teuchos::RCP<Vector<Real> > getLowerVectorRCP( void ) const {

    std::vector<Teuchos::RCP<V> > vec;

    for( uint k=0; k<dim_; ++k ) {
      vec.push_back( bnd_[k]->getLowerVectorRCP() );
    }
    
    Teuchos::RCP<V> lp = Teuchos::rcp( new PV( vec ) );
    return lp;
  }
       
  

  const Teuchos::RCP<Vector<Real> > getUpperVectorRCP( void ) const {

    std::vector<Teuchos::RCP<V> > vec;

    for( uint k=0; k<dim_; ++k ) {
      vec.push_back( bnd_[k]->getUpperVectorRCP() );
    }
    
    Teuchos::RCP<V> up = Teuchos::rcp( new PV( vec ) );
    return up;

  }



  /** \brief Set the input vector to the upper bound.

      This function sets the input vector \f$u\f$ to the upper bound \f$b\f$.
      @param[out]    u   is the vector to be set to the upper bound.
  */ 
  void setVectorToUpperBound( Vector<Real> &u ) {
    PV &upv = Teuchos::dyn_cast<PV>(u);

    for( uint k=0; k<dim_; ++k ) {
      Teuchos::RCP<V> uk = upv.get(k)->clone();
      bnd_[k]->setVectorToUpperBound(*uk);
      upv.set(k,*uk);
    }  
  }

  /** \brief Set the input vector to the lower bound.

      This function sets the input vector \f$l\f$ to the lower bound \f$a\f$.
      @param[out]    l   is the vector to be set to the lower bound.
  */ 
  void setVectorToLowerBound( Vector<Real> &l ) {
    PV &lpv = Teuchos::dyn_cast<PV>(l);

    for( uint k=0; k<dim_; ++k ) {
      Teuchos::RCP<V> lk = lpv.get(k)->clone();
      bnd_[k]->setVectorToLowerBound(*uk);
      lpv.set(k,*lk);
    }  
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

          PV &vpv = Teuchos::dyn_cast<PV>(v);
    const PV &xpv = Teuchos::dyn_cast<const PV>(x); 
 
    for( uint k=0; k<dim_; ++k ) {
      Teuchos::RCP<V> vk = vpv.get(k)->clone();
      bnd_[k]->pruneActive(*vk,*(xpv.get(k)),eps);
      vpv.set(k,*vk);
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
 
          PV &vpv = Teuchos::dyn_cast<PV>(v);
    const PV &gpv = Teuchos::dyn_cast<const PV>(g);
    const PV &xpv = Teuchos::dyn_cast<const PV>(x);

    for( uint k=0; k<dim_; ++k ) {
      if(bnd_[k]->isActive()) {
        Teuchos::RCP<V> vk = vpv.get(k)->clone();
        bnd_[k]->pruneActive(*vk,*(gpv.get(k),*(xpv.get(k))),eps;
        vpv.set(k,*vk);
      }
    }
  }

  /** \brief Check if the vector, v, is feasible.

      This function returns true if \f$v = P_{[a,b]}(v)\f$.
      @param[in]    v   is the vector to be checked.
  */
  bool isFeasible( const Vector<Real> &v ) { 
    bool feasible = true;
    const PV &vs = Teuchos::dyn_cast<const PV>(v);
    for( uint k=0; k<dim_; ++k ) {
      feasible = feasible && bnd_[k]_->isFeasible(*(vs.get(k)));
    }
    return feasible;
  }

}; // class BoundConstraint_Partitioned

} // namespace ROL

#endif
