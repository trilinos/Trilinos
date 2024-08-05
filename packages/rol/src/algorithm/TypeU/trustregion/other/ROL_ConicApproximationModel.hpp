// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_CONICAPPROXIMATIONMODEL_H
#define ROL_CONICAPPROXIMATIONMODEL_H

#include "ROL_Objective.hpp"
#include "ROL_VectorWorkspace.hpp"

/** @ingroup func_group
    \class ConicApproximationModel
    \brief Provides the interface to evaluate conic approximation function


    \f[ \psi_k(s) = f(x_k) + \frac{s^\top\nabla f(x_k)}{1-a^\top s}
                    +\frac{1}{2}\frac{s^\top \nabla^2 f(x_k) s}{(1-a^\top s)^2}
    \f]

    -----
*/


namespace ROL {

template <class Real>
class ConicApproximationModel : public Objective<Real> {

  using V   = Vector<Real>;
  using Obj = Objective<Real>;

private:
  Ptr<Obj>           obj_;
  const Ptr<const V> x_, a_;
  Ptr<V>             g_, s_, Hs_;
  Real               f_, gamma_, sHs_;
  
  VectorWorkspace<Real> workspace_;

public:

  virtual ~ConicApproximationModel() {}

  ConicApproximationModel( const Ptr<Obj>& obj, const Ptr<const V>& x, 
                           const Ptr<V>& s,     const Ptr<const V>& a ) :
    obj_( obj ), x_( x ),  a_( a ), g_( x_->dual().clone() ), s_( s ),
    Hs_( x->dual().clone() ) {
    Real tol = sqrt(ROL_EPSILON<Real>());
    gamma_ = 1.0-a_->dot(*s_);
    f_ = obj_->value( *x_, tol );
    obj_->gradient( *g_,*x,tol );
    obj_->hessVec( *Hs_, *s_, *x_, tol );
  }

  virtual void update( const V& s, bool flag=true, int iter=-1 ) override {
    Real tol = sqrt(ROL_EPSILON<Real>());
    s_->set(s);
    gamma_ = 1.0-a_->dot(*s_);
    obj_->hessVec( *Hs_, *s_, *x_, tol );
  }

  virtual Real value( const V& s, Real& tol ) override {
    return f_ + ( g_->dot(*s_) + 0.5*sHs_/gamma_ )/gamma_;
  }

  virtual void gradient( V &g, const V &s, Real &tol ) override {

    g.set( *g_ );         // g0
    g.scale( gamma_ );    // gamma*g0
    g.plus( *Hs_ );       // gamma*g0 + Hs

    auto u = workspace_.copy(*a_);
    u->scale( s_->dot(g) );
    g.scale( gamma_ );
    g.plus( *u );
    g.scale( std::pow(gamma_,-3) );

  }

  virtual void hessVec( V &hv, const V &v, const V &s, Real &tol ) override {

    auto u = workspace_.copy(v);
      
    u->scale( gamma_ );                 // gamma*v
    u->axpy( a_->dot(v), s );           // gamma*v + (a,v)*s
    obj_->hessVec( hv, *u, *x_, tol );  // gamma*Hv + (a,v)*Hs 
    hv.set(*u);
    hv.scale( gamma_ );
    hv.axpy(u->dot(s),*a_);
    hv.scale(std::pow( gamma_ ,-4));
  }
 
  virtual void invHessVec( V& hv, const V& v, const V& s, Real& tol ) override {

    auto u = workspace_.copy(v);
      
    u->axpy( -a_->dot(v), s );             //  v - (a,v)*s
    obj_->invHessVec( hv, *u, *x_, tol );  // Hv - (a,v)*Hs 
    hv.set(*u);
    hv.axpy(-u->dot(*s_),*a_);            
    hv.scale(std::pow(gamma_,2));
  }

  virtual void precond( V& Pv, const V& v, const V& s, Real &tol ) override {

    auto u = workspace_.copy(v);
      
    u->axpy( -a_->dot(v), *s_ );             //  v - (a,v)*s
    obj_->precond( Pv, *u, *x_, tol );  // Hv - (a,v)*Hs 
    Pv.set(*u);
    Pv.axpy(-u->dot(s),*a_);            
    Pv.scale(std::pow(gamma_,2));
}



}; // class ConicApproximationModel

} // namespace ROL


#endif
