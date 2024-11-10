// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_SEMISMOOTHNEWTONDUALMODEL_HPP
#define ROL_SEMISMOOTHNEWTONDUALMODEL_HPP

#include "ROL_TrustRegionModel.hpp"
#include "ROL_InactiveSetVector.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_VectorWorkspace.hpp"

/** @ingroup func_group
    \class ROL::SemismoothNewtonDualModel 
    \brief Implements the dual variable model function for a 
           semismooth Newton step.
    
    Reference: 
    Konstantin Pieper dissertation "Finite element discretization and efficient 
    numerical solution of elliptic and parabolic sparse control problems."
    http://nbn-resolving.de/urn/resolver.pl?urn:nbn:de:bvb:91-diss-20150420-1241413-1-4
    -----
*/


namespace ROL {


template<class Real>
class SemismoothNewtonDualModel : public TrustRegionModel<Real> {

  using V   = Vector<Real>;
  using VPrim = InactiveSet_PrimalVector<Real>;
  using VDual = InactiveSet_DualVector<Real>;

  using Obj = Objective<Real>;
  using Sec = Secant<Real>;
  using Bnd = BoundConstraint<Real>;
  
private:

  class ProjectedObjective : public Objective<Real> {
  private:
    Obj&   objPrimal_;
    Bnd&   bnd_;
    Ptr<V> primalVec_;

  public:
    ProjectedObjective( Obj& objPrimal, Bnd& bnd, const Ptr<V>& primalVec ) : 
      objPrimal_(objPrimal), bnd_(bnd), primalVec_( primalVec ) {}

      Real value( const V& p, Real& tol ) override {
        primalVec_->set(p);
        bnd_.project(*primalVec_);
        return objPrimal_->value(*primalVec_, tol);
      }
  
      void gradient( V& g, const V& p, Real& tol ) override {
        primalVec_->set(p);
        bnd_.project(*primalVec_);
        objPrimal_->gradient(g,*primalVec_, tol);
      }

      void hessVec( V& hv, const V& v, const V& p, Real& tol ) override {
        primalVec_->set(p);
        bnd_.project(*primalVec_);
        objPrimal_->hessVec(hv,v,*primalVec_, tol);
      }

  }; // ProjectedObjective 

  ProjectedObjective   projObj_;
  Bnd                  bnd_;
  Sec                  secant_;
  Ptr<V>               p_, g_, x_;
  Ptr<V>               ones_;
  Ptr<VPrim>           s_;    
  Real                 alpha_;

  VectorWorkspace<Real> workspace_;


public:

  SemismoothNewtonDualModel( Obj& obj, Bnd& bnd, const V& p, const V& g, const Real alpha ) :
    TrustRegionModel( obj, p, g, false ), bnd_( bnd ), 
    p_( p.clone() ), g_( p.dual().clone() ), x_( p.clone() ), ones_( p.clone() ),
    s_( p.clone(), ones_, p_, bnd_ ), projObj_( obj, bnd, p_ ), alpha_(alpha) {

    ones_->setScalar( Real(1.0) );
  }
    

  Real value( const V& s, Real& tol ) {

    auto hs = workspace_.clone(*g_);

    gradient(*g_,s,tol);
    hessVec(*hs,s,s,tol);
    hs->scale( 0.5 );
    hs->plus(*g_);
    s_->set(s);
    return s_->dot(*hs);                            
  } 

  void gradient( V& g, const V& s, Real& tol ) {
    projObj_->gradient(g,*p_,tol);
    g.axpy(alpha_,*p_);   
  }
  
  void hessVec( V& hv, const V& v, const V& s, Real& tol ) {
    auto vprune_ = workspace_.copy(v);
    bnd_->pruneActive( *vprune_, *p_ );
    projObj_->hessVec( hv, *vprune_, *p_, tol );
    hv.axpy(alpha_,v);
  }

  void update( const V& p, bool flag = true, int iter = -1 ) {
    p_->set(p);
    auto x = this->getIterate();
  }

} // namespace ROL

