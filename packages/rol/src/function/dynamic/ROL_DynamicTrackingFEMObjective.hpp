// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_DYNAMICTRACKINGFEMOBJECTIVE_HPP
#define ROL_DYNAMICTRACKINGFEMOBJECTIVE_HPP

#include "ROL_DynamicObjective.hpp"
#include "ROL_VectorWorkspace.hpp"


/** @ingroup func_group
    \class ROL::DynamicTrackingFEMObjective
    \brief Defines the time-local contribution to a quadratic tracking
           objective

           \f[ f_k(u,z) = \frac{1}{2} \int\limits_{t_{k-1}}^{t_k} \| u(t)-\tilde u(t)\|^2 
                                              + \alpha \|z(t)\^2\,\mathm{d}t \f]

           Currently approximates the state disparity norm using linear finite elements
           which couples the old and new state contributions
*/


namespace ROL {

template<typename Real>
class DynamicTrackingFEMObjective : public DynamicObjective<Real> {
public:

  using V         = Vector<Real>;
  using TS        = TimeStamp<Real>;
  using size_type = typename PartitionedVector<Real>::size_type;  

private:

  Ptr<PartitionedVector<Real>>  target_;

  size_type             Nt_;     // Number of time steps
  Real                  alpha_;  // Regularization parameter

  mutable VectorWorkspace<Real> workspace_;

public:

  DynamicTrackingFEMObjective( const Ptr<PartitionedVector<Real>>& target, Real alpha=0.0 ) : 
    target_(target), Nt_(target_->numVectors()), alpha_(alpha) {}
  
  virtual ~DynamicTrackingFEMObjective() {}  

  virtual Real value( const V& uo, const V& un, 
                      const V& z, const TS& timeStamp ) const override {

    Real dt = timeStamp.t.at(1)-timeStamp.t.at(0);
    Real w  = 2.0*dt/3.0;

    size_type k = timeStamp.k;     
 
    auto res_n = workspace_.copy(un);
    auto res_o = workspace_.copy(uo);

    Real result = 0.5*dt*alpha_*z.dot(z);

    res_n->set(un);
    res_n->axpy( -1.0, *(target_->get(k)) );
    result += w*res_n->dot(*res_n);

    if( k>0 ) {
      res_o->set(uo);
      res_o->axpy( -1, *(target_->get(k-1) ) );
      result += w*res_n->dot(*res_o);
      result += w*res_o->dot(*res_o);
    }

    return result;
  }

  //----------------------------------------------------------------------------
  // Gradient Terms
  virtual void gradient_uo( V& g, const V& uo, const V& un, 
                            const V& z, const TS& timeStamp ) const override {
    Real dt = timeStamp.t.at(1)-timeStamp.t.at(0);
    Real w  = dt/3.0;

    size_type k = timeStamp.k;     
 
    auto res_n = workspace_.copy(un);
    auto res_o = workspace_.copy(uo);

    


    else g.zero();
  }

  virtual void gradient_un( V& g, const V& uo, const V& un, 
                            const V& z, const TS& timeStamp ) const override {
    Real dt = timeStamp.t.at(1)-timeStamp.t.at(0);
    g.set(un);
    g.axpy(-1.0, *(target_->get(timeStamp.k)) );
    g.scale(0.5*dt);
  }

  virtual void gradient_z( V& g, const V& uo, const V& un, 
                           const V& z, const TS& timeStamp ) const override {
    Real dt = timeStamp.t.at(1)-timeStamp.t.at(0);
    g.set(z);
    g.scale(dt*alpha_);
  }

  //----------------------------------------------------------------------------
  // Hessian-Vector product terms
  virtual void hessVec_uo_uo( V& hv, const V& v, const V& uo, const V& un, 
                              const V& z, const TS& timeStamp ) const override {
    Real dt = timeStamp.t.at(1)-timeStamp.t.at(0);
    if( timeStamp.k>0 ) {
      hv.set(v);
      hv.axpy(-1.0, *(target_->get(timeStamp.k-1)) );
      hv.scale(0.5*dt);
    }
    else hv.zero();
  }

  virtual void hessVec_uo_un( V& hv, const V& v, const V& uo, const V& un, 
                              const V& z, const TS& timeStamp ) const override { 
    hv.zero(); 
  }

  virtual void hessVec_uo_z( V& hv, const V& v, const V& uo, const V& un, 
                             const V& z, const TS& timeStamp ) const override { 
    hv.zero(); 
  }

  virtual void hessVec_un_uo( V& hv, const V& v, const V& uo, const V& un, 
                              const V& z, const TS& timeStamp ) const override {
    hv.zero();
  }

  virtual void hessVec_un_un( V& hv, const V& v, const V& uo, const V& un, 
                              const V& z, const TS& timeStamp ) const override {
    hv.set(v);
    hv.scale(0.5*(timeStamp.t.at(1)-timeStamp.t.at(0)));
  }

  virtual void hessVec_un_z( V& hv, const V& v, const V& uo, const V& un, 
                             const V& z, const TS& timeStamp ) const override {
    hv.zero();
  }

  virtual void hessVec_z_uo( V& hv, const V& v, const V& uo, const V& un, 
                              const V& z, const TS& timeStamp ) const override {
    hv.zero();
  }

  virtual void hessVec_z_un( V& hv, const V& v, const V& uo, const V& un, 
                              const V& z, const TS& timeStamp ) const override {
    hv.zero();
  }

  virtual void hessVec_z_z( V& hv, const V& v, const V& uo, const V& un, 
                             const V& z, const TS& timeStamp ) const override {
    hv.set(v);
    hv.scale(alpha_*(timeStamp.t.at(1)-timeStamp.t.at(0)));
  }

}; // DynamicTrackingFEMObjective


template<typename Real>
inline Ptr<DynamicObjective<Real>> 
make_DynamicTrackingFEMObjective( const Ptr<PartitionedVector<Real>>& target, Real alpha=0.0 ) {
  Ptr<DynamicObjective<Real>> obj = makePtr<DynamicTrackingFEMObjective<Real>>(target,alpha);
  return obj;
}

} // namespace ROL

#endif  // ROL_DYNAMICTRACKINGFEMOBJECTIVE_HPP


