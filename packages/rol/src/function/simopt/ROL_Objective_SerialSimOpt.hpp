// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_OBJECTIVE_SERIALSIMOPT_HPP
#define ROL_OBJECTIVE_SERIALSIMOPT_HPP

#include "ROL_Objective_TimeSimOpt.hpp"

namespace ROL {

template<typename Real>
class Objective_SerialSimOpt : public Objective_SimOpt<Real> {

  using V   = Vector<Real>;
  using PV  = PartitionedVector<Real>;
  using Obj = Objective_TimeSimOpt<Real>;
 
  using size_type = typename PV<Real>::size_type;

private:

  const Ptr<Obj> obj_;
  const Ptr<V>   ui_;           // Initial condition
  
  VectorWorkspace<Real> workspace_; // Memory management

 
protected:

  PV& partition( V& x ) const { return static_cast<PV&>(x); }

  const PV& partition( const V& x ) const { return static_cast<const PV&>(x); }

public:

  Objective_SerialSimOpt( const Ptr<Obj>& obj, const V& ui ) : 
    obj_(obj), 
    ui_( workspace_.copy( ui ) ),
    u0_->zero(); z0_->zero();
  }

  
  virtual Real value( const V& u, const V& z, Real& tol ) override {

    auto& up = partition(u);
    auto& zp = partition(z);    

    // First step
    Real result = obj_->value( *ui_, *(up.get(0)), *(zp.get(0)), tol );
    
    for( size_type k=1; k<up.numVector(); ++k ) {
      result += obj_->value( *(up.get(k-1), *(up.get(k)), *(zp.get(k)), tol );    
    }
    return result;
  }

  virtual void gradient_1( V& g, const V &u, const V& z, Real& tol ) override {

    auto& up = partition(u);
    auto& zp = partition(z);    
    auto& gp = partition(g);

    // First step
    obj_->gradient_1( *(gp.get(0)), *ui_, *ui_, *(up.get(0)), *(zp.get(0)), tol );

    for( size_type k=1; k<up.numVector(); ++k ) {
      obj_->gradient_1( *(gp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), tol );
    }
  }

  virtual void gradient_2( V& g, const V& u, const V& z, Real& tol ) override {

    auto& up = partition(u);
    auto& zp = partition(z);    
    auto& gp = partition(g);

    // First step
    obj_->gradient_2( *(gp.get(0)), *ui_, *ui_, *(up.get(0)), *(zp.get(0)), tol );

    for( size_type k=1; k<up.numVector(); ++k ) {
      obj_->gradient_2( *(gp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), tol );
    }
  }

  virtual void hessVec_11( V& hv, const V& v, 
                           const V& u,  const V& z, Real& tol ) override {
  
    auto& hvp = partition(hv);
    auto& vp  = partition(v);
    auto& up  = partition(u);
    auto& zp  = partition(z);    
    auto& gp  = partition(g);

    // First step
    obj_->hessVec_11( *(hvp.get(0)), *(vp.get(0)), *ui_, *(up.get(0)), *(zp.get(0)), tol );

    for( size_type k=1; k<up.numVector(); ++k ) {
      obj_->hessVec_11( *(hvp.get(k)), *(vp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), tol );
    }
  }

  virtual void hessVec_12( V& hv, const V& v, const V&u, const V&z, Real &tol ) override {

    auto& hvp = partition(hv);
    auto& vp  = partition(v);
    auto& up  = partition(u);
    auto& zp  = partition(z);    
    auto& gp  = partition(g);

    // First step
    obj_->hessVec_12( *(hvp.get(0)), *(vp.get(0)), *ui_, *(up.get(0)), *(zp.get(0)), tol );

    for( size_type k=1; k<up.numVector(); ++k ) {
      obj_->hessVec_12( *(hvp.get(k)), *(vp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), tol );
    }
  }
 
  virtual void hessVec_21( V&hv, const V&v, const V&u, const V&z, Real &tol ) override {

    auto& hvp = partition(hv);
    auto& vp  = partition(v);
    auto& up  = partition(u);
    auto& zp  = partition(z);    
    auto& gp  = partition(g);

    // First step
    obj_->hessVec_21( *(hvp.get(0)), *(vp.get(0)), *ui_, *(up.get(0)), *(zp.get(0)), tol );

    for( size_type k=1; k<up.numVector(); ++k ) {
      obj_->hessVec_21( *(hvp.get(k)), *(vp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), tol );
    }
  }

  virtual void hessVec_22( V&hv, const V&v, const V&u,  const V&z, Real &tol ) override {

    auto& hvp = partition(hv);
    auto& vp  = partition(v);
    auto& up  = partition(u);
    auto& zp  = partition(z);    
    auto& gp  = partition(g);

    // First step
    obj_->hessVec_22( *(hvp.get(0)), *(vp.get(0)), *ui_, *(up.get(0)), *(zp.get(0)), tol );

    for( size_type k=1; k<up.numVector(); ++k ) {
      obj_->hessVec_22( *(hvp.get(k)), *(vp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), tol );
    }
  }

}; // class Objective_SerialSimOpt

} // namespace ROL


#endif // ROL_OBJECTIVE_SERIALSIMOPT_HPP

