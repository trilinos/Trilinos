// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_CONSTRAINT_SERIALSIMOPT_HPP
#define ROL_CONSTRAINT_SERIALSIMOPT_HPP

namespace ROL {

/** @ingroup func_group
    \class ROL::Constraint_SerialSimOpt
    \brief Unifies the constraint defined on a single time step that are
           defined through the Constraint_TimeSimOpt interface into a 
           SimOpt constraint for all time. Primarily intended for use
           in testing the parallel-in-time implementation.

           In contrast to Constraint_PinTSimOpt, where the argument vectors
           are PartitionedVectors with Vector_SimOpt, the approach here is
           to consider the state and control to be of PartitionedVector type
           directly.


   NOTE:   This class does not address non-autonomous systems as constraints.
           The formulas involving adjoint Jacobians are almost certainly 
           incorrect in such cases. This will need to be fixed if the 
           Constraint_TimeSimOpt interface can accommodate explicit functions 
           of time. Additionally, opt vectors are assumed to be constant in 
           time on a time step.

*/

#include "ROL_Constraint_TimeSimOpt.hpp"
#include "ROL_PartitionedVector.hpp"
#include "ROL_VectorWorkspace.hpp"


template<typename Real> 
class Constraint_SerialSimOpt : public Constraint_SimOpt<Real> {

  using V      = Vector<Real>;
  using PV     = PartitionedVector<Real>;
  using Con = Constraint_TimeSimOpt<Real>;
 
  using size_type = typename PV<Real>::size_type;

private:

  const Ptr<Con> con_;          // Constraint for single time step
  const Ptr<V>   ui_;           // Initial condition
  Ptr<V> u0_;                   // Zero sim vector on time step
  Pre<V> z0_;                   // Zero opt vector on time step
 
  VectorWorkspace<Real> workspace_; // Memory management

protected:

  PV& partition( V& x ) const { return static_cast<PV&>(x); }

  const V& partition( const V& x ) const { return static_cast<const PV&>(x); }


public:

  Constraint_SerialSimOpt( const Ptr<Con>& con, const V& ui, const V& zi ) :
    Constraint_SimOpt<Real>(), con_(con), 
    ui_( workspace_.copy( ui ) ),
    u0_( workspace_.clone( ui ) ),
    z0_( workspace_.clone( ui ) ) {
    u0_->zero(); z0_->zero();
  }

//  virtual void update_1( const V& u, const V& z, bool flag = true, int iter = -1 ) override {
//    
//  }
//  
//  virtual void update_2( const V& u, const V& z, bool flag = true, int iter = -1 ) override {
//
//  }
  
  virtual void value( V& c, const V& u, const V& z, Real& tol ) override {
  
    auto& cp = partition(c);
    auto& up = partition(u);
    auto& zp = partition(z);
    
    // First interval value uses initial condition
    con_->value( *(cp.get(0)), *ui_, *(zp.get(0)), tol );

    for( size_type k=1, k<cp.numVectors(), ++k ) {
      con_->value( *(cp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), tol );
    }
  }
  
  virtual void solve( V& c, V& u, const V& z, Real& tol ) override {

    auto& cp = partition(c);
    auto& up = partition(u);
    auto& zp = partition(z);
    
    // First interval value uses initial condition
    con_->solve( *(cp.get(0)), *ui_, *(zp.get(0)), tol );

    for( size_type k=1, k<cp.numVectors(), ++k ) {
      con_->solve( *(cp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), tol );
    }
  }
  
  virtual void applyJacobian_1( V& jv, const V& v, const V& u, const V& z, Real& tol ) override {
    
    auto& jvp = partition(jv);      auto& vp = partition(v);
    auto& up  = partition(u);       auto& zp = partition(z);

    auto tmp = workspace_.clone( up.get(0) ); 

    con_->applyJacobian_1_new( *(jvp.get(0)), *(vp.get(0)), *u0_, *(up.get(0)), *(zp.get(0)), tol );
    
    for( size_type k=1; k<jvp.numVectors(); ++jvp ) {
      con_->applyJacobian_1_new( *(jvp.get(k)), *(vp.get(k)), *(up.get(k)), *(up.get(k-1)), *(zp.get(k)), tol );
      con_->applyJacobian_1_old( *tmp, *(vp.get(k)), *(up.get(k)), *(up.get(k-1)), *(zp.get(k)), tol );
      jvp.plus(tmp);
    }
  } // applyJacobian_1
  
  virtual void applyJacobian_2( V& jv, const V& v, const V& u, const V& z, Real& tol)  override { 

    auto& jvp = partition(jv);      auto& vp = partition(v);
    auto& up  = partition(u);       auto& zp = partition(z);

    for( size_type k=0; k<jvp.numVectors(); ++k ) {
      con_->applyJacobian_2( *(jvp.get(k)), *(vp.get(k)), *(up.get(k)), *(up.get(k-1)), *(zp.get(k)), tol );
    }
  } // applyJacobian_2
  

  virtual void applyInverseJacobian_1( V& ijv, const V& v, const V& u, 
                                       const V& z, Real& tol ) override {

    auto& ijvp = partition(ijv);     auto& vp = partition(v);
    auto& up   = partition(u);       auto& zp = partition(z);

    // First step 
    con_->applyInverseJacobian_1_new( *(ijvp.get(0)), *(vp.get(0), *u0_, *(up.get(0)), *(zp.get(0)), tol );
    
    auto tmp = workspace.clone( vp.get(0) );

    for( size_type k=1; k<ijvp.numVectors(); ++k ) { 

      con_->applyJacobian_1_old( *tmp, *(ijvp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), tol );
      tmp->scale(-1.0);
      tmp->plus( *(vp.get(k) );
   
      con_->applyInverseJacobian_1_new( *(ijvp.get(k)), *tmp, *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), tol );
    }
  } // applyInverseJacobian_1


  
  virtual void applyAdjointJacobian_1( V& ajv, const V& v, const V& u, 
                                       const V& z, Real& tol) override {
 
    auto& ajvp = partition(ajv);     auto& vp = partition(v);
    auto& up   = partition(u);       auto& zp = partition(z);

    auto tmp = workspace.clone( ajvp.get(0) );

    size_type N = ajvp.numVectors()-1;

    // First step
    con_->applyAdjointJacobian_1_new( *(ajvp.get(0)), *(vp.get(0)), *u0_, *(up.get(0)), *(zp.get(0)) );
    con_->applyAdjointJacobian_1_old( *(ajvp.get(0)), *(vp.get(1)), *u0_, *(up.get(0)), *(zp.get(0)) );

    for( size_type k=1; k<N; ++k ) {
      con_->applyAdjointJacobian_1_new( *(ajvp.get(k)), *(vp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)) );
      con_->applyAdjointJacobian_1_old( *tmp, *(vp.get(k+1)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)) );
      ajvp.plus(*tmp);
    }

    // Last step
    con_->applyAdjointJacobian_1_new( *(ajvp.get(k)), *(vp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)) );

  } // applyAdjointJacobian_1
  


  virtual void applyAdjointJacobian_1( V& ajv, const V& v, const V& u, const V &z, 
                                       const V& dualv, Real& tol) override {

    auto& ajvp = partition(ajv);     auto& vp = partition(dualv);
    auto& up   = partition(u);       auto& zp = partition(z);
    auto tmp = workspace.clone( ajvp.get(0) );

    size_type N = ajvp.numVectors()-1;

    // First step
    con_->applyAdjointJacobian_1_new( *(ajvp.get(0)), *(vp.get(0)), *u0_, *(up.get(0)), *(zp.get(0)) );
    con_->applyAdjointJacobian_1_old( *(ajvp.get(0)), *(vp.get(1)), *u0_, *(up.get(0)), *(zp.get(0)) );

    for( size_type k=1; k<N; ++k ) {
      con_->applyAdjointJacobian_1_new( *(ajvp.get(k)), *(vp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)) );
      con_->applyAdjointJacobian_1_old( *tmp, *(vp.get(k+1)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)) );
      ajvp.plus(*tmp);
    }

    // Last step
    con_->applyAdjointJacobian_1_new( *(ajvp.get(k)), *(vp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)) );
 
  } // applyAdjointJacobian_1


  
  virtual void applyAdjointJacobian_2( V& ajv, const V& v, const V& u, 
                                       const V& z, Real& tol) override {
    auto& ajvp = partition(ajv);    auto& vp = partition(v);
    auto& up   = partition(u);      auto& zp = partition(z);

    for( size_type k=0; k<jvp.numVectors(); ++k ) {
      con_->applyAdjointJacobian_2_time( *(jvp.get(k)), *(vp.get(k)), *(up.get(k)), *(up.get(k-1)), *(zp.get(k)), tol );
    }
  } // applyAdjointJacobian_2


  virtual void applyAdjointJacobian_2( V& ajv, const V& v, const V& u, const V& z, 
                                       const V& dualv, Real& tol) override {

    auto& ajvp = partition(ajv);    auto& vp = partition(v);
    auto& up   = partition(u);      auto& zp = partition(z);

    for( size_type k=0; k<jvp.numVectors(); ++k ) {
      con_->applyAdjointJacobian_2_time( *(jvp.get(k)), *(vp.get(k)), *(up.get(k)), *(up.get(k-1)), *(zp.get(k)), tol );
    }
  }


  virtual void applyInverseAdjointJacobian_1( V& iajv, const V& v, const V& u, 
                                              const V& z, Real& tol) override {

    auto& iajvp = partition(iajv);   auto& vp = partition(v);
    auto& up    = partition(u);      auto& zp = partition(z);

    size_type N = iajvp.numVectors()-1;

    // Last Step
    con_->applyInverseAdjointJacobian_1_new( *(iajvp.get(N)), *(vp.get(N)), *(up.get(N-1)), *(up.get(N)), *(zp.get(N)), tol );
    
    auto tmp = workspace_.clone( vp.get(0) );

    for( size_type k=N-1; k>0; --k ) {
      con_->applyAdjointJacobian_1_old( *tmp, *(iajvp.get(k+1)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), tol );

      tmp->scale( -1.0 );
      tmp->plus( *(vp.get(k) );

      con_->applyAdjointJacobian_1_new( *(iajvp.get(k)), *tmp, *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), tol );
    }

    // "First step"
    con_->applyAdjointJacobian_1_old( *tmp, *(iajvp.get(1)), *u0_, *(up.get(0)), *(zp.get(0)), tol );

    tmp->scale( -1.0 );
    tmp->plus( *(vp.get(0) );

    con_->applyAdjointJacobian_1_new( *(iajvp.get(0)), *tmp, *u0_, *(up.get(0)), *(zp.get(0)), tol );
  }


  virtual void applyAdjointHessian_11( V& ahwv, const V& w, const V& v, 
                                       const V& u, const V& z, Real& tol) override { 
     // TODO 
  }


  virtual void applyAdjointHessian_11( V& ahwv, const V& w, const V& v, 
                                       const V& u, const V& z, Real& tol) override {
    ahwv.zero();
  }

  virtual void applyAdjointHessian_11( V& ahwv, const V& w, const V& v, 
                                       const V& u, const V& z, Real& tol) override { 
    ahwv.zero();
  }

  virtual void applyAdjointHessian_11( V& ahwv, const V& w, const V& v, 
                                       const V& u, const V& z, Real& tol) override {
    ahwv.zero();
  }


 

   

  
  


}; // class Constraint_SerialSimOpt



} // namespace ROL


#endif // ROL_CONSTRAINT_SERIALSIMOPT_HPP

