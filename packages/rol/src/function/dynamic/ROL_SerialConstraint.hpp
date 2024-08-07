// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_SERIALCONSTRAINT_HPP
#define ROL_SERIALCONSTRAINT_HPP

#include <type_traits>

#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_DynamicConstraint.hpp"
#include "ROL_SerialFunction.hpp"

/** @ingroup func_group
    \class ROL::SerialConstraint
    \brief Evaluates ROL::DynamicConstraint over a sequential set of time intervals

    ---
*/

namespace ROL {

template<typename Real>
class SerialConstraint : public Constraint_SimOpt<Real>,
                         public SerialFunction<Real> {
private:

  using PV = PartitionedVector<Real>;  
  using SerialFunction<Real>::ts;
  using SerialFunction<Real>::clone;

  Ptr<DynamicConstraint<Real>>  con_; // Constraint over a single time step

public:

  using size_type = typename std::vector<Real>::size_type;
  using SerialFunction<Real>::numTimeSteps;
  using SerialFunction<Real>::getZeroState;
  using SerialFunction<Real>::getInitialCondition;
  using SerialFunction<Real>::getSkipInitialCondition;

  SerialConstraint( const Ptr<DynamicConstraint<Real>>& con, 
                    const Vector<Real>& u_initial, 
                    const TimeStampsPtr<Real>& timeStampsPtr ) : 
    SerialFunction<Real>::SerialFunction( u_initial, timeStampsPtr ), 
    con_(con) {}

  virtual void solve(       Vector<Real>& c, 
                            Vector<Real>& u,
                      const Vector<Real>& z, 
                      Real& tol ) override {

    auto& cp = partition(c);
    auto& up = partition(u);
    auto& zp = partition(z);

    if( !getSkipInitialCondition() )
      con_->solve( cp[0], getInitialCondition(), up[0], zp[0], ts(0) );

    for( size_type k=1; k<numTimeSteps(); ++k ) 
       con_->solve( cp[k], up[k-1], up[k], zp[k], ts(k) );

  } // solve

  virtual void update( const Vector<Real> &u, const Vector<Real> &z, bool flag = true, int iter = -1 ) override {
    auto& up = partition(u);
    auto& zp = partition(z);

    if( !getSkipInitialCondition() )
      con_->update( getInitialCondition(), up[0], zp[0], ts(0) );

    for( size_type k=1; k<numTimeSteps(); ++k )
      con_->update( up[k-1], up[k], zp[k], ts(k) );
  }

   
  using Constraint_SimOpt<Real>::value;
  virtual void value(       Vector<Real>& c, 
                      const Vector<Real>& u, 
                      const Vector<Real>& z, 
                      Real& tol ) override {

    auto& cp = partition(c);
    auto& up = partition(u);
    auto& zp = partition(z);
  
    if( !getSkipInitialCondition() )
      con_->value( cp[0], getInitialCondition(), up[0], zp[0], ts(0) );
    
    for( size_type k=1; k<numTimeSteps(); ++k ) 
      con_->value( cp[k], up[k-1], up[k], zp[k], ts(k) );

  } // value


  virtual void applyJacobian_1(       Vector<Real>& jv, 
                                const Vector<Real>& v, 
                                const Vector<Real>& u, 
                                const Vector<Real>& z, 
                                Real& tol ) override {

    auto& jvp = partition(jv);   auto& vp = partition(v);
    auto& up  = partition(u);    auto& zp = partition(z);

    auto  tmp = clone(jvp[0]);
    auto& x   = *tmp; 

    if( !getSkipInitialCondition() )
      con_->applyJacobian_un( jvp[0],  vp[0], getZeroState(), up[0], zp[0], ts(0) );
    
    for( size_type k=1; k<numTimeSteps(); ++k ) {

      con_->applyJacobian_uo( x, vp[k-1], up[k-1], up[k], zp[k], ts(k) );
      con_->applyJacobian_un( jvp[k], vp[k], up[k-1], up[k], zp[k], ts(k) );
      jvp.get(k)->plus(x);

    } // end for

  } // applyJacobian_1
  

  virtual void applyInverseJacobian_1(       Vector<Real>& ijv, 
                                       const Vector<Real>& v, 
                                       const Vector<Real>& u, 
                                       const Vector<Real>& z, 
                                       Real& tol) override {

    auto& ijvp = partition(ijv);  auto& vp = partition(v);
    auto& up   = partition(u);    auto& zp = partition(z);

    auto  tmp  = clone(ijvp[0]);
    auto& x    = *tmp;

    if( !getSkipInitialCondition() ) 
      con_->applyInverseJacobian_un( ijvp[0], vp[0], getZeroState(), 
                                     up[0], zp[0], ts(0) );

    for( size_type k=1; k<numTimeSteps(); ++k ) {

      con_->applyJacobian_uo( x, ijvp[k-1], up[k-1], up[k], zp[k], ts(k) );
      x.scale(-1.0);
      x.plus( vp[k] );
      con_->applyInverseJacobian_un( ijvp[k], x, up[k-1], up[k], zp[k], ts(k) );

    } // end for

  } // applyInverseJacobian_1
  
 
  using Constraint_SimOpt<Real>::applyAdjointJacobian_1;
  virtual void applyAdjointJacobian_1(       Vector<Real>& ajv, 
                                       const Vector<Real>& v, 
                                       const Vector<Real>& u, 
                                       const Vector<Real>& z, 
                                       const Vector<Real>& dualv, 
                                       Real& tol) override {

    auto& ajvp  = partition(ajv);   auto& vp = partition(v);
    auto& up    = partition(u);     auto& zp = partition(z);

    auto  tmp  = clone(ajvp[0]); 
    auto& x    = *tmp; 

    if( !getSkipInitialCondition() )
      con_->applyAdjointJacobian_un( ajvp[0],  vp[0], getZeroState(), up[0], zp[0], ts(0) );
   
    for( size_type k=1; k<numTimeSteps(); ++k ) {

      con_->applyAdjointJacobian_un( ajvp[k], vp[k], up[k-1], up[k], zp[k], ts(k) );
      con_->applyAdjointJacobian_uo( x, vp[k], up[k-1], up[k], zp[k], ts(k) );
      ajvp[k-1].plus(x);

    } // end for

  } // applyAdjointJacobian_1

  void applyInverseAdjointJacobian_1(       Vector<Real>& iajv, 
                                      const Vector<Real>& v, 
                                      const Vector<Real>& u, 
                                      const Vector<Real>& z, 
                                      Real& tol) override {

    auto& iajvp = partition(iajv);  auto& vp = partition(v);
    auto& up    = partition(u);     auto& zp = partition(z);

    auto  tmp  = clone(iajvp[0]);
    auto& x    = *tmp; 

    size_type k = numTimeSteps()-1;

    con_->applyInverseAdjointJacobian_un( iajvp[k], vp[k], up[k-1], up[k], zp[k], ts(k) );

    for( size_type k=numTimeSteps()-2; k>0; --k ) {

      con_->applyAdjointJacobian_uo( x, iajvp[k+1], up[k], up[k+1], zp[k+1], ts(k+1) );
      x.scale(-1.0);
      x.plus( vp[k] );
      con_->applyInverseAdjointJacobian_un( iajvp[k], x, up[k-1], up[k], zp[k], ts(k) );             

    } // end for

    con_->applyAdjointJacobian_uo( x, iajvp[1], up[0], up[1], zp[1], ts(1) );
    x.scale(-1.0);
    x.plus( vp[0] );

    if( !getSkipInitialCondition() ) 
      con_->applyInverseAdjointJacobian_un( iajvp[0], x, getZeroState(), up[0], zp[0], ts(0) );           

    // this weird condition places iajvp in the final vector slot    
    else iajvp[0].set(x);
     
  } // applyInverseAdjointJacobian_1


  virtual void applyJacobian_2(       Vector<Real>& jv, 
                                const Vector<Real>& v, 
                                const Vector<Real>& u, 
                                const Vector<Real>& z, 
                                Real &tol ) override { 
 
    auto& jvp = partition(jv);   auto& vp = partition(v);
    auto& up  = partition(u);    auto& zp = partition(z);
  
    if( !getSkipInitialCondition() )
      con_->applyJacobian_z( jvp[0], vp[0], getInitialCondition(), up[0], zp[0], ts(0) );
  
    for( size_type k=1; k<numTimeSteps(); ++k ) 
      con_->applyJacobian_z( jvp[k], vp[k], up[k-1], up[k], zp[k], ts(k) );

  } // applyJacobian_2

 
  using Constraint_SimOpt<Real>::applyAdjointJacobian_2;
  virtual void applyAdjointJacobian_2(       Vector<Real>& ajv, 
                                       const Vector<Real>& v, 
                                       const Vector<Real>& u, 
                                       const Vector<Real>& z, 
                                       Real& tol ) override {

    auto& ajvp  = partition(ajv);   auto& vp = partition(v);
    auto& up    = partition(u);     auto& zp = partition(z);

    if( !getSkipInitialCondition() )
      con_->applyAdjointJacobian_z( ajvp[0], vp[0], getInitialCondition(), up[0], zp[0], ts(0) );

    for( size_type k=1; k<numTimeSteps(); ++k ) 
      con_->applyAdjointJacobian_z( ajvp[k], vp[k], up[k-1], up[k], zp[k], ts(k) );

  } // applyAdjointJacobian_2



  virtual void applyAdjointHessian_11(       Vector<Real>& ahwv, 
                                       const Vector<Real>& w, 
                                       const Vector<Real>& v, 
                                       const Vector<Real>& u, 
                                       const Vector<Real>& z, 
                                       Real& tol) override {

    auto& ahwvp = partition(ahwv);    auto& wp = partition(w);
    auto& vp    = partition(v);       auto& up = partition(u);
    auto& zp    = partition(z);

    auto  tmp  = clone(ahwvp[0]); 
    auto& x    = *tmp; 

   if( !getSkipInitialCondition() ) {

     con_->applyAdjointHessian_un_un( ahwvp[0], wp[0], vp[0], getZeroState(), up[0], zp[0], ts(0) );
     con_->applyAdjointHessian_un_uo( x, wp[1], vp[1], up[0], up[1], zp[1], ts(1) );
     ahwvp[0].plus(x);
     con_->applyAdjointHessian_uo_uo( x, wp[1], vp[0], up[0], up[1], zp[1], ts(1) );
     ahwvp[0].plus(x);

   }

    for( size_type k=1; k<numTimeSteps(); ++k ) {

      con_->applyAdjointHessian_un_un( ahwvp[k], wp[k], vp[k], up[k-1], up[k], zp[k], ts(k) );
      con_->applyAdjointHessian_uo_un( x, wp[k], vp[k-1], up[k-1], up[k], zp[k], ts(k) );
      ahwvp[k].plus(x);

      if( k < numTimeSteps()-1 ) {
        con_->applyAdjointHessian_un_uo( x, wp[k+1], vp[k+1], up[k], up[k+1], zp[k+1], ts(k+1) );
        ahwvp[k].plus(x);
        con_->applyAdjointHessian_uo_uo( x, wp[k+1], vp[k], up[k], up[k+1], zp[k+1], ts(k+1) );
        ahwvp[k].plus(x);
      } // endif
    } // endfor

  } // applyAdjointHessian_11


  // TODO: Implement off-diagonal blocks  
  virtual void applyAdjointHessian_12(       Vector<Real>& ahwv,
                                       const Vector<Real>& w,
                                       const Vector<Real>& v,
                                       const Vector<Real>& u,
                                       const Vector<Real>& z,
                                       Real& tol) override {
  }


  virtual void applyAdjointHessian_21(       Vector<Real>& ahwv,
                                       const Vector<Real>& w,
                                       const Vector<Real>& v,
                                       const Vector<Real>& u,
                                       const Vector<Real>& z,
                                       Real& tol) override { 
  }

  virtual void applyAdjointHessian_22(       Vector<Real>& ahwv, 
                                       const Vector<Real>& w, 
                                       const Vector<Real>& v,
                                       const Vector<Real>& u, 
                                       const Vector<Real>& z, 
                                       Real& tol) override {

    auto& ahwvp = partition(ahwv); auto& vp = partition(v);   auto& wp = partition(w);
    auto& up    = partition(u);    auto& zp = partition(z);
  
    con_->applyAdjointHessian_z_z( ahwvp[0], wp[0], vp[0], getInitialCondition(), 
                                   up[0], zp[0], ts(0) );
  
    for( size_type k=1; k<numTimeSteps(); ++k ) {
      con_->applyAdjointHessian_z_z( ahwvp[k], wp[k], vp[k], up[k-1], 
                                     up[k], zp[k], ts(k) );
    }

  } // applyAdjointHessian_22


}; // SerialConstraint

// Helper function to create a new SerialConstraint

template<typename DynCon, typename Real, typename P = Ptr<SerialConstraint<Real>>>
inline typename std::enable_if<std::is_base_of<DynamicConstraint<Real>,DynCon>::value,P>::type
make_SerialConstraint( const Ptr<DynCon>& con,
                       const Vector<Real>& u_initial,
                       const TimeStampsPtr<Real> timeStampsPtr ) {
  return makePtr<SerialConstraint<Real>>(con,u_initial,timeStampsPtr);
}






} // namespace ROL


#endif // ROL_SERIALCONSTRAINT_HPP
